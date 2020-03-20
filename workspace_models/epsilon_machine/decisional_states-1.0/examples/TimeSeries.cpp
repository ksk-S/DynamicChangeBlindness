/*
This file is part of the decisional state reconstruction algorithm
technique exposed in "Decisional States", by Nicolas Brodu.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free
    Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
    MA  02110-1301  USA

See http://nicolas.brodu.numerimoire.net/en/programmation/decisional_states/index.html
for more information and possibly updates.

Copyright holder: Nicolas Brodu <nicolas.brodu@numerimoire.net>
File release Date: February 09
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>

#include <boost/array.hpp>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/unordered_set.hpp>
#include <boost/lexical_cast.hpp>

#include "decisional_states.hpp"
#include "decisional_states_joint_data_manager.hpp"
#include "decisional_states_optimisers.hpp"

#include <stdlib.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace boost;
using namespace decisional_states;

#include <stdlib.h>
#include <signal.h>
#include <libgen.h>

typedef float FloatType;

#ifndef N_PAST_SAMPLES
#define N_PAST_SAMPLES 5
#endif

#ifndef N_FUTURE_SAMPLES
#define N_FUTURE_SAMPLES 1
#endif

#ifndef TIMELAG
#define TIMELAG 1
#endif

static const int n_past_samples = N_PAST_SAMPLES;
static const int n_future_samples = N_FUTURE_SAMPLES;

// Taken's style time-lag reconstruction:
// past = [ x(t-timelag*(M-1)), ..., x(t-timelag), x(t)]  with M = n_past_samples
// future = [ x(t+timelag), x(t+timelag*2), ..., x(t+timelag*M)] with M = n_future_samples
// In other words, take all possible subsamplings dividing by timelag in the series
// If n_past_samples = K * n_future_samples, we have a series of the variable
// that would be produced by Taken's theorem. Otherwise we "just" have subsamplings.
// timelag = 1 is the original series.
// From a signal processing POV, subsampling introduces spurious frequencies if a
// filter is not applied before to remove all freqs above the Nyquist freq of the result series.
// How that interferes with Taken's reconstruction would be an interesting topic... mail me if
// you wish to investigate that, I'll be interested :)
static const int timelag = TIMELAG;

// number of data per line = dimension of one vector per line in the file
// the series length (in time) is the number of lines in the file
static const int base_sample_dim = 1;

static int spacialVariance = 10000;
//static int numSpacialSamples = 1000;

template <int NPAST, int NFUTURE>
struct JointData {
    typedef array<FloatType, NPAST*base_sample_dim> DataType;
    typedef array<FloatType, NFUTURE*base_sample_dim>  PredictionType;
    static const int static_size = (NPAST + NFUTURE) * base_sample_dim;
    DataType data;
    PredictionType prediction;
    static const int baselag = timelag * base_sample_dim;

    /// Joint concept API
    const DataType& getData() const {
        return data;
    }
    const PredictionType& getPrediction() const {
        return prediction;
    }
    void setData(const DataType& data) {
        this->data = data;
    }
    void setPrediction(const PredictionType& prediction) {
        this->prediction = prediction;
    }
    // optional function, used if specified to handle NaNs
    static bool isMissing(const DataType& d) {
        for (int i=0; i<NPAST; ++i) if (!isfinite(d[i])) return true;
        return false;
    }
    static bool isMissing(const PredictionType& p) {
        for (int i=0; i<NFUTURE; ++i) if (!isfinite(p[i])) return true;
        return false;
    }

    /// Extra utility for this time series
    /// valid from sampleidx = NPAST-1 to series.size()-static_size+NPAST-1= series.size()-1 - NFUTURE
    static JointData at(int sampleidx, const std::vector<FloatType>& series) {
        JointData ret;
        for (int i=0; i<NPAST; ++i) for (int j=0; j<base_sample_dim; ++j) ret.data[i*base_sample_dim + j] = series[sampleidx - (NPAST-1-i) * baselag + j];
        for (int i=0; i<NFUTURE; ++i) for (int j=0; j<base_sample_dim; ++j) ret.prediction[i*base_sample_dim + j] = series[sampleidx + (1+i) * baselag + j];
        return ret;
    }

    /// first index in the time series on which the "at" method can be called
    static int begin(const std::vector<FloatType>& series) {
        // from above: sampleidx - (NPAST-1-i) * baselag + j = 0,  with i=0 and j=0
        return (NPAST-1) * baselag;
    }

    /// just after the last index in the time series on which the "at" method can be called
    static int end(const std::vector<FloatType>& series) {
        // from above: sampleidx + (1+i) * baselag + j = series.size(),  with i=NFUTURE-1 and j=base_sample_dim-1
        return series.size() - ((NFUTURE) * baselag + base_sample_dim-1);
    }

    inline int size() const {return static_size;}
    inline const FloatType& operator[](int idx) const {return (idx<NPAST*base_sample_dim) ? data[idx] : prediction[idx-NPAST*base_sample_dim];}
    inline FloatType& operator[](int idx) {return (idx<NPAST*base_sample_dim) ? data[idx] : prediction[idx-NPAST*base_sample_dim];}
    bool operator==(const JointData& other) const {
        return (data == other.data) && (prediction == other.prediction);
    }

};
namespace boost {
template<int NPAST, int NFUTURE>
size_t hash_value(const JointData<NPAST,NFUTURE>& a) {
    size_t seed = 0;
    for(int i=0; i<JointData<NPAST,NFUTURE>::static_size; ++i) boost::hash_combine(seed, a[i]);
    return seed;
}
}

template<typename JointType, typename FloatType>
struct SeriesJointSampler {
    // sampler in series index space
    SimpleGaussianKernel<FloatType> spacialKernel;

    template<class RNG>
    SeriesJointSampler(int sampleidx, const std::vector<FloatType>& series, int _nsamples, RNG& rng) : spacialKernel(spacialVariance), nsamples(_nsamples), samples(new JointType[_nsamples]) {
        for (int i=0; i<nsamples; ++i) {
            int idx = -1;
            do {
                FloatType s = -1;
                spacialKernel.sample((FloatType)sampleidx+0.5f, s, rng);
                idx = (int)floorf(s);
            } while (idx<JointType::begin(series) || idx>=JointType::end(series));
            samples[i] = JointType::at(idx,series);
        }
    }
    int nsamples;
    boost::shared_array<JointType> samples;

    typedef JointType* iterator;
    typedef const JointType* const_iterator;
    iterator begin() {return &samples[0];}
    iterator end() {return &samples[nsamples];}
    const_iterator begin() const {return &samples[0];}
    const_iterator end() const {return &samples[nsamples];}
    JointType& operator[](unsigned int i) {
        return samples[i];
    }
    const JointType& operator[](unsigned int i) const {
        return samples[i];
    }
    unsigned int size() const {return nsamples;}
};


namespace boost {
template<std::size_t N>
size_t hash_value(const array<FloatType, N>& a) {
    size_t seed = 0;
    for(std::size_t i=0; i<N; ++i) boost::hash_combine(seed, a[i]);
    return seed;
}
}

template<typename T>
struct NeighbourHoodHash : public std::unary_function<T, std::size_t> {
    std::size_t operator()(T const& a) const {
        std::size_t seed = 0;
        for(int i=0; i<T::static_size; ++i) boost::hash_combine(seed, a[i]);
        return seed;
    }
};

struct TypeTraits {
    typedef ::FloatType FloatType;
    typedef JointData<n_past_samples,n_future_samples> JointType;
    typedef JointType::DataType DataType;
    typedef JointType::PredictionType PredictionType;

    typedef GaussianMixtureSampler<PredictionType,FloatType> PredictionSampler;

    typedef SimpleGaussianKernel<FloatType> JointKernel;

    typedef NearTreeNeighbourhoodFinder<DataType, FloatType> DataNeighbourhoodFinder;
    typedef NearTreeNeighbourhoodFinder<PredictionType, FloatType> PredictionNeighbourhoodFinder;
    typedef NearTreeNeighbourhoodFinder<JointType, FloatType> JointNeighbourhoodFinder;

    typedef IndexedDistributionStorage<FloatType> DistributionStorage;
};

typedef TypeTraits::JointType JointType;
typedef TypeTraits::DataType DataType;
typedef TypeTraits::PredictionType PredictionType;
typedef TypeTraits::PredictionSampler PredictionSampler;
typedef TypeTraits::JointKernel JointKernel;
typedef TypeTraits::JointNeighbourhoodFinder JointNeighbourhoodFinder;
typedef TypeTraits::DataNeighbourhoodFinder DataNeighbourhoodFinder;
typedef TypeTraits::PredictionNeighbourhoodFinder PredictionNeighbourhoodFinder;

typedef ExhaustiveOptimiser<PredictionSampler> Optimiser;
typedef JointDataManager<TypeTraits>::DataSet DataSet;
typedef JointDataManager<TypeTraits>::DistributionManager DistributionManager;
typedef DecisionalStatesAnalyser<DistributionManager, DataSet, Optimiser> DSA;

typedef DistributionManager::Distribution Distribution;

template<int dim> struct Utility {
    FloatType operator()(const PredictionType& p1, const PredictionType& p2) {
        FloatType negssq = 0;
        for (int i=0; i<PredictionType::static_size; ++i) negssq -= (p1[i] - p2[i])*(p1[i] - p2[i]);
        return negssq;
    }
};

template<> struct Utility<1> {
    FloatType operator()(const PredictionType& p1, const PredictionType& p2) {
        return -(p1[0] - p2[0])*(p1[0] - p2[0]);
    }
};

template<class Vector>
void vector_add(Vector& result, const Vector& v1, const Vector& v2) {
    for (int i = 0; i<(int)result.size(); ++i) result[i] = v1[i] + v2[i];
}

template<class Vector, typename Scalar>
void vector_scale(Vector& result, const Scalar& s, const Vector& v) {
    for (int i = 0; i<(int)result.size(); ++i) result[i] = v[i] * s;
}

template<class State>
struct StateSorter {
    bool operator()(typename boost::shared_ptr<State> a, typename boost::shared_ptr<State> b) {
        // less count = more complex
        // equal counts but more members
        return ((a->count) < (b->count)) || (((a->count) == (b->count)) && ((a->members->size()) > (b->members->size())));
    }
};

bool readDouble(istream& is, double& value) {
    string x;
    if (!(is >> x)) return false;
    try {
        value = boost::lexical_cast<double>(x);
    } catch(bad_lexical_cast &) {
        value = numeric_limits<double>::quiet_NaN();
    }
    return true;
}

int main(int argc, char**argv) {
    cout << "Using " << n_past_samples << " past samples, " << n_future_samples << " future, " << timelag << " subsampling factor." << endl;
    if (argc<5) {
        cout << "Arguments expected:  kernel_width_or_0  flag  output_dir  input_file1  [input_file2 [...]]" << endl;
        cout << "Negative kernel width values are scaling factors for the default adaptative kernel width." << endl;
        cout << "Flags (maskable): 0: no preprocessing. 1: take first differences. 2: compute the average kernel width over all input files, and use it as factor for the first argument. 4: stop after computing the kwidth." << endl;
        return 1;
    }
    
    
    int argi = 0;
    double ksizeparam = atof(argv[++argi]);
    int flag = atoi(argv[++argi]);

    string output_dir = argv[++argi];

    string filename = output_dir+"/statC.txt";
    ofstream statCfile(filename.replace(filename.end()-3, filename.end(), "txt").c_str());

    filename = output_dir+"/statD.txt";
    ofstream statDfile(filename.replace(filename.end()-3, filename.end(), "txt").c_str());

    filename = output_dir+"/statP.txt";
    ofstream statPfile(filename.replace(filename.end()-3, filename.end(), "txt").c_str());

    filename = output_dir+"/statU.txt";
    ofstream statUfile(filename.replace(filename.end()-3, filename.end(), "txt").c_str());

    timeval tv; gettimeofday(&tv,0);
    double totalTime = tv.tv_sec + tv.tv_usec * 1e-6;

    boost::mt19937 rng;
    //rng.seed(time(0));
    rng.seed(42);

    double avgkw = 0; // flag mode = 2
    
    const int nsamples = 1200;
    PredictionType nullPrediction;
    for (int i=0; i<PredictionType::static_size; ++i) nullPrediction[i] = 0;

    int argifiles = ++argi;

    // Auto-kernel estimation on all files
    if ((flag&2)!=0) {
        cout << "Computing average kernel width on all files: " << flush;
        for (argi = argifiles; argi<argc; ++argi) {
            string infilename = argv[argi];
            ifstream seriesfile(argv[argi]);
            if (!seriesfile) {
                cerr << "cannot read file: " << infilename << endl;
                return 1;
            }
            vector<FloatType> series;
            double number;
            array<double, base_sample_dim> prev;
            // take first differences to remove the global trend if requested
            if ((flag&1)!=0) for (int i=0; i<base_sample_dim; ++i) if (!readDouble(seriesfile,prev[i])) {
                cerr << "file does not contain even one sample: " << infilename << endl;
                continue;
            }
            int idx_sample = 0;
            while(readDouble(seriesfile,number)) {
                double value = number;
                if ((flag&1)!=0) value -= prev[idx_sample];
                series.push_back(value);
                if ((flag&1)!=0) {
                    prev[idx_sample] = number;
                    ++idx_sample; if (idx_sample==base_sample_dim) idx_sample=0;
                }
            }
            if (series.empty()) {
                cerr << "file results in an empty series: " << infilename << endl;
                continue;
            }
            vector<JointType> allJointData(JointType::end(series) - JointType::begin(series));
            helpers::ContainerMemberAdaptor<std::vector<JointType>, PredictionType, &JointType::prediction> allPredictions(allJointData);
            for (int i = JointType::begin(series); i<JointType::end(series); ++i) {
                allJointData[i-JointType::begin(series)] = JointType::at(i,series);
            }
            SimpleGaussianKernel<FloatType> jkernel(ksizeparam,JointType::static_size,1e-9f);
            JointNeighbourhoodFinder jnf;
            avgkw += jkernel.setSizeFromSamples(allJointData, jnf);
        }
        avgkw /= (argc - argifiles);
        cout << avgkw << endl;
        if ((flag&4)!=0) return 0;
        if (ksizeparam<0) ksizeparam = -ksizeparam;
        if (ksizeparam==0) ksizeparam = 1;
        ksizeparam = avgkw * ksizeparam;
    }

    for (argi = argifiles; argi<argc; ++argi) {
        double ksizej = ksizeparam;

        string infilename = argv[argi];

        cout << "Processing file " << infilename <<"..." << endl;

        ifstream seriesfile(argv[argi]);
        if (!seriesfile) {
            cerr << "cannot read file: " << infilename << endl;
            return 1;
        }


        // rudimentary file reading routine
        vector<FloatType> series;
        double number;
        array<double, base_sample_dim> prev;
        // take first differences to remove the global trend if requested
        if ((flag&1)!=0) for (int i=0; i<base_sample_dim; ++i) if (!readDouble(seriesfile,prev[i])) {
            cerr << "file does not contain even one sample: " << infilename << endl;
            continue;
        }
        int idx_sample = 0;
        while(readDouble(seriesfile,number)) {
            double value = number;
            if ((flag&1)!=0) value -= prev[idx_sample];
            series.push_back(value);
            if ((flag&1)!=0) {
                prev[idx_sample] = number;
                ++idx_sample; if (idx_sample==base_sample_dim) idx_sample=0;
            }
        }
        if (series.empty()) {
            cerr << "file results in an empty series: " << infilename << endl;
            continue;
        }

        // prepropressing to determine a reasonable kernel size value
        vector<JointType> allJointData(JointType::end(series) - JointType::begin(series));
        helpers::ContainerMemberAdaptor<std::vector<JointType>, PredictionType, &JointType::prediction> allPredictions(allJointData);
        for (int i = JointType::begin(series); i<JointType::end(series); ++i) {
            allJointData[i-JointType::begin(series)] = JointType::at(i,series);
        }
        SimpleGaussianKernel<FloatType> jkernel(ksizej,JointType::static_size,1e-9f);
        if (ksizej<=0) {
            cout << "computing the kernel size in joint data space..." << endl;
            JointNeighbourhoodFinder jnf;
            double kauto = jkernel.setSizeFromSamples(allJointData, jnf);
            if (ksizej<0) jkernel.setSize(kauto = -ksizej * kauto); // scaling factor
            ksizej= kauto;
        }
        cout << "Using a joint kernel size = " << ksizej << endl;

        // product of gaussians with var_dim variance
        SimpleGaussianKernel<FloatType> pkernel(ksizej,PredictionType::static_size,1e-9f);

        PredictionSampler psampler(
            nsamples,
            allPredictions,
            pkernel,
            rng,
            vector_add,
            vector_scale,
            nullPrediction
            //,1000 // EM steps - 0 is default = just data sampling
        );
        cout << "Prediction sampler is built" << endl;

        // Now, prepare a data set using this sampler for conditional distributions p(pred|data)
        // The sampler is fixed for the whole data set, so distributions are comparable
        JointDataManager<TypeTraits> jdm(psampler);

        SimpleGaussianKernel<FloatType> dkernel(ksizej,DataType::static_size,1e-9f);
        jdm.addSeparable(allJointData, dkernel, pkernel, allJointData);

        // done building the distributions
        cout << "Distributions are built" << endl;

        Optimiser optimiser(psampler);
        DataSet& dataset = jdm.asDataSet();
        DistributionManager& distributionManager = jdm.asDistributionManager();
        DSA analyser(distributionManager, dataset, optimiser);

        gettimeofday(&tv,0);
        double theTime = tv.tv_sec + tv.tv_usec * 1e-6;

        cout << "Computing causal states... " << flush;
        //analyser.computeCausalStates( DistributionManager::DistributionMatcher(0.95), AggregateClustering<RandomProvider>(3, RandomProvider(42)) );
        analyser.computeCausalStates();
        cout << "Num causal states: " << analyser.causalStates.size();
        analyser.buildCausalStateGraph();    
        cout << ". Number of recurrent states: " << analyser.causalStatesCounts.nrecurrent << endl;
        
        gettimeofday(&tv,0);
        theTime = tv.tv_sec + tv.tv_usec * 1e-6 - theTime;
        cout << "time spent for computing causal states: " << theTime << endl;

        gettimeofday(&tv,0);
        theTime = tv.tv_sec + tv.tv_usec * 1e-6;

        cout << "Applying utility..." << flush;
        analyser.applyUtility( Utility<PredictionType::static_size>() );
        cout << " done!" << endl;

        gettimeofday(&tv,0);
        theTime = tv.tv_sec + tv.tv_usec * 1e-6 - theTime;
        cout << "time spent for applying utility: " << theTime << endl;

        analyser.computeIsoUtilityStates();
        cout << "Num iso utility states: " << analyser.isoUtilityStates.size();
        analyser.buildIsoUtilityStateGraph();
        cout << ". Number of recurrent states: " << analyser.isoUtilityStatesCounts.nrecurrent << endl;
        
        analyser.computeIsoPredictionStates();
        cout << "Num iso prediction states: " << analyser.isoPredictionStates.size();
        analyser.buildIsoPredictionStateGraph();
        cout << ". Number of recurrent states: " << analyser.isoPredictionStatesCounts.nrecurrent << endl;

        analyser.computeDecisionalStates();
        cout << "Num decisional states: " << analyser.decisionalStates.size();
        analyser.buildDecisionalStateGraph();
        cout << ". Number of recurrent states: " << analyser.decisionalStatesCounts.nrecurrent << endl;

        double C = analyser.statisticalComplexity();
        cout << "Statistical complexity of series "<<infilename<<": " << C << endl;
        statCfile << C << endl;

        double U = analyser.isoUtilityComplexity();
        cout << "Iso-Utility complexity of series "<<infilename<<": " << U << endl;
        statUfile << U << endl;

        double P = analyser.isoPredictionComplexity();
        cout << "Iso-Prediction complexity of series "<<infilename<<": " << P << endl;
        statPfile << P << endl;

        double D = analyser.decisionalComplexity();
        cout << "Decisional complexity of series "<<infilename<<": " << D << endl;
        statDfile << D << endl;

        char* tmpcopy = strdup(infilename.c_str());
        string inbasename = basename(tmpcopy);
        free(tmpcopy);
        
        ofstream dot;
        
        filename = output_dir+"/cgraph_"+inbasename;
        dot.open(filename.replace(filename.end()-3, filename.end(), "dot").c_str());
        analyser.writeCausalStateGraph(dot);
        dot.close();
        filename = output_dir+"/cgraph_transients_"+inbasename;
        dot.open(filename.replace(filename.end()-3, filename.end(), "dot").c_str());
        analyser.writeCausalStateGraph(dot,false);
        dot.close();

        filename = output_dir+"/ugraph_"+inbasename;
        dot.open(filename.replace(filename.end()-3, filename.end(), "dot").c_str());
        analyser.writeIsoUtilityStateGraph(dot);
        dot.close();
        filename = output_dir+"/ugraph_transients_"+inbasename;
        dot.open(filename.replace(filename.end()-3, filename.end(), "dot").c_str());
        analyser.writeIsoUtilityStateGraph(dot,false);
        dot.close();

        filename = output_dir+"/pgraph_"+inbasename;
        dot.open(filename.replace(filename.end()-3, filename.end(), "dot").c_str());
        analyser.writeIsoPredictionStateGraph(dot);
        dot.close();
        filename = output_dir+"/pgraph_transients_"+inbasename;
        dot.open(filename.replace(filename.end()-3, filename.end(), "dot").c_str());
        analyser.writeIsoPredictionStateGraph(dot,false);
        dot.close();
        
        filename = output_dir+"/dgraph_"+inbasename;
        dot.open(filename.replace(filename.end()-3, filename.end(), "dot").c_str());
        analyser.writeDecisionalStateGraph(dot);
        dot.close();
        filename = output_dir+"/dgraph_transients_"+inbasename;
        dot.open(filename.replace(filename.end()-3, filename.end(), "dot").c_str());
        analyser.writeDecisionalStateGraph(dot,false);
        dot.close();
        
        // Possible: for 0 to <JointType::begin(series) fill output with 0
        
        filename = output_dir+"/localC_"+inbasename;
        ofstream outfile(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for(DataSet::iterator datait = dataset.begin(); datait != dataset.end(); ++datait) {
            FloatType lc = analyser.localStatisticalComplexity(dataset.data(datait));
            outfile << lc << endl;
        }
        outfile.close();

        filename = output_dir+"/localU_"+inbasename;
        outfile.open(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for(DataSet::iterator datait = dataset.begin(); datait != dataset.end(); ++datait) {
            FloatType lu = analyser.localIsoUtilityComplexity(dataset.data(datait));
            outfile << lu << endl;
        }
        outfile.close();
        
        filename = output_dir+"/localP_"+inbasename;
        outfile.open(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for(DataSet::iterator datait = dataset.begin(); datait != dataset.end(); ++datait) {
            FloatType lp = analyser.localIsoPredictionComplexity(dataset.data(datait));
            outfile << lp << endl;
        }
        outfile.close();
        
        filename = output_dir+"/localD_"+inbasename;
        outfile.open(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for(DataSet::iterator datait = dataset.begin(); datait != dataset.end(); ++datait) {
            FloatType ld = analyser.localDecisionalComplexity(dataset.data(datait));
            outfile << ld << endl;
        }
        outfile.close();

        // compute state members while outputing series
        map<int, set<DataType> > cmembers, umembers, pmembers, dmembers;
        
        filename = output_dir+"/causal_states_series_"+inbasename;
        outfile.open(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for(DataSet::iterator datait = dataset.begin(); datait != dataset.end(); ++datait) {
            DataType data = dataset.data(datait);
            int index = analyser.getCausalState(data)->index;
            cmembers[index].insert(data);
            outfile << index << endl;
        }
        outfile.close();
        
        filename = output_dir+"/iso_utility_states_series_"+inbasename;
        outfile.open(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for(DataSet::iterator datait = dataset.begin(); datait != dataset.end(); ++datait) {
            DataType data = dataset.data(datait);
            int index = analyser.getIsoUtilityState(data)->index;
            umembers[index].insert(data);
            outfile << index << endl;
        }
        outfile.close();
        
        filename = output_dir+"/iso_prediction_states_series_"+inbasename;
        outfile.open(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for(DataSet::iterator datait = dataset.begin(); datait != dataset.end(); ++datait) {
            DataType data = dataset.data(datait);
            int index = analyser.getIsoPredictionState(data)->index;
            pmembers[index].insert(data);
            outfile << index << endl;
        }
        outfile.close();
        
        filename = output_dir+"/decisional_states_series_"+inbasename;
        outfile.open(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for(DataSet::iterator datait = dataset.begin(); datait != dataset.end(); ++datait) {
            DataType data = dataset.data(datait);
            int index = analyser.getDecisionalState(data)->index;
            dmembers[index].insert(data);
            outfile << index << endl;
        }
        outfile.close();

        // now output state members
        // members of each state.
        filename = output_dir+"/causal_states_members_"+inbasename;
        outfile.open(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for (map<int, set<DataType> >::iterator it=cmembers.begin(); it!=cmembers.end(); ++it) {
            if (it!=cmembers.begin()) outfile << endl;
            outfile << "State " << it->first << ":" << endl;
            for (set<DataType>::iterator sit = it->second.begin(); sit != it->second.end(); ++sit) {
                for (unsigned int i=0; i<sit->size(); ++i) {
                    if (i>0) outfile << " ";
                    outfile << (*sit)[i];
                }
                outfile << endl;
            }
        }
        outfile.close();
        
        filename = output_dir+"/iso_utility_states_members_"+inbasename;
        outfile.open(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for (map<int, set<DataType> >::iterator it=umembers.begin(); it!=umembers.end(); ++it) {
            if (it!=umembers.begin()) outfile << endl;
            outfile << "State " << it->first << ":" << endl;
            for (set<DataType>::iterator sit = it->second.begin(); sit != it->second.end(); ++sit) {
                for (unsigned int i=0; i<sit->size(); ++i) {
                    if (i>0) outfile << " ";
                    outfile << (*sit)[i];
                }
                outfile << endl;
            }
        }
        outfile.close();
                
        filename = output_dir+"/iso_prediction_states_members_"+inbasename;
        outfile.open(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for (map<int, set<DataType> >::iterator it=pmembers.begin(); it!=pmembers.end(); ++it) {
            if (it!=pmembers.begin()) outfile << endl;
            outfile << "State " << it->first << ":" << endl;
            for (set<DataType>::iterator sit = it->second.begin(); sit != it->second.end(); ++sit) {
                for (unsigned int i=0; i<sit->size(); ++i) {
                    if (i>0) outfile << " ";
                    outfile << (*sit)[i];
                }
                outfile << endl;
            }
        }
        outfile.close();
                
        filename = output_dir+"/decisional_states_members_"+inbasename;
        outfile.open(filename.replace(filename.end()-3, filename.end(), "txt").c_str());
        for (map<int, set<DataType> >::iterator it=dmembers.begin(); it!=dmembers.end(); ++it) {
            if (it!=dmembers.begin()) outfile << endl;
            outfile << "State " << it->first << ":" << endl;
            for (set<DataType>::iterator sit = it->second.begin(); sit != it->second.end(); ++sit) {
                for (unsigned int i=0; i<sit->size(); ++i) {
                    if (i>0) outfile << " ";
                    outfile << (*sit)[i];
                }
                outfile << endl;
            }
        }
        outfile.close();

    }
    
    statCfile.close();
    statUfile.close();
    statPfile.close();
    statDfile.close();


    gettimeofday(&tv,0);
    totalTime = tv.tv_sec + tv.tv_usec * 1e-6 - totalTime;
    cout << "Total time spent: " << totalTime << endl;

    return 0;
}
