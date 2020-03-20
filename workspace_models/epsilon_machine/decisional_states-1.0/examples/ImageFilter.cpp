#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <boost/gil/gil_all.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/gil/extension/io/jpeg_dynamic_io.hpp>
#include <boost/gil/extension/io/png_dynamic_io.hpp>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/unordered_set.hpp>

#include "decisional_states.hpp"
#include "decisional_states_joint_data_manager.hpp"
#include "decisional_states_optimisers.hpp"
#include "decisional_states_connected_clustering.hpp"

#include <stdlib.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_ANN
#include "ANN_neighbour_finder.hpp"
#endif

using namespace std;
using namespace decisional_states;
using namespace boost::gil;

typedef float FloatType;
int shiftgray = 1;

namespace boost {
template<std::size_t N>
size_t hash_value(const array<uint8_t, N>& a) {
    size_t seed = 0;
    for(std::size_t i=0; i<N; ++i) boost::hash_combine(seed, a[i]);
    return seed;
}
}

#define DEFINE_BOOST_HASH(Neigh) \
size_t hash_value(const Neigh& a) { \
    size_t seed = 0; \
    for(int i=0; i<Neigh::static_size; ++i) boost::hash_combine(seed, a[i]); \
    return seed; \
}

template<int arraySize>
struct DataWrapper {
    typedef boost::array<uint8_t,arraySize> ArrayType;
    ArrayType joint;
    static const int static_size = arraySize - 1;
    inline int size() const {return arraySize-1;}
    inline const uint8_t& operator[](int idx) const {return joint[idx];}
    inline uint8_t& operator[](int idx) {return joint[idx];}
    bool operator==(const DataWrapper& other) const {
        for (int i=0; i<arraySize-1; ++i) if (other.joint[i] != joint[i]) return false;
        return true;
    }
};

template<int arraySize>
size_t hash_value(const DataWrapper<arraySize>& a) {
    size_t seed = 0;
    for (int i=0; i<(arraySize/2); ++i) boost::hash_combine(seed, a.joint[i]);
    for (int i=(arraySize/2)+1; i<arraySize; ++i) boost::hash_combine(seed, a.joint[i]);
    return seed;
}

using namespace boost;

/// Base class for all joint neighbourhoods, using the Mixin pattern
template<class _Derived, int jointsize>
struct JointNeighbourHood {
    static const int static_size = jointsize;
    static const int element_size = 1; // uint8_t = 1 byte
    static const int predPos = jointsize / 2;

    typedef _Derived Derived;
    typedef DataWrapper<jointsize> DataType;
    typedef uint8_t PredictionType;
    typedef array<uint8_t,static_size> ArrayType;
    ArrayType joint;

    static Derived at(int x, int y, gray8_view_t& imgview) {
        Derived ret;
        ret.at_impl(x,y,imgview);
        swap(ret.joint[predPos],ret.joint[jointsize-1]);
        if (shiftgray==1) {
            int minvalue = 256;
            for (int i=0; i<static_size-1; ++i) minvalue = min(minvalue, (int)ret.joint[i]);
            for (int i=0; i<static_size; ++i) ret.joint[i] -= minvalue;
        }
        return ret;
    }

    // utility for derived classes, should be protected
    uint8_t getPixel(int x, int y, gray8_view_t& imgview) {
        if (x<0) x = -x;
        if (y<0) y = -y;
        if (x>=imgview.width()) x = (imgview.width()-1)*2 - x;
        if (y>=imgview.height()) y = (imgview.height()-1)*2 - y;
        return *imgview.xy_at(x,y);
    }

    /// Joint concept API
    const DataType& getData() const {
        // same memory layout: joint and first field in struct
        return *reinterpret_cast<const DataType*>(&joint);
    }

    const PredictionType& getPrediction() const {
        return joint[jointsize-1];
    }

    void setData(const DataType& data) {
        for (int i=0; i<jointsize-1; ++i) joint[i] = (*data.array)[i];
    }
    void setPrediction(const PredictionType& prediction) {
        joint[jointsize-1] = prediction;
    }

    inline int size() const {return static_size;}
    inline const uint8_t& operator[](int idx) const {return joint[idx];}
    inline uint8_t& operator[](int idx) {return joint[idx];}
    bool operator==(const JointNeighbourHood& other) const {
        for (int i=0; i<static_size; ++i) if (other[i] != joint[i]) return false;
        return true;
    }
};

/** Horizontal and vertical immediate neighbourhood, no diagonal */
struct SmallNeighbourHood : public JointNeighbourHood<SmallNeighbourHood,5> {
    void at_impl(int x, int y, gray8_view_t& imgview) {
        joint[0] = getPixel(x, y-1, imgview);
        joint[1] = getPixel(x-1, y, imgview);
        joint[2] = getPixel(x, y, imgview);
        joint[3] = getPixel(x+1, y, imgview);
        joint[4] = getPixel(x, y+1, imgview);
    }
};
DEFINE_BOOST_HASH(SmallNeighbourHood)

/** Square neighbourhood 3x3 pixels */
struct SquareNeighbourHood : public JointNeighbourHood<SquareNeighbourHood,9> {
    void at_impl(int x, int y, gray8_view_t& imgview) {
        for (int j=0; j<3; ++j) for (int i=0; i<3; ++i) {
            joint[j*3+i] = getPixel(x+i-1, y+j-1, imgview);
        }
    }
};
DEFINE_BOOST_HASH(SquareNeighbourHood)

/** Rounded neighbourhood 5x5 pixels without corners */
struct RoundedNeighbourHood : public JointNeighbourHood<RoundedNeighbourHood,21> {
    void at_impl(int x, int y, gray8_view_t& imgview) {
        int idx = 0;
        for (int j=0; j<5; ++j) for (int i=0; i<5; ++i) {
            if ((i==0 || i==4) && (j==0 || j==4)) continue; // corners
            joint[idx++] = getPixel(x+i-2, y+j-2, imgview);
        }
    }
};
DEFINE_BOOST_HASH(RoundedNeighbourHood)

/** Rounded neighbourhood 7x7 pixels without corners
x x + + + x x
x + + + + + x
+ + + + + + +
+ + + o + + +
+ + + + + + +
x + + + + + x
x x + + + x x
*/
struct LargeRoundedNeighbourHood : public JointNeighbourHood<LargeRoundedNeighbourHood,37> {
    void at_impl(int x, int y, gray8_view_t& imgview) {
        int idx = 0;
        for (int i=-1; i<=1; ++i) joint[idx++] = getPixel(x+i, y-3, imgview);
        for (int i=-2; i<=2; ++i) joint[idx++] = getPixel(x+i, y-2, imgview);
        for (int j=-1; j<=1; ++j) for (int i=-3; i<=3; ++i) joint[idx++] = getPixel(x+i, y+j, imgview);
        for (int i=-2; i<=2; ++i) joint[idx++] = getPixel(x+i, y+2, imgview);
        for (int i=-1; i<=1; ++i) joint[idx++] = getPixel(x+i, y+3, imgview);
    }
};
DEFINE_BOOST_HASH(LargeRoundedNeighbourHood)

template<typename T>
struct NeighbourHoodHash : public std::unary_function<T, std::size_t> {
    std::size_t operator()(T const& a) const {
        std::size_t seed = 0;
        for(int i=0; i<T::static_size; ++i) boost::hash_combine(seed, a[i]);
        return seed;
    }
};


template<class VecType>
inline FloatType vector_dist2(const VecType& a, const VecType& b) {
    FloatType ret1 = 0, ret2 = 0;
    for (int i=0; i<(VecType::static_size&(-2)); i+=2) {
        ret1 += (a[i] - b[i])*(a[i] - b[i]);
        ret2 += (a[i+1] - b[i+1])*(a[i+1] - b[i+1]);
    }
    if (VecType::static_size&1) ret1+=(a[VecType::static_size-1] - b[VecType::static_size-1])*(a[VecType::static_size-1] - b[VecType::static_size-1]);
    return ret1+ret2;
}

template<class VecType>
inline FloatType vector_distance(const VecType& a, const VecType& b) {
    return sqrtf(vector_dist2(a,b));
}


static int spacialVariance = 10000;

namespace boost {
size_t hash_value(const boost::array<FloatType, 2>& a) {
    size_t seed = 0;
    boost::hash_combine(seed, a[0]);
    boost::hash_combine(seed, a[1]);
    return seed;
}
}

template<typename JointType, typename FloatType>
struct ImageJointSampler {
    // sampler in image space
    SimpleGaussianKernel<FloatType> spacialKernel;

    template<class RNG>
    ImageJointSampler(int x, int y, int _nsamples, RNG& rng, gray8_view_t& imgview) : spacialKernel(spacialVariance, 2), nsamples(_nsamples) {
        samples = boost::shared_array<JointType>(new JointType[nsamples]);
        boost::array<FloatType, 2> center; center[0] = x+0.5f; center[1] = y+0.5f;
        boost::array<FloatType, 2> xysample;
        boost::unordered_set<boost::array<FloatType, 2> > samplesloc;
        for (int i=0; i<nsamples; ++i) {
            do {
                spacialKernel.sample(center, xysample, rng);
            } while (xysample[0] < 0 || xysample[0] >= imgview.width() || xysample[1] < 0 || xysample[1] >= imgview.height() || (samplesloc.find(xysample)!=samplesloc.end()));
            samplesloc.insert(xysample);
            samples[i] = JointType::at((int)floorf(xysample[0]),(int)floorf(xysample[1]),imgview);
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


template<class DataSet, typename DataType>
struct CustomTransitionFeeder {
    DataSet& dataset;
    int the_size;
    int nhoriz;
    int nvert;
    
    CustomTransitionFeeder(DataSet& _dataset, int _nhoriz, int _nvert) : dataset(_dataset), nhoriz(_nhoriz), nvert(_nvert) {
        // left-to-right, right-to-left, top-to-bottom, bottom-to-top
        the_size = ((nhoriz-1) * nvert)*2 + (nhoriz * (nvert-1))*2;
    }
    
    int size() const {return the_size;}
    
    typedef int iterator;
    iterator begin() {return 0;}
    iterator end() {return the_size;}
    
    DataType getDataBeforeTransition(iterator idx) {
        // left-to-right
        if (idx < (nhoriz-1) * nvert) {
            int v = idx / (nhoriz-1);
            int h = idx % (nhoriz-1);
            // before transition = left
            return dataset.data(v * nhoriz + h);
        }
        idx -= (nhoriz-1) * nvert;
        // right-to-left
        if (idx < (nhoriz-1) * nvert) {
            int v = idx / (nhoriz-1);
            int h = idx % (nhoriz-1);
            // before transition = right
            return dataset.data(v * nhoriz + h + 1);
        }
        idx -= (nhoriz-1) * nvert;
        // top-to-bottom
        if (idx < nhoriz * (nvert-1)) {
            int v = idx / nhoriz;
            int h = idx % nhoriz;
            // before transition = top
            return dataset.data(v * nhoriz + h);
        }
        idx -= nhoriz * (nvert-1);
        // bottom-to-top
        int v = idx / nhoriz;
        int h = idx % nhoriz;
        // before transition = bottom
        return dataset.data((v+1) * nhoriz + h);
    }
    
    DataType getDataAfterTransition(iterator idx) {
        // left-to-right
        if (idx < (nhoriz-1) * nvert) {
            int v = idx / (nhoriz-1);
            int h = idx % (nhoriz-1);
            // after transition = right
            return dataset.data(v * nhoriz + h + 1);
        }
        idx -= (nhoriz-1) * nvert;
        // right-to-left
        if (idx < (nhoriz-1) * nvert) {
            int v = idx / (nhoriz-1);
            int h = idx % (nhoriz-1);
            // after transition = left
            return dataset.data(v * nhoriz + h);
        }
        idx -= (nhoriz-1) * nvert;
        // top-to-bottom
        if (idx < nhoriz * (nvert-1)) {
            int v = idx / nhoriz;
            int h = idx % nhoriz;
            // after transition = bottom
            return dataset.data((v+1) * nhoriz + h);
        }
        idx -= nhoriz * (nvert-1);
        // bottom-to-top
        int v = idx / nhoriz;
        int h = idx % nhoriz;
        // after transition = top
        return dataset.data(v * nhoriz + h);
    }

};


struct TypeTraits {
    typedef ::FloatType FloatType;

    //typedef SquareNeighbourHood JointType;
    typedef RoundedNeighbourHood JointType;

    typedef JointType::DataType DataType;
    typedef JointType::PredictionType PredictionType;

    typedef RegularGridSampler<FloatType,PredictionType> PredictionSampler;

    typedef SimpleGaussianKernel<FloatType> JointKernel;

#ifdef USE_ANN
    typedef ANN_NeighbourhoodFinder<DataType> DataNeighbourhoodFinder;
    typedef ANN_NeighbourhoodFinder<PredictionType> PredictionNeighbourhoodFinder;
    typedef ANN_NeighbourhoodFinder<JointType> JointNeighbourhoodFinder;
#else
    typedef NearTreeNeighbourhoodFinder<DataType, FloatType, &vector_distance<DataType>, &vector_dist2<DataType> > DataNeighbourhoodFinder;
    typedef NearTreeNeighbourhoodFinder<PredictionType, FloatType > PredictionNeighbourhoodFinder;
    typedef NearTreeNeighbourhoodFinder<JointType, FloatType, &vector_distance<JointType>, &vector_dist2<JointType> > JointNeighbourhoodFinder;
#endif

    //typedef IndexRangeDistributionStorage<FloatType> DistributionStorage;
    typedef IndexedDistributionStorage<FloatType> DistributionStorage;
    
    typedef ExhaustiveOptimiser<PredictionSampler> Optimiser;
    typedef Optimiser::PredictionSet PredictionSet;
};

typedef TypeTraits::JointType JointType;
typedef TypeTraits::DataType DataType;
typedef TypeTraits::PredictionType PredictionType;
typedef TypeTraits::PredictionSampler PredictionSampler;
typedef TypeTraits::JointKernel JointKernel;
typedef TypeTraits::JointNeighbourhoodFinder JointNeighbourhoodFinder;
typedef TypeTraits::DataNeighbourhoodFinder DataNeighbourhoodFinder;
typedef TypeTraits::PredictionNeighbourhoodFinder PredictionNeighbourhoodFinder;
typedef TypeTraits::Optimiser Optimiser;

typedef JointDataManager<TypeTraits>::DataSet DataSet;
typedef JointDataManager<TypeTraits>::DistributionManager DistributionManager;

typedef CustomTransitionFeeder<DataSet, DataType> TransitionFeeder;
typedef DecisionalStatesAnalyser<DistributionManager, DataSet, Optimiser, TransitionFeeder, TypeTraits> DSA;

typedef DistributionManager::Distribution Distribution;


int gray_tolerance = 5;

struct Utility {
    template <typename P>
    typename boost::disable_if<boost::is_arithmetic<P>, FloatType>::type
    operator()(const P& p1, const P& p2) {
        FloatType ret = 0;
        for (int i=0; i<(int)p1.size(); ++i) {
            FloatType val = fabs((FloatType)p1[i] - (FloatType)p2[i]) - gray_tolerance;
            if (val <0) val = 0;
            ret += val;
        }
        return -ret;
    }

    template <typename P>
    typename boost::enable_if<boost::is_arithmetic<P>, FloatType>::type
    operator()(const P& p1, const P& p2) {
        FloatType ret = fabs((FloatType)p1 - (FloatType)p2) - gray_tolerance;
        if (ret <0) ret = 0;
        return -ret;
    }
};
/*
struct SQ_Utility {
    template <typename P> FloatType operator()(const P& p1, const P& p2) {
        return -vector_dist2<P>(p1,p2);
    }
};
*/


template<class State>
struct StateSorter {
    bool operator()(typename boost::shared_ptr<State> a, typename boost::shared_ptr<State> b) {
        // less count = more complex
        // equal counts but more members
        return ((a->count) < (b->count)) || (((a->count) == (b->count)) && ((a->members->size()) > (b->members->size())));
    }
};

struct Identity {
template<class T> inline T operator()(T t) const {return t;}
};


struct StateUserInfo {
    int rank;
    int num;
    int avgPix;
    StateUserInfo(int r, int n) : rank(r), num(n) {}
};

FloatType kernelVar = 0;


void setHuePixel(float hue, rgb8_pixel_t& pixel) {
    hue = hue - floorf(hue); // 0 <= hue < 1
    int r,g,b;
    if (hue*6.0f < 1.0f) {
        r=255; b=0;
        g = (int)(255.99f * hue * 6.0f);
    }
    else if (hue*6.0f < 2.0f) {
        g=255; b=0;
        r = (int)(255.99f * (1.0f-hue) * 6.0f);
    }
    else if (hue*6.0f < 3.0f) {
        g=255; r=0;
        b = (int)(255.99f * (hue-2.0f) * 6.0f);
    }
    else if (hue*6.0f < 4.0f) {
        b=255; r=0;
        g = (int)(255.99f * (3.0f-hue) * 6.0f);
    }
    else if (hue*6.0f < 5.0f) {
        b=255; g=0;
        r = (int)(255.99f * (hue-4.0f) * 6.0f);
    }
    else {
        r=255; g=0;
        b = (int)(255.99f * (5.0f-hue) * 6.0f);
    }
    pixel[0] = r;
    pixel[1] = g;
    pixel[2] = b;
}

int main(int argc, char**argv) {
    if (argc<6) {
        cout << "Arguments expected:  kernel_width_or_0  gray_tolerance shift_gray_or_not output_dir  input_file1  [input_file2 [...]]" << endl;
        return 1;
    }

    timeval tv; gettimeofday(&tv,0);
    double totalTime = tv.tv_sec + tv.tv_usec * 1e-6;

    int argi = 0;
    kernelVar = atof(argv[++argi]);
    gray_tolerance = atoi(argv[++argi]);
    shiftgray = atoi(argv[++argi]);
    string output_dir = argv[++argi];

    boost::mt19937 rng;
    rng.seed(42);

    PredictionSampler psampler(0,255,256); // fixed
    JointDataManager<TypeTraits> jdm(psampler);
    int argimg = ++argi;
    
    gettimeofday(&tv,0);
    double theTime = tv.tv_sec + tv.tv_usec * 1e-6;

    int nhoriz = 0, nvert = 0;
    
    cout << "Loading images and building probability distributions..." << endl;
    for (argi = argimg; argi<argc; ++argi) {
        typedef mpl::vector<gray8_image_t, gray16_image_t, rgb8_image_t, rgb16_image_t> my_img_types;
        any_image<my_img_types> runtime_image;
        bool ispng = true;
        try {
            png_read_image(argv[argi], runtime_image);
        } catch(...) {
            try {
                jpeg_read_image(argv[argi], runtime_image);
                ispng = false;
            } catch(...) {
                cerr << "Cannot load " << argv[argi] << " as either png or jpeg image. Check it exists and it is 8 or 16 gray or rgb format." << endl;
                return 2;
            }
        }
        gray8_image_t grayimage(runtime_image.dimensions());
        gray8_view_t grayview = view(grayimage);
        copy_pixels(color_converted_view<gray8_pixel_t>(const_view(runtime_image)), grayview);
        JointKernel jointKernel(0,JointType::static_size,1e-5f); // variance is adapted locally

        nhoriz = runtime_image.width();
        nvert = runtime_image.height();
        
        spacialVariance = std::min(runtime_image.width(), runtime_image.height()) / 2;
        spacialVariance *= spacialVariance;

        int npix = runtime_image.height() * runtime_image.width();

        std::vector<JointType> allPoints(npix);
        for (int y=0; y<runtime_image.height(); ++y) for (int x=0; x<runtime_image.width(); ++x) allPoints[y*runtime_image.width()+x] = JointType::at(x,y,grayview);

        // Or adapt variance globally, but only once ? TODO: param
        if (kernelVar>0) jointKernel.setSize(kernelVar);
        else {
            JointNeighbourhoodFinder jnf;
            kernelVar = jointKernel.setSizeFromSamples(allPoints, jnf);
            cout << "kvar = " << kernelVar << endl;
        }
        
        SimpleGaussianKernel<FloatType> pkernel(kernelVar,1,1e-9f);
        SimpleGaussianKernel<FloatType> dkernel(kernelVar,DataType::static_size,1e-9f);
        
        jdm.addSeparable(allPoints,dkernel,pkernel,allPoints);

    }
    cout << "Building analyser" << endl;

    // OK, now we have a fully loaded JointDataManager
    // perform the clusterings to get the states
    Optimiser optimiser(psampler);
    
    // Note: feeder for only one image. TODO.
    TransitionFeeder transitionFeeder(jdm.asDataSet(), nhoriz, nvert);
    DSA analyser(jdm.asDistributionManager(), transitionFeeder, optimiser);
    
    gettimeofday(&tv,0);
    theTime = tv.tv_sec + tv.tv_usec * 1e-6 - theTime;
    cout << "time spent for building stats: " << theTime << endl;

    bool consistent;
    cout << "Computing causal states... " << flush;
    gettimeofday(&tv,0);
    theTime = tv.tv_sec + tv.tv_usec * 1e-6;
    
    consistent = analyser.computeCausalStates();
    
    cout << "Num causal states: " << analyser.causalStates.size() << ", consistent flag=" << consistent << endl;

    gettimeofday(&tv,0);
    theTime = tv.tv_sec + tv.tv_usec * 1e-6 - theTime;
    cout << "time spent for computing causal states: " << theTime << endl;

    gettimeofday(&tv,0);
    theTime = tv.tv_sec + tv.tv_usec * 1e-6;

    cout << "Applying utility" << flush;
        cout << " (gray_tolerance="<<gray_tolerance<<")..." << flush;
        analyser.applyUtility( Utility() );
    cout << " done!" << endl;

    gettimeofday(&tv,0);
    theTime = tv.tv_sec + tv.tv_usec * 1e-6 - theTime;
    cout << "time spent for applying utility: " << theTime << endl;

    consistent = analyser.computeIsoPredictionStates();
    cout << "Num iso prediction states: " << analyser.isoPredictionStates.size() << ", consistent flag=" << consistent << endl;

    consistent = analyser.computeIsoUtilityStates();
    cout << "Num iso utility states: " << analyser.isoUtilityStates.size() << ", consistent flag=" << consistent << endl;

    analyser.computeDecisionalStates();
    cout << "Num decisional states: " << analyser.decisionalStates.size() << endl;

    //analyser.countClusters();
    analyser.buildCausalStateGraph();

    ofstream dot((output_dir+"/causal_states_graph.dot").c_str());
    analyser.writeCausalStateGraph(dot);
    dot.close();
    dot.open((output_dir+"/causal_states_graph_with_transients.dot").c_str());
    analyser.writeCausalStateGraph(dot,false);
    dot.close();
    
    analyser.buildIsoPredictionStateGraph();

    dot.open((output_dir+"/iso_prediction_states_graph.dot").c_str());
    analyser.writeIsoPredictionStateGraph(dot);
    dot.close();
    dot.open((output_dir+"/iso_prediction_states_graph_with_transients.dot").c_str());
    analyser.writeIsoPredictionStateGraph(dot,false);
    dot.close();
    
    analyser.buildIsoUtilityStateGraph();

    dot.open((output_dir+"/iso_utility_states_graph.dot").c_str());
    analyser.writeIsoUtilityStateGraph(dot);
    dot.close();
    dot.open((output_dir+"/iso_utility_states_graph_with_transients.dot").c_str());
    analyser.writeIsoUtilityStateGraph(dot,false);
    dot.close();
    
    analyser.buildDecisionalStateGraph();

    dot.open((output_dir+"/decisional_states_graph.dot").c_str());
    analyser.writeDecisionalStateGraph(dot);
    dot.close();
    dot.open((output_dir+"/decisional_states_graph_with_transients.dot").c_str());
    analyser.writeDecisionalStateGraph(dot,false);
    dot.close();
    

    cout << "Outputting result images in the directory: " << output_dir << endl;

    // shared ptr: user field modified in original and copy
    int curcount = -1;
    DSA::CausalStates statesC = analyser.causalStates;
    sort(statesC.begin(), statesC.end(), StateSorter<DSA::CausalState>());
    long maxrankC = -1;
    for (int i=0; i<(int)statesC.size(); ++i) {
        if (statesC[i]->count != curcount) {
            curcount = statesC[i]->count;
            ++maxrankC;
        }
        statesC[i]->user = (void*)maxrankC;
    }
    curcount = -1;
    DSA::DecisionalStates statesD = analyser.decisionalStates;
    sort(statesD.begin(), statesD.end(), StateSorter<DSA::DecisionalState>());
    long maxrankD = -1;
    for (int i=0; i<(int)statesD.size(); ++i) {
        if (statesD[i]->count != curcount) {
            curcount = statesD[i]->count;
            ++maxrankD;
        }
        statesD[i]->user = (void*)maxrankD;
    }
    curcount = -1;
    DSA::IsoUtilityStates statesU = analyser.isoUtilityStates;
    sort(statesU.begin(), statesU.end(), StateSorter<DSA::IsoUtilityState>());
    long maxrankU = -1;
    for (int i=0; i<(int)statesU.size(); ++i) {
        if (statesU[i]->count != curcount) {
            curcount = statesU[i]->count;
            ++maxrankU;
        }
        statesU[i]->user = (void*)maxrankU;
    }
    curcount = -1;
    DSA::IsoPredictionStates statesP = analyser.isoPredictionStates;
    sort(statesP.begin(), statesP.end(), StateSorter<DSA::IsoPredictionState>());
    long maxrankP = -1;
    for (int i=0; i<(int)statesP.size(); ++i) {
        if (statesP[i]->count != curcount) {
            curcount = statesP[i]->count;
            ++maxrankP;
        }
        statesP[i]->user = (void*)maxrankP;
    }

    DataSet::iterator datait = jdm.asDataSet().begin();

    for (argi = argimg; argi<argc; ++argi) {
        typedef mpl::vector<gray8_image_t, gray16_image_t, rgb8_image_t, rgb16_image_t> my_img_types;
        any_image<my_img_types> runtime_image;
        bool ispng = true;
        try {
            png_read_image(argv[argi], runtime_image);
        } catch(...) {
            try {
                jpeg_read_image(argv[argi], runtime_image);
                ispng = false;
            } catch(...) {
                cerr << "Cannot load " << argv[argi] << " as either png or jpeg image. Check it exists and it is 8 or 16 gray or rgb format." << endl;
                return 2;
            }
        }
        gray8_image_t grayimage(runtime_image.dimensions());
        gray8_view_t grayview = view(grayimage);
        copy_pixels(color_converted_view<gray8_pixel_t>(const_view(runtime_image)), grayview);

        gray8_image_t Cimage(runtime_image.dimensions());
        gray8_view_t Cview = view(Cimage);
        gray8_image_t Dimage(runtime_image.dimensions());
        gray8_view_t Dview = view(Dimage);
        gray8_image_t Uimage(runtime_image.dimensions());
        gray8_view_t Uview = view(Uimage);
        gray8_image_t Pimage(runtime_image.dimensions());
        gray8_view_t Pview = view(Pimage);
        rgb8_image_t CimageColor(runtime_image.dimensions());
        rgb8_view_t CviewColor = view(CimageColor);
        rgb8_image_t DimageColor(runtime_image.dimensions());
        rgb8_view_t DviewColor = view(DimageColor);
        rgb8_image_t UimageColor(runtime_image.dimensions());
        rgb8_view_t UviewColor = view(UimageColor);
        rgb8_image_t PimageColor(runtime_image.dimensions());
        rgb8_view_t PviewColor = view(PimageColor);
        gray8_image_t Cimage_rank(runtime_image.dimensions());
        gray8_view_t Cview_rank = view(Cimage_rank);
        gray8_image_t Dimage_rank(runtime_image.dimensions());
        gray8_view_t Dview_rank = view(Dimage_rank);
        gray8_image_t Uimage_rank(runtime_image.dimensions());
        gray8_view_t Uview_rank = view(Uimage_rank);
        gray8_image_t Pimage_rank(runtime_image.dimensions());
        gray8_view_t Pview_rank = view(Pimage_rank);

        FloatType logCmax = helpers::log2((FloatType)analyser.causalStatesCounts.maxcount);
        FloatType deltalogC = logCmax - helpers::log2((FloatType)analyser.causalStatesCounts.mincount);
        FloatType logDmax = helpers::log2((FloatType)analyser.decisionalStatesCounts.maxcount);
        FloatType deltalogD = logDmax - helpers::log2((FloatType)analyser.decisionalStatesCounts.mincount);
        FloatType logUmax = helpers::log2((FloatType)analyser.isoUtilityStatesCounts.maxcount);
        FloatType deltalogU = logUmax - helpers::log2((FloatType)analyser.isoUtilityStatesCounts.mincount);
        FloatType logPmax = helpers::log2((FloatType)analyser.isoPredictionStatesCounts.maxcount);
        FloatType deltalogP = logPmax - helpers::log2((FloatType)analyser.isoPredictionStatesCounts.mincount);

        for (int y=0; y<runtime_image.height(); ++y) for (int x=0; x<runtime_image.width(); ++x) {
            assert(datait!=jdm.asDataSet().end());

            DataType d = jdm.asDataSet().data(datait);

            shared_ptr<DSA::CausalState> causalState = analyser.getCausalState(d);
            shared_ptr<DSA::DecisionalState> decisionalState = analyser.getDecisionalState(d);
            shared_ptr<DSA::IsoUtilityState> isoUtilityState = analyser.getIsoUtilityState(d);
            shared_ptr<DSA::IsoPredictionState> isoPredictionState = analyser.getIsoPredictionState(d);

            FloatType complexityC = (logCmax-helpers::log2((FloatType)causalState->count)) / deltalogC;
            FloatType complexityD = (logDmax-helpers::log2((FloatType)decisionalState->count)) / deltalogD;
            FloatType complexityU = (logUmax-helpers::log2((FloatType)isoUtilityState->count)) / deltalogU;
            FloatType complexityP = (logPmax-helpers::log2((FloatType)isoPredictionState->count)) / deltalogP;

            *Cview.xy_at(x, y) = max(0,min(255, (int)((1.0-complexityC)*255.99)));
            *Dview.xy_at(x, y) = max(0,min(255, (int)((1.0-complexityD)*255.99)));
            *Uview.xy_at(x, y) = max(0,min(255, (int)((1.0-complexityU)*255.99)));
            *Pview.xy_at(x, y) = max(0,min(255, (int)((1.0-complexityP)*255.99)));

            setHuePixel((complexityC + 1) * 2 / 3, *CviewColor.xy_at(x, y));
            setHuePixel((complexityD + 1) * 2 / 3, *DviewColor.xy_at(x, y));
            setHuePixel((complexityU + 1) * 2 / 3, *UviewColor.xy_at(x, y));
            setHuePixel((complexityP + 1) * 2 / 3, *PviewColor.xy_at(x, y));

            *Cview_rank.xy_at(x, y) = max(0,min(255, (int)((long)causalState->user/(FloatType)maxrankC*255.99)));
            *Dview_rank.xy_at(x, y) = max(0,min(255, (int)((long)decisionalState->user/(FloatType)maxrankD*255.99)));
            *Uview_rank.xy_at(x, y) = max(0,min(255, (int)((long)isoUtilityState->user/(FloatType)maxrankU*255.99)));
            *Pview_rank.xy_at(x, y) = max(0,min(255, (int)((long)isoPredictionState->user/(FloatType)maxrankP*255.99)));


            ++datait;
        }

        string filename = output_dir+"/compC_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),Cview);
        filename = output_dir+"/compD_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),Dview);
        filename = output_dir+"/compU_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),Uview);
        filename = output_dir+"/compP_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),Pview);

        filename = output_dir+"/compC_color_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),CviewColor);
        filename = output_dir+"/compD_color_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),DviewColor);
        filename = output_dir+"/compU_color_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),UviewColor);
        filename = output_dir+"/compP_color_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),PviewColor);

        filename = output_dir+"/compC_rank_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),Cview_rank);
        filename = output_dir+"/compD_rank_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),Dview_rank);
        filename = output_dir+"/compU_rank_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),Uview_rank);
        filename = output_dir+"/compP_rank_"+argv[argi];
        png_write_view(filename.replace(filename.end()-3, filename.end(), "png"),Pview_rank);

    }

    gettimeofday(&tv,0);
    totalTime = tv.tv_sec + tv.tv_usec * 1e-6 - totalTime;
    cout << "Total time spent in main processing: " << totalTime << endl;

    return 0;
}



