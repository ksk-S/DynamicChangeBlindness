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

#ifndef DECISIONAL_STATES_SAMPLERS_H
#define DECISIONAL_STATES_SAMPLERS_H

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/key_extractors.hpp>
#include <boost/random.hpp>

#include <vector>

// DEBUG only ?
#include <iostream>
#include <limits>

#include "decisional_states_helpers.hpp"

namespace decisional_states {

// Sampling on regularly spaced points in each dimension
template<class ElementType, class PredictionType = ElementType>
struct RegularGridSampler {
    typedef PredictionType SampleType;

    boost::shared_array<PredictionType> samples;
    unsigned int ssize;

    // constructor for arithmetic prediction types (float, int...)
    template <class T>
    void fillSamples(ElementType min, ElementType max, int n, int dim, PredictionType emptyPrediction, typename boost::enable_if<boost::is_arithmetic<T> >::type* sfinae = 0) {
        ElementType range = max - min;
        ssize = n;
        samples = boost::shared_array<PredictionType>(new PredictionType[ssize]);
        int nm1 = n-1;
        assert(nm1>0);
        for (int i=0; i<n; ++i) samples[i] = min + i * range / nm1;
    }
    // constructor for more elaborate types using [] to subscript dimensions
    template <class T>
    void fillSamples(ElementType min, ElementType max, int n, int dim, PredictionType emptyPrediction, typename boost::disable_if<boost::is_arithmetic<T> >::type* sfinae = 0) {
        ElementType range = max - min;
        int nm1 = n-1;
        assert(nm1>0);
        ssize = n;
        for (int i=1; i<dim; ++i) ssize *= n;
        samples = boost::shared_array<PredictionType>(new PredictionType[ssize]);
        PredictionType currentSample = emptyPrediction;
        for(int d=0; d<dim; ++d) currentSample[d] = min;
        std::vector<int> positions(dim,0);
        for (unsigned int i=0; i<ssize; ++i) {
            // modify only what's necessary in current sample
            for(int d=0; d<dim; ++d) {
                ++positions[d];
                currentSample[d] = min + positions[d] * range / nm1;
                if (positions[d] < n) break; // done, still in the same dimension
                positions[d] = 0;            // and increase the next dimension
                currentSample[d] = min;
            }
            samples[i] = currentSample;
        }
    }

    RegularGridSampler(ElementType min, ElementType max, int n, int dim = 1, PredictionType emptyPrediction = PredictionType()) {
        fillSamples<PredictionType>(min,max,n,dim,emptyPrediction);
        endit = samples.get()+ssize;
        endcit = endit;
    }

    typedef PredictionType* iterator;
    typedef const PredictionType* const_iterator;

    iterator endit;
    const_iterator endcit;
    iterator begin() {return samples.get();}
    iterator end() {return endit;}
    const_iterator begin() const {return samples.get();}
    const_iterator end() const {return endcit;}
    PredictionType& operator[](unsigned int i) {
        return samples[i];
    }
    const PredictionType& operator[](unsigned int i) const {
        return samples[i];
    }
    unsigned int size() const {return ssize;}
};

// Sample points using a kernel around the base points in the data set
// TODO: check the implementation, possibly make a correct one
template<class PredictionType>
struct MonteCarloSampler {

    typedef PredictionType SampleType;

    struct Z {
        PredictionType* zptr;
        Z(PredictionType* _zptr) : zptr(_zptr) {}
        PredictionType& deref() const {return *zptr;}
    };
    typedef boost::multi_index_container<
        Z,
        boost::multi_index::indexed_by<
            boost::multi_index::hashed_unique< boost::multi_index::const_mem_fun<Z,PredictionType&,&Z::deref> >,
            boost::multi_index::random_access<>
        >
    > ZContainer;

    unsigned int nsamples;
    boost::shared_array<PredictionType> samples;
    boost::shared_ptr<ZContainer> basePoints;
    boost::shared_ptr<boost::mt19937> rng;
    boost::shared_ptr<boost::uniform_int<> > randomIndexDist;
    boost::shared_ptr<boost::variate_generator<boost::mt19937&, boost::uniform_int<> > > randomIndex;

    template<class SamplingKernel, class DataSet>
    MonteCarloSampler(int _nsamples, SamplingKernel kernel, DataSet& dataset, unsigned int seed = 42) : nsamples(_nsamples), samples(new PredictionType[_nsamples]) {
        // sample uniformly amongst all known prediction values
        // the data set distribution of these predictions is taken into account in the main algorithm
        basePoints = boost::shared_ptr<ZContainer>(new ZContainer);
        for (typename DataSet::iterator it = dataset.begin(); it != dataset.end(); ++it) {
            basePoints->insert(Z(&helpers::DataSetIteratorUnrefHelper<DataSet>::getPrediction(dataset,it)));
        }
        rng = boost::shared_ptr<boost::mt19937>(new boost::mt19937(seed));
        randomIndexDist = boost::shared_ptr<boost::uniform_int<> >(new boost::uniform_int<>(0,basePoints->template get<1>().size()-1));
        randomIndex = boost::shared_ptr<boost::variate_generator<boost::mt19937&, boost::uniform_int<> > >(new boost::variate_generator<boost::mt19937&, boost::uniform_int<> >(*rng, *randomIndexDist));

        for (unsigned int i=0; i<nsamples; ++i) {
            getSample(samples[i],kernel);
/*std::cout << "b: " << dataBase[0] << " " << dataBase[1] << std::endl;
std::cout << "s: " << samples[i][0] << " " << samples[i][1] << std::endl;
std::cout << "k: " << kernel(dataBase,samples[i]) << std::endl;
*/
        }
        endit = &samples[0]+nsamples;
        endcit = endit;
    }
    typedef PredictionType* iterator;
    typedef const PredictionType* const_iterator;

    iterator endit;
    const_iterator endcit;
    iterator begin() {return &samples[0];}
    iterator end() {return endit;}
    const_iterator begin() const {return &samples[0];}
    const_iterator end() const {return endcit;}
    PredictionType& operator[](unsigned int i) {
        return samples[i];
    }
    const PredictionType& operator[](unsigned int i) const {
        return samples[i];
    }
    unsigned int size() const {return nsamples;}

    // Expose extra utility function
    template<class SamplingKernel>
    void getSample(PredictionType& sample, SamplingKernel& kernel) {
        // choose a point at random as a base.
        const PredictionType& dataBase = *basePoints->template get<1>()[(*randomIndex)()].zptr;
        kernel.sample(dataBase,sample,*rng);
    }
};

namespace helpers {
template<class _Vector, class _Scalar>
struct VectorSpaceTrait {
    typedef _Vector Vector;
    typedef _Scalar Scalar;
    typedef void (*VectorAddFunctor)(Vector&, const Vector&,const Vector&);
    typedef void (*ScalarMultiplyFunctor)(Vector&, const Scalar&,const Vector&);
};
}
template <typename FloatType>
struct CWI {
    FloatType contrib;
    FloatType weight;
    int idx;
};

// For N sample points and some real data, this sampler adapts a N-components Gaussian mixture so it fits best the data, using the EM algorithm
template<class _SampleType, typename FloatType>
struct GaussianMixtureSampler {
    typedef _SampleType SampleType;
    int nsamples;

    boost::shared_array<SampleType> samples;

    /// The constructor performs the EM algo, so the container need not be stored and just the retained samples
    /// The random number generator is used to determine an initial position for the samples
    template<class Container, class RNG, class SamplingKernel>
    GaussianMixtureSampler(
        int _nsamples,
        Container& container,
        SamplingKernel kernel,
        RNG& rng,
        typename helpers::VectorSpaceTrait<SampleType, FloatType>::VectorAddFunctor add_functor,
        typename helpers::VectorSpaceTrait<SampleType, FloatType>::ScalarMultiplyFunctor scale_functor,
        SampleType nullSample,
        int maxEMSteps = 0
    ) : nsamples(_nsamples), samples(new SampleType[_nsamples]) {

        int ndata = container.size();

        // Just pick points distributed evenly amongst data
        // Many other mechanisms were tried, various EM initialisation techniques (including KKZ)
        // sampling relative to contributions of the data, constraint sampling
        // to the subspace explored by EM, modifying the kernel size in the EM,
        // checking that points are not too close by even during EM, doing cross-validation
        // by computing contributions from other points than the current one, etc...
        // NONE of these improve the situation, and often produce worse results.
        // See the history in the GIT repository for code snippets of all these initialisations.
        for (int j=0; j<nsamples; ++j) {
            int didx = (ndata-1) * j / (nsamples-1);
            samples[j] = container[didx];
        }

        // Now perform the EM algorithm - if requested by the user
        if (maxEMSteps==0) return;

        // the means obtained by EM will be the sample points
        std::cout << "Performing EM steps for adapting sample points to data" << std::endl;

        boost::shared_array<FloatType> normalizationFactors = boost::shared_array<FloatType>(new FloatType[ndata]);
        for (int i=0; i<ndata; ++i) normalizationFactors[i] = 0;
        for (int emstep=0; emstep<maxEMSteps; ++emstep) {
            for (int l=0; l<nsamples; ++l) {
                FloatType sump = 0;
                SampleType avgs = nullSample;
                for (int i=0; i<ndata; ++i) {
                    FloatType p = kernel(samples[l], container[i]);
                    sump += p;
                    normalizationFactors[i] += p;
                    SampleType scaled = nullSample;
                    scale_functor(scaled, p, container[i]);
                    SampleType tmp = nullSample;
                    add_functor(tmp, avgs, scaled);
                    avgs = tmp;
                }
                if (sump!=0) scale_functor(samples[l], FloatType(1.0f/sump), avgs);
            }
            // compute log likelihood up to a constant
            FloatType logl = 0;
            for (int i=0; i<ndata; ++i) {
                if (normalizationFactors[i]!=0) logl += std::log(normalizationFactors[i]);
            }
            std::cout << "EM step " << emstep << ", log likelihood (up to a constant) = " << logl << std::endl;
        }

    }

    typedef SampleType* iterator;
    typedef SampleType* const_iterator;
    iterator begin() {return &samples[0];}
    iterator end() {return &samples[nsamples-1];}
    const_iterator begin() const {return &samples[0];}
    const_iterator end() const {return &samples[nsamples-1];}
    SampleType& operator[](unsigned int i) {
        return samples[i];
    }
    const SampleType& operator[](unsigned int i) const {
        return samples[i];
    }
    unsigned int size() const {return nsamples;}

};


// Sample points around a given point using a kernel
// actually just a wrapper around the kernel "sample" method
template<class _SampleType>
struct KernelSampler {
    typedef _SampleType SampleType;

    int nsamples;

    boost::shared_array<SampleType> samples;

    template<class KernelType, class RNG>
    KernelSampler(SampleType center, KernelType kernel, int _nsamples, RNG& rng) : nsamples(_nsamples) {
        samples = boost::shared_array<SampleType>(new SampleType[nsamples]);
        for (int i=0; i<nsamples; ++i) kernel.sample(center, samples[i], rng);
    }

    typedef SampleType* iterator;
    typedef const SampleType* const_iterator;

    iterator begin() {return &samples[0];}
    iterator end() {return &samples[nsamples];}
    const_iterator begin() const {return &samples[0];}
    const_iterator end() const {return &samples[nsamples];}
    SampleType& operator[](unsigned int i) {
        return samples[i];
    }
    const SampleType& operator[](unsigned int i) const {
        return samples[i];
    }
    unsigned int size() const {return nsamples;}
    
};


}

#endif
