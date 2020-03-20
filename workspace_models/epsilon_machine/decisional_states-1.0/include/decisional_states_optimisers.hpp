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

#ifndef DECISIONAL_STATES_OPTIMISERS_H
#define DECISIONAL_STATES_OPTIMISERS_H

#include "decisional_states_samplers.hpp"
#include "decisional_states_kernels.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/random.hpp>

#include <set>

#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace decisional_states {

template <class _PredictionType>
struct BasicPredictionSet {
    typedef _PredictionType PredictionType;
    BasicPredictionSet() : ptr(new std::set<PredictionType>()) {}
    BasicPredictionSet(std::set<PredictionType>* s) : ptr(s) {}
    BasicPredictionSet(const BasicPredictionSet& p) {ptr = p.ptr;}
    boost::shared_ptr< std::set<PredictionType> > ptr;
    // specialisation: compare contents, not pointers
    bool operator==(const BasicPredictionSet& p) const {return *ptr == *p.ptr;}
    std::set<PredictionType>* operator->() const {return ptr.operator->();}
    typedef typename std::set<PredictionType>::iterator iterator;
    typedef typename std::set<PredictionType>::const_iterator const_iterator;
    iterator begin() {return ptr->begin();}
    iterator end() {return ptr->end();}
    const_iterator begin() const {return ptr->begin();}
    const_iterator end() const {return ptr->end();}
    std::size_t size() const {return ptr->size();}
    void aggregate(const BasicPredictionSet& other) {
        ptr->insert(other.ptr->begin(), other.ptr->end());
    }
};

/// Exhaustive Optimiser tries all possible values in the given range/sampler
/// Parallel computations are supported for multi-cores with OpenMP
template <class Sampler>
struct ExhaustiveOptimiser {

    typedef typename helpers::SampleTypeWrapper<Sampler>::SampleType PredictionType;

    Sampler& range;
    ExhaustiveOptimiser(Sampler& _range) : range(_range) {}


    typedef BasicPredictionSet<PredictionType> PredictionSet;

    template<class Functor> double optimise(Functor f, PredictionSet& res) {
        double bestvalue = -std::numeric_limits<double>::max();
        res = PredictionSet(new std::set<PredictionType>());
        // exhaustive search in range
        int iend = range.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i=0; i<iend; ++i) {
            // costly part out of critical
            double result = f(range[i]);
            // critical part
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                if (result > bestvalue) {
                    res->clear();
                    res->insert(range[i]);
                    bestvalue = result;
                }
                else if (result == bestvalue) {
                    res->insert(range[i]);
                }
            }
        }
        return bestvalue;
    }
};

/// Simple Genetic algorithm optimiser:
/// - Use the provided monte-carlo sampler for the initial population
/// - Mutation is done by gaussian kernel around individuals
/// - Crossover is done by averaging
/// - Some of the best individuals can be (default 5%) preserved each generation untouched
/// - The first top of the population (1/3 default) mates and the worse individuals are replaced (fixed population)
/// - The mutation/kernel size can optionally (default) slowly decrease, like simulated annealing
/// Parallel computations are supported for multi-cores with OpenMP
template <class _PredictionType>
struct GAOptimiser {

    typedef _PredictionType PredictionType;

    MonteCarloSampler<PredictionType>& sampler;
    double mutationSize;
    int pop_size;
    int ngeneration;
    // mutation size has decayed by that ratio at the last generation: last_gen_mut = ori_mut * (1.0-mutationDecay)
    int turnover;
    int elite;
    double mutationDecay;
    boost::uniform_int<> mateDist;
    int individual_dimension;


    template <class T, class Enable = void>
    struct PredictionHelper {
        static void crossover(T& child, const T& parent1, const T& parent2) {
            for (unsigned int j=0; j<parent1.size(); ++j) child[j] = (parent1[j] + parent2[j]) * 0.5;
        }
        static unsigned int dim(const T& t) {return t.size();}
    };

    template <class T>
    struct PredictionHelper<T, typename boost::enable_if<boost::is_arithmetic<T> >::type> {
        static void crossover(T& child, const T& parent1, const T& parent2) {
            child = (parent1 + parent2) * 0.5;
        }
        static unsigned int dim(const T&) {return 1;}
    };


    GAOptimiser(MonteCarloSampler<PredictionType>& _sampler, double _mutationSize, int _pop_size = 80, int _ngeneration = 40, double _turnoverRatio = 0.3333333, double _eliteRatio = 0.05, double _mutationDecay = 0.1)
    : sampler(_sampler), mutationSize(_mutationSize), pop_size(_pop_size), ngeneration(_ngeneration), turnover((int)(_turnoverRatio * _pop_size)), elite((int)(_eliteRatio * _pop_size)), mutationDecay(_mutationDecay), mateDist(0,(int)(_turnoverRatio * _pop_size)) {
        assert(turnover < pop_size / 2);
        assert(elite < turnover);
        PredictionType p;
        individual_dimension = PredictionHelper<PredictionType>::dim(p);
    }

    typedef BasicPredictionSet<PredictionType> PredictionSet;

    struct IndexedValue {
        int idx;
        double value;
        bool operator>(const IndexedValue& other) const {
            return value > other.value;
        }
    };

    template<class Functor> double optimise(Functor f, PredictionSet& res) {
        // create initial population
        PredictionType* population = new PredictionType[pop_size];
        IndexedValue* scores = new IndexedValue[pop_size];
        SimpleGaussianKernel<double> orikernel(mutationSize,individual_dimension);
        boost::variate_generator<boost::mt19937&, boost::uniform_int<> > mateSelector(*sampler.rng,mateDist);
        for (int i=0; i<pop_size; ++i) {
            sampler.getSample(population[i],orikernel);
            scores[i].idx = i;
        }
        for (int g=0; g<ngeneration; ++g) {
            int ibeg = 0;
            // Only for non-init generations
            if (g>0) {
                //std::set<MatingPair> mates; // never mind dup pairs, count on mutation
                // Step 1: crossover for best individuals
                for (int i=0; i<turnover; ++i) {
                    int m;
                    do { m = mateSelector(); } while(m==i);
                    PredictionHelper<PredictionType>::crossover(population[scores[pop_size-1-i].idx], population[scores[i].idx], population[scores[m].idx]);
                }
                // Step 2: mutation
                SimpleGaussianKernel<double> kernel(mutationSize*(1.0-g*mutationDecay/ngeneration),individual_dimension);
                PredictionType worker;
                for (int i=elite; i<pop_size; ++i) {
                    kernel.sample(population[scores[i].idx],worker,*sampler.rng);
                    population[scores[i].idx] = worker;
                }
                // elites did not mutate: don't compute their score again
                ibeg = elite;
            }
            // Step 3: evaluation
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i=ibeg; i<pop_size; ++i) {
                scores[i].value = f(population[scores[i].idx]);
            }
            // Step 4: selection
            std::sort(scores, scores+pop_size, std::greater<IndexedValue>());
        }
        int nbest = 1;
        double bestres = scores[0].value;
        for (int i=1; i<pop_size && bestres==scores[i].value; ++i) ++nbest;
        res = PredictionSet(new std::set<PredictionType>());
        for (int i=0; i<nbest; ++i) {
//std::cout << i << "   " << population[scores[i].idx] << endl;
//            for (int j=0; j<future_samples; ++j) population[scores[i].idx][j] = (int)round(population[scores[i].idx][j]);
            res->insert(population[scores[i].idx]); // set ensures that dups in the populations are represented only once
        }
        delete [] scores;
        delete [] population;
        return bestres;
    }
};


/// NULL Optimiser in case we are only computing the epsilon-machine, it is not used
struct NullOptimiser {
    typedef int PredictionSet;
// force error if using the function
//    template<class Functor> double optimise(Functor f, PredictionSet& res) {return 0;}
};

}

#endif
