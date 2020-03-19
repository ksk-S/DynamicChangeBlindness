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

// TODO too : rename old discrete dist manager to sparse sth

#ifndef DECISIONAL_STATES_DISCRETE_DISTRIBUTION_MANAGER_H
#define DECISIONAL_STATES_DISCRETE_DISTRIBUTION_MANAGER_H

#include "decisional_states_helpers.hpp"

#include <vector>
//#include <map>

#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>

#include <math.h>

namespace decisional_states {

template<class _DataSet, class TypeTraits = helpers::DefaultTypeTraitsDS<_DataSet> > struct DiscreteDistributionManager {
    typedef _DataSet DataSet;
    typedef typename TypeTraits::FloatType FloatType;
    typedef typename TypeTraits::PredictionType PredictionType;
    typedef typename TypeTraits::DataType DataType;

    struct Distribution {
        // pointer is either in the globalDistCount mem zone, then not handled here
        // or it is external (ex: aggregate dist), then handled here
        DiscreteDistributionManager* manager;
        int* counts;
        void* info;
        int index;
        Distribution() : manager(0), counts(0), info(0), index(-1) {}
        
        bool ownMem() const {
            return !manager || (counts<&manager->globalDistCount[0]) || (counts>&manager->globalDistCount.back());
        }
        
        protected:
        void internal_copy(const Distribution& d) {
            manager = d.manager;
            info = d.info;
            index = d.index;
            // shared global vector if possible
            if (!d.ownMem()) counts = d.counts;
            else if (manager) {
                // dup otherwise
                counts = new int[manager->nz];
                for (int i=0; i<manager->nz; ++i) counts[i] = d.counts[i];
            } else counts = 0;
        }
        public:
        
        // copy constructor and operator=
        Distribution(const Distribution& d) {
            internal_copy(d);
        }
        Distribution& operator=(const Distribution& d) {
            if (ownMem()) delete[] counts;
            internal_copy(d);
            return *this;
        }
        
        ~Distribution() {
            if (ownMem()) { delete[] counts; counts = 0; }
        }
        
        // Allows clustering algorithm to fusion distributions
        void aggregate(const Distribution& d) {
            if (counts==0 && d.counts==0) return; // nothing to aggregate
            if (manager==0 && d.manager==0) return; // nothing to aggregate
            if (manager==0) {
                internal_copy(d);
                return;
            }
            if (d.manager==0 && d.counts==0) return;
            assert(manager==d.manager);
            // detach from global vector if not already the case
            if (counts==0 || !ownMem()) {
                int* ncounts = new int[manager->nz];
                for (int i=0; i<manager->nz; ++i) ncounts[i] = counts[i];
                counts = ncounts;
            }
            for (int i=0; i<manager->nz; ++i) counts[i] += d.counts[i];
        }
    };
    
    // map => unicity, and the int argument is the ordered index in map order
    boost::unordered_map<DataType,int> uniqueX;
    boost::unordered_map<PredictionType,int> uniqueZ;
    std::vector<int> globalDistCount;
    std::vector<Distribution> distribs;
    int nx, nz;

    size_t size() {return nx;}
    typedef typename std::vector<Distribution>::iterator iterator;
    iterator begin() { return distribs.begin(); }
    iterator end() { return distribs.end(); }
    // iterator dereferences to a Distribution
    typedef Distribution value_type;
    
    DiscreteDistributionManager(DataSet& dataset) {
        // First pass: get unique X and Z in the data set
        for (typename DataSet::iterator datait = dataset.begin(); datait != dataset.end(); ++datait) {
            uniqueX[helpers::DataSetIteratorUnrefHelper<DataSet,TypeTraits>::getData(dataset,datait)] = 0;
            uniqueZ[helpers::DataSetIteratorUnrefHelper<DataSet,TypeTraits>::getPrediction(dataset,datait)] = 0;
        }
        uniqueX.max_load_factor(1.0f);
        uniqueZ.max_load_factor(1.0f);
        uniqueX.rehash(uniqueX.bucket_count());
        uniqueZ.rehash(uniqueZ.bucket_count());
        
        // compute indices
        nx = 0;
        for (typename boost::unordered_map<DataType,int>::iterator it = uniqueX.begin(); it != uniqueX.end(); ++it) {
            it->second = nx++;
        }
        nz = 0;
        for (typename boost::unordered_map<PredictionType,int>::iterator it = uniqueZ.begin(); it != uniqueZ.end(); ++it) {
            it->second = nz++;
        }

        // now work only on the indices
        // memory-efficient distributions for non-sparse data = all in one big vector
        // so as to minimise extra variables
        globalDistCount.resize(nx * nz, 0);
        
        // Second pass: accumulate distributions in the global vector
        for (typename DataSet::iterator datait = dataset.begin(); datait != dataset.end(); ++datait) {
            int xidx = uniqueX[helpers::DataSetIteratorUnrefHelper<DataSet,TypeTraits>::getData(dataset,datait)];
            int zidx = uniqueZ[helpers::DataSetIteratorUnrefHelper<DataSet,TypeTraits>::getPrediction(dataset,datait)];
            ++globalDistCount[xidx * nz + zidx];
        }
        
        // distributions are just placeholders for extra info
        distribs.resize(nx);
        for (int i=0; i<nx; ++i) {
            distribs[i].manager = this;
            distribs[i].counts = &globalDistCount[i*nz];
            distribs[i].info = 0;
            distribs[i].index = i;
        }
    }
    Distribution nulldist;

    Distribution& distribution(const DataType& x) {
        typename boost::unordered_map<DataType,int>::iterator xit = uniqueX.find(x);
        if (xit != uniqueX.end()) return distribs[xit->second];
        return nulldist;
    }
    Distribution& operator[](int index) {
        return distribs[index];
    }
    
    const DataType* getData(int index) {
        for (typename boost::unordered_map<DataType,int>::iterator it = uniqueX.begin(); it != uniqueX.end(); ++it) {
            if (index==it->second) return &it->first;
        }
        return 0;
    }
    const DataType* getData(const Distribution& dist) {
        return getData(dist.index);
    }

    void* getInfo(const Distribution& dist) {
        return dist.info;
    }
    void setInfo(Distribution& dist, void* info) {
        dist.info = info;
    }

    int getCount(const Distribution& dist) {
        if (dist.manager==0) return 0;
        int ret = 0;
        for (int i=0; i<dist.manager->nz; ++i) ret += dist.counts[i];
        return ret;
    }

    template<class Functor>
    FloatType expectation(const Distribution& dist, Functor f) {
        if (dist.manager!=this) return 0;
        FloatType result = 0;
        int totalcounts = 0;
        for (typename boost::unordered_map<PredictionType,int>::iterator it = uniqueZ.begin(); it != uniqueZ.end(); ++it) {
            result += f(it->first) * dist.counts[it->second];
            totalcounts += dist.counts[it->second];
        }
        return result / totalcounts;
    }

    // Test if the 2 distributions match, perform a Chi-square computation.
    // Return true if the chi-square result is above the given required level
        
    struct ChiSquareMatcher {
        float matchLevel;
        ChiSquareMatcher(float m = 0.95) : matchLevel(m) {}

        bool operator()(const Distribution& a, const Distribution& b) {
            if (a.manager==0) return 0;
            if (a.manager!=b.manager) return 0;
            FloatType n1 = a.manager->getCount(a);
            FloatType n2 = b.manager->getCount(b);
            FloatType n1n2 = n1 * n2;
            FloatType n1overn2 = n1 / n2;
            FloatType n2overn1 = n2 / n1;
            FloatType chisq = 0.0f;
            unsigned int ndegree = 0;

            for (int zidx = 0; zidx < a.manager->nz; ++zidx) {
                if (a.counts[zidx]==0 && b.counts[zidx]==0) continue;
                ++ndegree;
                if (b.counts[zidx]==0) {
                    chisq += n2overn1 * a.counts[zidx];
                    ++ndegree;
                }
                else if (a.counts[zidx]==0) {
                    chisq += n1overn2 * b.counts[zidx];
                    ++ndegree;
                }
                else {
                    FloatType tmp = n2 * a.counts[zidx] - n1 * b.counts[zidx];
                    chisq += tmp * tmp / (n1n2 * (a.counts[zidx] + b.counts[zidx]));
                    ++ndegree;
                }
            }
            
            return helpers::ChiSquareMatch(ndegree, chisq, matchLevel);
        }

    };
    
    // For default arguments
    typedef ChiSquareMatcher DistributionMatcher;
    
};


}

#endif
