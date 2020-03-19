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

#ifndef DECISIONAL_STATES_JOINT_DATA_SET_H
#define DECISIONAL_STATES_JOINT_DATA_SET_H

#include "decisional_states_helpers.hpp"

#include <limits>
#include <vector>
#include <map>
#include <cstring>

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace decisional_states {

namespace helpers {
template<typename FloatType>
struct IndexedContribution {
    int index;
    FloatType contribution;
    IndexedContribution(int _index, FloatType _contribution) : index(_index), contribution(_contribution) {}
    IndexedContribution() : index(-1), contribution(-1) {}
};
template<typename FloatType>
struct IndexedContributionVector {
    boost::shared_array<IndexedContribution<FloatType> > contributions;
    int size;
    IndexedContributionVector clone() const {
        IndexedContributionVector cloned;
        cloned.size = size;
        cloned.contributions = boost::shared_array<IndexedContribution<FloatType> >(new IndexedContribution<FloatType>[size]);
        std::memcpy(&cloned.contributions[0], &contributions[0], size*sizeof(IndexedContribution<FloatType>));
        return cloned;
    }
};
}

/**
Base class helper for making data sets and distribution manager, where the data and prediction types are derived from a "joint" type that is more natural, on which we can apply the random Markov field theory.
Ex: prediction = pixel, data = neighbourhood of that pixel in an image (minus center pixel) => joint type = whole zone (ex: 3x3)
Ex2: data = past 5 points in a time series, prediction = next 3 points => joint type = the 8 points, now easier to define on a sliding window

There are two ways of using this class:

Push mode:
- Build an object, then add() as many joint observations as desired. The object is ready to be used 
*/

/**
    Contiguous range of non-null indices
    Might save memory compared to storing the indices one by one, especially when applying kernels
    based on the indexing space => the range of non-null values is well localised
*/
template<typename FloatType>
struct IndexRangeDistributionStorage {

    // indexes into the prediction sampler
    int minNonNullContrib, maxNonNullContrib;
    // array of size maxNonNullContrib - minNonNullContrib + 1, when min <= max, or null pointer otherwise
    FloatType* contribs;

    IndexRangeDistributionStorage()
    : minNonNullContrib(std::numeric_limits<int>::max()), maxNonNullContrib(-1), contribs(0) {
    }
    ~IndexRangeDistributionStorage() {
        if (contribs) delete [] contribs;
        contribs = 0;
    }
    IndexRangeDistributionStorage& operator=(const IndexRangeDistributionStorage& other) {
        minNonNullContrib = other.minNonNullContrib;
        maxNonNullContrib = other.maxNonNullContrib;
        if (contribs) delete [] contribs;
        if (other.contribs==0) contribs = 0;
        else {
            int ncontribs = maxNonNullContrib - minNonNullContrib + 1;
            contribs = new FloatType[ncontribs];
            for(int i=0; i<ncontribs; ++i) contribs[i] = other.contribs[i];
        }
        return *this;
    }
    IndexRangeDistributionStorage(const IndexRangeDistributionStorage& other) : contribs(0) {*this = other;}

    // store the indices and their associated contributions in the internal format
    template<class Container>
    void store(Container& contributions) {
        if (contribs) delete [] contribs;
        minNonNullContrib = 0;
        maxNonNullContrib = contributions.size()-1;
        while(minNonNullContrib<(int)contributions.size() && (contributions[minNonNullContrib]==0)) ++minNonNullContrib;
        while(maxNonNullContrib>=0 && (contributions[maxNonNullContrib]==0)) --maxNonNullContrib;
        if (minNonNullContrib <= maxNonNullContrib) {
            int ncontribs = maxNonNullContrib - minNonNullContrib + 1;
            contribs = new FloatType[ncontribs];
            for (int i=0; i<ncontribs; ++i) contribs[i] = contributions[i+minNonNullContrib];
        }
    }

    void store(const helpers::IndexedContributionVector<FloatType>& icv) {
        if (contribs) delete [] contribs;
        // The icv are not sorted by index => must look for min/max
        minNonNullContrib = std::numeric_limits<int>::max();
        maxNonNullContrib = std::numeric_limits<int>::min();
        for (int i = 0; i < icv.size; ++i) {
            int idx = icv.contributions[i].index;
            minNonNullContrib = std::min(idx, minNonNullContrib);
            maxNonNullContrib = std::max(idx, maxNonNullContrib);
        }
        if (minNonNullContrib <= maxNonNullContrib) {
            int ncontribs = maxNonNullContrib - minNonNullContrib + 1;
            contribs = new FloatType[ncontribs];
            for (int i=0; i<ncontribs; ++i) contribs[i] = 0;
            for (int i=0; i<icv.size; ++i) contribs[icv.contributions[i].index - minNonNullContrib] = icv.contributions[i].contribution;
        }
    }

    int size() const {
        if (contribs==0) return 0;
        if (maxNonNullContrib < minNonNullContrib) return 0;
        return maxNonNullContrib - minNonNullContrib + 1;
    }
    bool empty() const { return size()==0; }

    // pseudo-container with accessors
    typedef int iterator;
    iterator begin() {return 0;}
    iterator end() {return size();}
    // index in the original container
    /// Indices MUST be sorted during iterator traversal
    int index(iterator it) {return it + minNonNullContrib;}
    FloatType contribution(iterator it) {return contribs[it];}

    // Aggregate directly at the storage level, called by the distributions
    void aggregate(const IndexRangeDistributionStorage& other) {
        if (empty()) {
            *this = other;
            return;
        }
        if (other.empty()) return; // nothing to aggregate

        int newmin = min(minNonNullContrib, other.minNonNullContrib);
        int newmax = max(maxNonNullContrib, other.maxNonNullContrib);
        assert(newmin<=newmax); // cannot happen if both are non-empty

        FloatType* newcontribs = new FloatType[newmax - newmin + 1];
        for (int i=newmin; i<=newmax; ++i) newcontribs[i-newmin]=0;
        for (int i=minNonNullContrib; i<=maxNonNullContrib; ++i) newcontribs[i-newmin]+=contribs[i-minNonNullContrib];
        for (int i=other.minNonNullContrib; i<=other.maxNonNullContrib; ++i) newcontribs[i-newmin]+=other.contribs[i-other.minNonNullContrib];
        delete[] contribs;
        contribs = newcontribs;
        minNonNullContrib = newmin;
        maxNonNullContrib = newmax;
    }


    // TODO: better handling of this case, directly aggregating with internal storage format
    //       yet the above code is not trivial, and copy/paste is not a good idea for maintainance
    void aggregate(const helpers::IndexedContributionVector<FloatType>& icv) {
        IndexRangeDistributionStorage other;
        other.store(icv);
        aggregate(other);
    }

};

/// Store the non-null indices explicitly, unlike the range
/// This however takes (possibly much) less memory when the indices of non-null contributions are non-contiguous
/// Ex: indices taken from any kind of randomized sampler
/// For contiguous ranges... use the range storage.
template<typename FloatType>
struct IndexedDistributionStorage {

    int nelements;
    int* indices;
    FloatType* contribs;

    IndexedDistributionStorage() : nelements(0), indices(0), contribs(0) {}

    ~IndexedDistributionStorage() {
        if (contribs) delete [] contribs;
        if (indices) delete [] indices;
        contribs = 0;
        indices = 0;
        nelements = 0;
    }
    IndexedDistributionStorage& operator=(const IndexedDistributionStorage& other) {
        nelements = other.nelements;
        if (contribs) delete [] contribs;
        if (indices) delete [] indices;
        if (nelements==0) {
            contribs = 0;
            indices = 0;
        }
        else {
            contribs = new FloatType[nelements];
            indices = new int[nelements];
            for(int i=0; i<nelements; ++i) {
                contribs[i] = other.contribs[i];
                indices[i] = other.indices[i];
            }
        }
        return *this;
    }
    IndexedDistributionStorage(const IndexedDistributionStorage& other) : indices(0), contribs(0) {*this = other;}

    /// store the indices and their associated contributions in the internal format
    template<class Container>
    void store(Container& contributions) {
        if (contribs) delete [] contribs;
        if (indices) delete [] indices;
        nelements = 0;
        for (int i=0; i<(int)contributions.size(); ++i) if (contributions[i]!=0) ++nelements;
        if (nelements==0) {
            contribs = 0;
            indices = 0;
            return;
        }
        contribs = new FloatType[nelements];
        indices = new int[nelements];
        nelements = 0;
        for (int i=0; i<(int)contributions.size(); ++i) if (contributions[i]!=0) {
            indices[nelements] = i;
            contribs[nelements] = contributions[i];
            ++nelements;
        }
    }

    void store(const helpers::IndexedContributionVector<FloatType>& icv) {
        if (contribs) delete [] contribs;
        if (indices) delete [] indices;
        nelements = icv.size;
        contribs = new FloatType[icv.size];
        indices = new int[icv.size];
        for (int i=0; i<icv.size; ++i) {
            indices[i] = icv.contributions[i].index;
            contribs[i] = icv.contributions[i].contribution;
        }
    }

    int size() const {return nelements;}
    bool empty() const { return nelements==0; }

    // pseudo-container with accessors
    typedef int iterator;
    iterator begin() {return 0;}
    iterator end() {return nelements;}
    // index in the original container
    /// Indices MUST be sorted during iterator traversal
    int index(iterator it) {return indices[it];}
    FloatType contribution(iterator it) {return contribs[it];}

    // Aggregate directly at the storage level, called by the distributions
    void aggregate(const IndexedDistributionStorage& other) {
        if (nelements==0) {(*this) = other; return;}
        if (other.nelements==0) return;
        // map preserves the indices order
        std::map<int, FloatType> mergedContribs;
        // merge all entries in the map, knowing map elements are default-initialized hence 0
        for (int i=0; i<nelements; ++i) mergedContribs[indices[i]] += contribs[i];
        for (int i=0; i<other.nelements; ++i) mergedContribs[other.indices[i]] += other.contribs[i];
        // copy the map back to our own arrays
        nelements = mergedContribs.size();
        assert(nelements >0); // since nullity case was tested above
        if (contribs) delete [] contribs;
        contribs = new FloatType[nelements];
        if (indices) delete [] indices;
        indices = new int[nelements];
        int i=0;
        for (typename std::map<int,FloatType>::iterator it = mergedContribs.begin(); it != mergedContribs.end(); ++it, ++i) {
            indices[i] = it->first;
            contribs[i] = it->second;
        }
    }

    // aggregate, using the internal storage format
    void aggregate(const helpers::IndexedContributionVector<FloatType>& icv) {
        if (nelements==0) {store(icv); return;}
        if (icv.size==0) return;
        // map preserves the indices order
        std::map<int, FloatType> mergedContribs;
        // merge all entries in the map, knowing map elements are default-initialized hence 0
        for (int i=0; i<nelements; ++i) mergedContribs[indices[i]] += contribs[i];
        for (int i=0; i<icv.size; ++i) mergedContribs[icv.contributions[i].index] += icv.contributions[i].contribution;
        // copy the map back to our own arrays
        nelements = mergedContribs.size();
        assert(nelements >0); // since nullity case was tested above
        if (contribs) delete [] contribs;
        contribs = new FloatType[nelements];
        if (indices) delete [] indices;
        indices = new int[nelements];
        int i=0;
        for (typename std::map<int,FloatType>::iterator it = mergedContribs.begin(); it != mergedContribs.end(); ++it, ++i) {
            indices[i] = it->first;
            contribs[i] = it->second;
        }
    }
};


/// Container view to get the data and predictions (DataSet API)
template<class TypeTraits>
struct DistDataSet {
    typedef typename TypeTraits::FloatType FloatType;
    typedef typename TypeTraits::DataType DataType;
    typedef typename TypeTraits::PredictionSampler PredictionSampler;
    typedef typename TypeTraits::PredictionType PredictionType;

    typedef typename helpers::IndexedContribution<FloatType> IndexedContribution;
    typedef typename helpers::IndexedContributionVector<FloatType> IndexedContributionVector;

    struct DistributionBase : public TypeTraits::DistributionStorage {
        // number of aggregated data from the data set
        int count;
        void* info;
        FloatType sumContribs;
        DistributionBase() : count(0), info(0), sumContribs(0) {}
    };

    struct Distribution : public boost::shared_ptr<DistributionBase> {

        Distribution() {}
        
        Distribution clone() {
            Distribution ret;
            *static_cast<boost::shared_ptr<DistributionBase>*>(&ret) = boost::shared_ptr<DistributionBase>(new DistributionBase(*this->get()));
            return ret;
        }

        Distribution(const IndexedContributionVector& icv) : boost::shared_ptr<DistributionBase>(new DistributionBase()) {
            (*this)->store(icv);
            (*this)->count = 1; // corresponds to one observed data by default
            (*this)->sumContribs = 0;
            for (int i=0; i<icv.size; ++i) (*this)->sumContribs += icv.contributions[i].contribution;
        }

        template<class Container>
        Distribution(const Container& contributions)
        : boost::shared_ptr<DistributionBase>(new DistributionBase()) {
            (*this)->store(contributions);
            (*this)->count = 1; // corresponds to one observed data by default
            (*this)->sumContribs = 0;
            for(typename DistributionBase::iterator it = (*this)->begin(); it != (*this)->end(); ++it) {
                (*this)->sumContribs += (*this)->contribution(it);
            }
        }

        void aggregate(const Distribution& other) {
            // MUST make a copy, otherwise further aggregate will modify both this and the original
            if (this->get()==0) (*(boost::shared_ptr<DistributionBase>*)this) = boost::shared_ptr<DistributionBase>(new DistributionBase());
            if (other.get()==0 || other->empty()) return;
            (*this)->aggregate(*other);
            (*this)->sumContribs += other->sumContribs;
            (*this)->count += other->count;
        }

        void aggregate(const IndexedContributionVector& icv) {
            if (this->get()==0) (*(boost::shared_ptr<DistributionBase>*)this) = boost::shared_ptr<DistributionBase>(new DistributionBase());
            if (icv.size==0) return;
            (*this)->aggregate(icv);
            for (int i=0; i<icv.size; ++i) (*this)->sumContribs += icv.contributions[i].contribution;
        }

        int getMaxContribIndex() {
            FloatType maxc = -std::numeric_limits<float>::max();
            int maxidx = 0;
            for(typename DistributionBase::iterator it = (*this)->begin(); it != (*this)->end(); ++it) {
                if ( (*this)->contribution(it) > maxc ) {
                    maxc = (*this)->contribution(it);
                    maxidx = (*this)->index(it);;
                }
            }
            return maxidx;
        }

        typedef typename DistributionBase::iterator iterator;
        iterator begin() {return (*this)->begin();}
        iterator end() {return (*this)->end();}
        int size() {return (*this)->size();}
    };

    // x => distrib, aggregated for similar x, just in case, and also to keep all known x
    struct DataDistPair {
        DataType data;
        mutable Distribution distribution;
        DataDistPair(DataType _data, Distribution _distribution) : data(_data), distribution(_distribution) {}
    };

    typedef boost::multi_index_container<
        DataDistPair,
        boost::multi_index::indexed_by<
            boost::multi_index::hashed_unique< boost::multi_index::member<DataDistPair,DataType,&DataDistPair::data> >
            //,boost::multi_index::random_access<>
        >
    > DistMap;
    DistMap distmap;

    // pointers to elements in the map are never invalidated (dixit boost doc)
    std::vector<DataType*> dataSequence;

    // The PredictionSampler provides the indices common to all distributions
    PredictionSampler predictionSampler;

    DistDataSet(PredictionSampler _predictionSampler) : predictionSampler(_predictionSampler) {}

    // no sequence, just aggregation
    void aggregate(const DataType& data, const IndexedContributionVector& icv) {
        typename DistMap::iterator it = distmap.find(data);
        if (it==distmap.end()) {
            distmap.insert(DataDistPair(data,Distribution(icv)));
            return;
        }
        it->distribution.aggregate(icv);
    }

    std::vector<DataType> missingData;
    
    void prepareSequence(int seqsize) {
        dataSequence.resize(seqsize);
    }
    
    // add a full distribution, maintain sequence numbers
    void addDist(const DataType& x, Distribution d, int seqnum=-1) {
        // simulate push_back
        if (seqnum==-1) {seqnum = dataSequence.size(); dataSequence.resize(seqnum+1);}
        // missing data are kept in seq, but the null distribution is not added
        if (helpers::IsMissingHelper<typename TypeTraits::JointType,DataType>::isMissing(x) || d.get()==0 || d->sumContribs==0) {
            missingData.push_back(x);
            dataSequence[seqnum] = &missingData.back();
            return;
        }
        typename DistMap::iterator it = distmap.find(x);
        if (it==distmap.end()) it = distmap.insert(DataDistPair(x,d.clone())).first;
        else it->distribution.aggregate(d);
        dataSequence[seqnum] = const_cast<DataType*>(&it->data);
    }

    typedef typename std::vector<DataType*>::iterator iterator;

    DataType& data(iterator it) {
        return **it;
    }
    
    DataType& data(int index) {
        return *dataSequence[index];
    }

    // just to complete the API, might not be the original prediction in case of aggregation of conflicting predictions for the same data
    PredictionType& prediction(iterator it) {
        DataType& x = **it;
        typename DistMap::iterator dmit = distmap.find(x);
        //if (it==distmap.end()) return nullPrediction; // shall not happen, actually assert it 
        assert(dmit!=distmap.end());  // use isMissing before calling this function
        int idx = dmit->distribution.getMaxContribIndex();
        return predictionSampler[idx];
    }

    iterator begin() {return dataSequence.begin();}
    iterator end() {return dataSequence.end();}
    std::size_t size() {return dataSequence.size();}

    static bool isMissing(DataType& d) {
        return helpers::IsMissingHelper<typename TypeTraits::JointType,DataType>::isMissing(d);
    }

    static bool isMissing(PredictionType& p) {
        return helpers::IsMissingHelper<typename TypeTraits::JointType,DataType>::isMissing(p);
    }

};

/// Container view to get the distributions (DistributionManager API)

template<class TypeTraits>
struct BaseDistributionManager : public DistDataSet<TypeTraits> {
    typedef typename DistDataSet<TypeTraits>::Distribution Distribution;
    typedef typename DistDataSet<TypeTraits>::DataDistPair DataDistPair;
    typedef typename DistDataSet<TypeTraits>::DistMap DistMap;
    typedef typename DistDataSet<TypeTraits>::IndexedContribution IndexedContribution;
    typedef typename DistDataSet<TypeTraits>::IndexedContributionVector IndexedContributionVector;
    typedef typename TypeTraits::FloatType FloatType;
    typedef typename TypeTraits::DataType DataType;
    typedef Distribution value_type;
    Distribution nullDist;

    struct DistMapToDist {
        Distribution& operator()(const DataDistPair& p) const {
            return p.distribution;
        }
        typedef Distribution& result_type;
    };

    // nth_index<0> with random access
    typedef helpers::ContainerAdaptor<typename DistMap::template nth_index<0>::type, Distribution, DistMapToDist> DistAdaptor;
    DistAdaptor adaptedMap;


    BaseDistributionManager(typename TypeTraits::PredictionSampler _predictionSampler) : DistDataSet<TypeTraits>(_predictionSampler), adaptedMap(static_cast<DistDataSet<TypeTraits>*>(this)->distmap.template get<0>()) {}
    // get<1> with random access

    typedef typename DistAdaptor::iterator iterator;
    iterator begin() {return adaptedMap.begin();}
    iterator end() {return adaptedMap.end();}
    std::size_t size() {return this->dataSequence.size();}

    Distribution& operator[](int index) {
        typename DistMap::iterator it = static_cast<DistDataSet<TypeTraits>*>(this)->distmap.find(*this->dataSequence[index]);
        if (it==static_cast<DistDataSet<TypeTraits>*>(this)->distmap.end()) return nullDist;
        return it->distribution;
//        return adaptedMap[index];
    }

    Distribution& distribution(const DataType& x) {
        if (helpers::IsMissingHelper<typename TypeTraits::JointType,DataType>::isMissing(x)) {
            nullDist = Distribution();
            return nullDist;
        }
        typename DistMap::iterator it = static_cast<DistDataSet<TypeTraits>*>(this)->distmap.find(x);
        if (it==static_cast<DistDataSet<TypeTraits>*>(this)->distmap.end()) return nullDist;
        return static_cast<DistDataSet<TypeTraits>*>(this)->distmap.find(x)->distribution;
    }

    bool isNull(Distribution& dist) {
        return dist.get()==0 || dist->empty() || dist->sumContribs==0;
    }

    void* getInfo(Distribution& dist) {
        if (dist.get()==0 || dist->empty()) return 0;
        return dist->info;
    }
    void setInfo(Distribution& dist, void* info) {
        if (dist.get()==0 || dist->empty()) return;
        dist->info = info;
    }
    int getCount(Distribution& dist) {
        if (dist.get()==0 || dist->empty()) return 0;
        return dist->count;
    }

    template<class Functor>
    FloatType expectation(const Distribution& dist, Functor f) {
        if (dist.get()==0 || dist->empty()) return 0;
        FloatType ret = 0;
        for(typename Distribution::iterator it = dist->begin(); it != dist->end(); ++it) {
            ret += f(static_cast<DistDataSet<TypeTraits>*>(this)->predictionSampler[dist->index(it)]) * dist->contribution(it);
        }
        return ret / dist->sumContribs;
    }

    struct BhattacharyyaMatcher {
        FloatType matchlevel;

        BhattacharyyaMatcher(FloatType _matchlevel = 0.95) : matchlevel(_matchlevel * 2 - _matchlevel * _matchlevel) {}

        FloatType operator()(const Distribution& a, const Distribution& b) const {
            if (a.get()==0 || b.get()==0) return -1;
            if (a->empty() || b->empty()) return -1;

            FloatType bhat = 0;
            typename Distribution::iterator ait = a->begin();
            typename Distribution::iterator aend = a->end();
            typename Distribution::iterator bit = b->begin();
            typename Distribution::iterator bend = b->end();

            // indices are sorted by contract of the distribution storage
            // so we can advance iterators on missing indices
            while(true) {
                int aidx = a->index(ait);
                int bidx = b->index(bit);
                if (aidx < bidx) {
                    ++ait; if (ait==aend) break;
                    continue;
                }
                if (bidx < aidx) {
                    ++bit; if (bit==bend) break;
                    continue;
                }
                // both match, we can add non-null sqrt
                bhat += helpers::sqrt(a->contribution(ait) * b->contribution(bit));
                ++ait; if (ait==aend) break;
                ++bit; if (bit==bend) break;
            }
            // if more elements remain, ignore and do not add 0
            // return >0 is match
            // b = sum_i sqrt((pza_i/sumca) * (pzb_i/sumcb))
            // b = sum_i sqrt(pza_i * pzb_i/(sumca*sumcb))
            // b = sum_i sqrt(pza_i * pzb_i)/sqrt(sumca*sumcb)
            // b = (sum_i sqrt(pza_i * pzb_i) )/sqrt(sumca*sumcb)
            // b > matchlevel => sum_i sqrt(pza_i * pzb_i) > matchlevel * sqrt(sumca*sumcb)
            return bhat - matchlevel * helpers::sqrt(a->sumContribs * b->sumContribs);
        }
    };
    typedef BhattacharyyaMatcher DistributionMatcher;

};

/// Main object, joint API, that can be viewed as either of the above

template<class TypeTraits>
struct JointDataManager : public BaseDistributionManager<TypeTraits> {

    typedef typename TypeTraits::FloatType FloatType;
    typedef typename TypeTraits::JointType JointType;
    typedef typename TypeTraits::DataType DataType;
    typedef typename TypeTraits::PredictionType PredictionType;
    typedef typename TypeTraits::PredictionSampler PredictionSampler;
    typedef typename TypeTraits::JointKernel JointKernel;
    typedef typename TypeTraits::JointNeighbourhoodFinder JointNeighbourhoodFinder;
    typedef typename TypeTraits::DataNeighbourhoodFinder DataNeighbourhoodFinder;
    typedef typename TypeTraits::PredictionNeighbourhoodFinder PredictionNeighbourhoodFinder;

    typedef typename DistDataSet<TypeTraits>::Distribution Distribution;
    typedef typename DistDataSet<TypeTraits>::IndexedContribution IndexedContribution;
    typedef typename DistDataSet<TypeTraits>::IndexedContributionVector IndexedContributionVector;

    /// Derived types
    typedef DistDataSet<TypeTraits> DataSet;
    typedef BaseDistributionManager<TypeTraits> DistributionManager;

    DataSet& asDataSet() {return *static_cast<DataSet*>(this);}
    DistributionManager& asDistributionManager() {return *static_cast<DistributionManager*>(this);}

    JointDataManager(PredictionSampler _predictionSampler) : DistributionManager(_predictionSampler) {}

    template<class JointType, class SpacialJointSampler>
    struct SumContribAction {
        JointType& workjt;
        SpacialJointSampler& sampler;
        JointKernel& kernel;
        FloatType& contrib;
        SumContribAction(JointType& _workjt, SpacialJointSampler& _sampler, JointKernel& _kernel, FloatType& _contrib) : workjt(_workjt), sampler(_sampler), kernel(_kernel), contrib(_contrib) {}
        template<typename Iterator>
        void operator()(std::size_t index, const Iterator& begin, const Iterator& end) {
            for (Iterator it = begin; it != end; ++it) {
                // single element container wrapper => ignore index == 0
                contrib += kernel(workjt, sampler[*it]);
            }
        }
    };

    template<class PredictionKernel>
    struct ComputePredictionFactorAction {
        IndexedContributionVector& icv;
        const PredictionType& prediction;
        PredictionKernel& predictionKernel;
        PredictionSampler& predictionSampler;
        ComputePredictionFactorAction(IndexedContributionVector& _icv, const PredictionType& _prediction, PredictionKernel& _predictionKernel, PredictionSampler& _predictionSampler) : icv(_icv), prediction(_prediction), predictionKernel(_predictionKernel), predictionSampler(_predictionSampler) {}
        template<typename Iterator>
        void operator()(std::size_t index, const Iterator& begin, const Iterator& end) {
            icv.size = end - begin;
            icv.contributions = boost::shared_array<IndexedContribution>(new IndexedContribution[icv.size]);
            for (Iterator it = begin; it != end; ++it) {
                icv.contributions[it-begin].index = *it;
                icv.contributions[it-begin].contribution = predictionKernel(prediction, predictionSampler[*it]);
            }
        }
    };

    template<class DataKernel, class DataSampler>
    struct ComputeDataFactorAction {
        IndexedContributionVector& icv;
        const DataType& data;
        DataKernel& dataKernel;
        const DataSampler& dataSampler;
        ComputeDataFactorAction(IndexedContributionVector& _icv, const DataType& _data, DataKernel& _dataKernel, const DataSampler& _dataSampler) : icv(_icv), data(_data), dataKernel(_dataKernel), dataSampler(_dataSampler) {}
        template<typename Iterator>
        void operator()(std::size_t index, const Iterator& begin, const Iterator& end) {
            icv.size = end - begin;
            icv.contributions = boost::shared_array<IndexedContribution>(new IndexedContribution[icv.size]);
            for (Iterator it = begin; it != end; ++it) {
                icv.contributions[it-begin].index = *it;
                icv.contributions[it-begin].contribution = dataKernel(data, dataSampler[*it]);
            }
        }
    };

    struct JointToDataConverter {
        const DataType& operator()(const JointType& j) const {return j.getData();}
        typedef const DataType& result_type;
    };


    template<class JointContainer, class DataKernel, class PredictionKernel, class JointSampler>
    void addSeparable(JointType jt, DataKernel dataKernel, PredictionKernel predictionKernel, JointSampler& jsampler) {
        addSeparable(helpers::SingleElementContainerWrapper<JointType>(jt), dataKernel, predictionKernel, jsampler);
    }

    template<class JointContainer, class DataKernel, class PredictionKernel, class JointSampler>
    void addSeparable(JointContainer jointContainer, DataKernel dataKernel, PredictionKernel predictionKernel, JointSampler& jsampler) {
        // syntactic alias for clarity
        PredictionSampler& predictionSampler = static_cast<DistDataSet<TypeTraits>*>(this)->predictionSampler;
        PredictionNeighbourhoodFinder pnf;
        // neighbour finders ignore NaN values
        pnf.setup(predictionSampler);
        // TODO: boolean flags for trying to gather similar predictions or not
        //       This is very good for small discrete space, as it avoids recomputing the neighbourhoods for the same predictions
        //       when many joint samples projected on prediction space give the same predictions.
        //       This is quite bad for large/continuous spaces where the probability of having exactly the same prediction tends to 0
        //       and we're just wasting memory, possibly even too much to fit in physical ram.
        boost::unordered_map<PredictionType, IndexedContributionVector> uniquePredictionContributions;
        for (int i=0; i<(int)jsampler.size(); ++i) {
            const PredictionType& p = jsampler[i].getPrediction();
            if (helpers::IsMissingHelper<JointType,PredictionType>::isMissing(p)) continue;
            // do computations only once
            if (uniquePredictionContributions.find(p) != uniquePredictionContributions.end()) continue;
            IndexedContributionVector icv;
            pnf.find(
                helpers::SingleElementContainerWrapper<PredictionType>(p),
                predictionKernel.neighbourhoodSize(),
                // icv (and other args) passed by ref
                ComputePredictionFactorAction<PredictionKernel>(icv, p, predictionKernel, predictionSampler)
            );
            // icv was updated, use it
            uniquePredictionContributions.insert(std::make_pair(p,icv));
        }

        DataNeighbourhoodFinder dnf;
        helpers::ContainerAdaptor<JointContainer, DataType, JointToDataConverter> adaptedSamplesContainer(jsampler);
        dnf.setup(adaptedSamplesContainer);

        int jcsize = jointContainer.size();
        this->prepareSequence(jcsize);
        
        // Implementation note: openmp was observed to give better results in the present case
        // at this outer loop level, despite the critical sync.
        // An inner parallel search seems too costly in terms of thread setup
        
#if defined(_OPENMP)
#pragma omp parallel for
#endif
        for (int jci = 0; jci < jcsize; ++jci) {
            const DataType& currentData = jointContainer[jci].getData();

            if (helpers::IsMissingHelper<JointType,DataType>::isMissing(jointContainer[jci],currentData)) {
#if defined(_OPENMP)
#pragma omp critical
#endif
{
                addDist(currentData, Distribution(), jci);
}
            } else {

            // find neighbours in samplers based only on data part
            IndexedContributionVector sampleDataContribs; // index in ssampler, data kernel contrib
            dnf.find(
                helpers::SingleElementContainerWrapper<DataType>(currentData),
                dataKernel.neighbourhoodSize(),
                ComputeDataFactorAction<DataKernel, helpers::ContainerAdaptor<JointContainer, DataType, JointToDataConverter> >(sampleDataContribs, currentData, dataKernel, adaptedSamplesContainer)
            );

            std::vector<FloatType> contributions(predictionSampler.size(), FloatType(0));
            for (int i=0; i<sampleDataContribs.size; ++i)
            {
                IndexedContribution& ic = sampleDataContribs.contributions[i];
                const PredictionType& p = jsampler[ic.index].getPrediction();
                if (helpers::IsMissingHelper<JointType,PredictionType>::isMissing(p)) continue;
                FloatType dataContrib = ic.contribution;
                const IndexedContributionVector& icv = uniquePredictionContributions[p];
                for (int j=0; j<icv.size; ++j) contributions[icv.contributions[j].index] += icv.contributions[j].contribution * dataContrib;
            }
#if defined(_OPENMP)
#pragma omp critical
#endif
{
            Distribution currentDistribution(contributions);
            addDist(currentData,currentDistribution,jci);
}

            }
        }
    }


    template<class SpacialJointSampler>
    void add(JointType jt, JointKernel jointKernel, SpacialJointSampler& ssampler) {
        add(helpers::SingleElementContainerWrapper<JointType>(jt), jointKernel, ssampler);
    }

    // Non-separable kernel case
    // assume the space of JointType is so large that it is not efficient to gather similar values + count memory-wise (this is only useful if count > 1 on avg)
    template<class JointContainer, class SpacialJointSampler>
    void add(JointContainer jointContainer, JointKernel jointKernel, SpacialJointSampler& ssampler) {
        JointNeighbourhoodFinder jnf;
        jnf.setup(ssampler);
        std::vector<FloatType> contributions(static_cast<DistDataSet<TypeTraits>*>(this)->predictionSampler.size(), FloatType(0));
        for (typename JointContainer::iterator it = jointContainer.begin(); it != jointContainer.end(); ++it) {
            const DataType& currentData = it->getData();
            if (helpers::IsMissingHelper<JointType,DataType>::isMissing(*it,currentData)) continue;
#if defined(_OPENMP)
#pragma omp parallel for
#endif
            for (int j=0; j<(int)contributions.size(); ++j) {
                // inside the loop to make it local for omp
                JointType workjt;
                workjt.setData(currentData);
                workjt.setPrediction(static_cast<DistDataSet<TypeTraits>*>(this)->predictionSampler[j]);
                contributions[j] = 0;
                if (helpers::IsMissingHelper<JointType,PredictionType>::isMissing(static_cast<DistDataSet<TypeTraits>*>(this)->predictionSampler[j])) continue;
                SumContribAction<JointType, SpacialJointSampler> sca(workjt, ssampler, jointKernel, contributions[j]);
                // find neighbours for workjt
                jnf.find(
                    helpers::SingleElementContainerWrapper<JointType>(workjt),
                    jointKernel.neighbourhoodSize(),
                    sca
                );
            }
            Distribution currentDistribution(contributions);
            addDist(currentData,currentDistribution);
        }
    }

    template<class SpacialJointSampler>
    void addBruteForce(JointType jt, JointKernel jointKernel, SpacialJointSampler& ssampler) {
        addBruteForce(helpers::SingleElementContainerWrapper<JointType>(jt), jointKernel, ssampler);
    }

    // Non-separable kernel case, brute-force = no neighbourhood search (faster when most points are within kernel distance of each other)
    template<class JointContainer, class SpacialJointSampler>
    void addBruteForce(JointContainer jointContainer, JointKernel jointKernel, SpacialJointSampler& ssampler) {
        std::vector<FloatType> contributions(static_cast<DistDataSet<TypeTraits>*>(this)->predictionSampler.size(), FloatType(0));
        for (typename JointContainer::iterator it = jointContainer.begin(); it != jointContainer.end(); ++it) {
            const DataType& currentData = it->getData();
            if (helpers::IsMissingHelper<JointType,DataType>::isMissing(*it,currentData)) continue;
#if defined(_OPENMP)
#pragma omp parallel for
#endif
            for (int j=0; j<(int)contributions.size(); ++j) {
                JointType workjt;
                workjt.setData(currentData);
                workjt.setPrediction(static_cast<DistDataSet<TypeTraits>*>(this)->predictionSampler[j]);
                FloatType contrib = 0;
                if (!helpers::IsMissingHelper<JointType,PredictionType>::isMissing(static_cast<DistDataSet<TypeTraits>*>(this)->predictionSampler[j]))
                for (int i=0; i<(int)ssampler.size(); ++i) {
                    contrib += jointKernel(workjt, ssampler[i]);
                }
                contributions[j] = contrib;
            }
            Distribution currentD(contributions);
            addDist(currentData,currentD);
        }
    }

};

}

#endif
