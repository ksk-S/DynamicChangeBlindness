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

#ifndef DECISIONAL_STATES_KERNEL_DISTRIBUTION_MANAGER_H
#define DECISIONAL_STATES_KERNEL_DISTRIBUTION_MANAGER_H

#include "decisional_states_helpers.hpp"
#include "decisional_states_kernels.hpp"
#include "decisional_states_samplers.hpp"
#include "decisional_states_indexed_contribution_distribution.hpp"

#include <boost/range.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/key_extractors.hpp>
#include <boost/unordered_map.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/random.hpp>

#include <vector>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif


#include <sys/time.h>

namespace decisional_states {

template<class _DataSet, class _DataKernel, class _PredictionKernel, class _PredictionSampleContainer = RegularGridSampler<typename helpers::DefaultTypeTraitsDS<_DataSet>::PredictionType>, class TypeTraits = helpers::DefaultTypeTraitsDS<_DataSet> >
struct KernelDistributionManager {
    typedef _DataSet DataSet;
    typedef _DataKernel DataKernel;
    typedef _PredictionKernel PredictionKernel;
    // ForwardTraversalReadableRange, with optionally a size() function defined
    typedef _PredictionSampleContainer PredictionSampleContainer;
    typedef typename TypeTraits::FloatType FloatType;
    typedef typename TypeTraits::PredictionType PredictionType;
    typedef typename TypeTraits::DataType DataType;
    typedef typename TypeTraits::DataNeighbourhoodFinder DataNeighbourhoodFinder;
    typedef typename TypeTraits::PredictionNeighbourhoodFinder PredictionNeighbourhoodFinder;

    // Made the distribution a standalone class so its type is not clobbered
    typedef IndexedContributionDistribution<FloatType> Distribution;

    // Important: the first field in this structure is the array
    // This allows to find back distributions even from an incomplete KDM type
    helpers::SizedArray<Distribution> distributions;

    FloatType probaThreshold;

    //typedef boost::multi_array<FloatType, 2> KZArray;

    // ref to the data set only used in the constructor
    // further access are redirections to unique x/z maps
    DataSet& dataset;
    DataKernel dataKernel;
    PredictionKernel predictionKernel;
    PredictionSampleContainer samples;
    int ndata;
    int nsamples;

    /**
        NOTE ON DATA STRUCTURES:
        - In a first version performance was favored at all costs over memory.
        - This was not viable, memory usage blowed up
        - The current version tries to make a better compromise.
    */

    struct IndexedSizedArrayExtractor {
        typedef helpers::SizedArray<int> result_type;
        result_type& operator()(boost::shared_ptr<helpers::SizedArray<int> > ptr) const {return *ptr;}
    };

    // we need to store known index arrays for sharing them between distributions
    // storing the shared_ptrs in a container impose all objects to have the same
    // life time as the container
    // weak_ptr do not work: the object is destroyed but not removed from the container
    // using intrusive_ptr and provide a deallocation function could to the trick
    // but here we simply only store arrays that have the same life time
    // as the distribution manager, and no other (see aggregate in Distribution for example)
    // Use multi-index advanced features:
    // - store pointers but map on content
    // - possibility to lookup (find) based on compatible key, need not create object to find it!
    typedef boost::multi_index_container<
        boost::shared_ptr<helpers::SizedArray<int> >,
        boost::multi_index::indexed_by<
            boost::multi_index::hashed_unique< IndexedSizedArrayExtractor >
        >
    > AllSampleIdxContainer;
    AllSampleIdxContainer allSampleIdx;


    // Moving expectation out of distribution makes the distribution a standalone object
    // that does not need a ref back to its manager (or to samples)
    template<class Functor>
    FloatType expectation(const Distribution& dist, Functor f) {
        FloatType ret = 0;
        int numentries = dist->idxarray->size();
#ifdef _OPENMP
#pragma omp parallel for reduction(+:ret)
#endif
        for(int i = 0; i<numentries; ++i) {
            ret += f(samples[(*dist->idxarray)[i]]) * dist->contributions[i];
        }
        return ret / dist->sumcontribs;
    }

    struct IndexedContribution {
        int index;
        FloatType contribution;
    };

    // utility for main program file
    struct GreaterContribution {
        inline bool operator()(const IndexedContribution& a, const IndexedContribution& b) const {return a.contribution > b.contribution;}
    };

    // Unique X => influencial neighbors & distribution
    struct XEntry {
        // key: unique X that is pointed to
        DataType x;
        mutable int count;
        // Distribution for this x: usage on training set
        //mutable Distribution* distribution;  // non-key
        // Other information needed to generalise/build distribution for unknown X values
        // pairs (i, K(x_i, x) ) of influencial neighbors for this x. i = index in unique x
        // NOTE: unordered indices here
        mutable helpers::SizedArray<IndexedContribution> influencialNeighbours;
        // NOTE: ordered indices here
        mutable helpers::SizedArray< std::pair<int,int> > zidxcount;  // zidx in unique Z / number of times it appeared
        // default-construct distribution - x map built in a first pass
        XEntry(const DataType& _x) : x(_x), count(1) {}
    };

    // multi-index container trick to benefit from random access
    // memory: better to add an index here than storing iterators in xentry (size int < size iterator, count also that entries are unique)
    // perf: better to store iterators
    typedef boost::multi_index_container<
        XEntry,
        boost::multi_index::indexed_by<
            boost::multi_index::hashed_unique< boost::multi_index::member<XEntry,DataType,&XEntry::x> >,
            boost::multi_index::random_access<>
        >
    > UniqueXContainer;
    UniqueXContainer uniqueX;

    // unique Z in data set. See comments in X structure
    struct ZEntry {
        PredictionType z;
        // see Distribution
        mutable boost::shared_ptr<helpers::SizedArray<int> > idxarray;
        mutable boost::shared_array<FloatType> contributions;
        //mutable helpers::SizedArray<IndexedContribution> influencedSamples;
        ZEntry(const PredictionType& _z) : z(_z), contributions(0) {}
    };
    typedef boost::multi_index_container<
        ZEntry,
        boost::multi_index::indexed_by<
            boost::multi_index::hashed_unique< boost::multi_index::member<ZEntry,PredictionType,&ZEntry::z> >,
            boost::multi_index::random_access<>
        >
    > UniqueZContainer;
    UniqueZContainer uniqueZ;

    // Functor for taking into account points found in the nearest neighbours search
    // might be called in parallel
    struct XNeighboursFoundAction {
        KernelDistributionManager* dm;
        inline XNeighboursFoundAction(KernelDistributionManager* _dm) : dm(_dm) {}
        // allows for passing plain pointers in an array
        // requires random access iterator for efficiency
        template<typename Iterator>
        void operator()(std::size_t index, const Iterator& begin, const Iterator& end) {
            // might be called in //, but each index is unique here so we can safely modify its entry
            const DataType& x = dm->uniqueX.template get<1>()[index].x;
            helpers::SizedArray<IndexedContribution>& contribs = dm->uniqueX.template get<1>()[index].influencialNeighbours;
            int numentries = end - begin;
            //helpers::SizedArray<IndexedContribution>(numentries).swap(contribs);
            contribs.reset(numentries);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(int i = 0; i<numentries; ++i) {
                Iterator b = begin + i;
                contribs[i].index = *b;
                contribs[i].contribution = dm->dataKernel(dm->uniqueX.template get<1>()[*b].x, x);
            }

            // keep only contributions up to probaThreshold proba
            FloatType sumcontrib = 0; int n = 0;
            for(; n<numentries; ++n) {
                sumcontrib += contribs[n].contribution;
                if (contribs[n].contribution < dm->probaThreshold * sumcontrib) break;
            }
            if (n<numentries) {
                helpers::SizedArray<IndexedContribution> tmp(n);
                for (int i=0; i<n; ++i) tmp[i] = contribs[i];
                contribs.swap(tmp);
            }
        }
    };

    struct ZNeighboursFoundAction {
        KernelDistributionManager* dm;
        inline ZNeighboursFoundAction(KernelDistributionManager* _dm) : dm(_dm) {}
        template<typename Iterator>
        void operator()(std::size_t index, const Iterator& begin, const Iterator& end) {
            const ZEntry& ze = dm->uniqueZ.template get<1>()[index];
            // We first need to sort out the indices for comparisons later on
            // copy, then sort
            boost::shared_ptr<helpers::SizedArray<int> > array = boost::shared_ptr<helpers::SizedArray<int> >(
                // copy-construct the range into an array
                new helpers::SizedArray<int>(boost::iterator_range<Iterator>(begin,end))
            );
            std::sort(array->begin(), array->end());
// Actions might be called concurently
// We need to serialise access to the multi_index object as it is modified, unlike X case above
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                typename AllSampleIdxContainer::iterator sit = dm->allSampleIdx.find(*array);
                if (sit == dm->allSampleIdx.end()) sit = dm->allSampleIdx.insert(array).first;
                // indices are shared
                ze.idxarray = *sit;
            }
            int numentries = ze.idxarray->size();
            // contributions are unique to this unique z
            ze.contributions = boost::shared_array<FloatType>(new FloatType[numentries]);
            const PredictionType& z = ze.z;
// When actions are not parallelised, we might still do so here. Otherwise, negligible impact on perf
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(int i = 0; i<numentries; ++i) {
                ze.contributions[i] = dm->predictionKernel(dm->samples[(*ze.idxarray)[i]], z);
            }
/*if (index < 20) {
for(int i = 0; i<numentries; ++i) std::cout << ze.contributions[i] << " ";
std::cout << std::endl;
}*/
        }
    };

    int distCount;
    FloatType distSum;
    struct XNearestNeighbourFoundAction {
        KernelDistributionManager* dm;
        DataNeighbourhoodFinder* dnf;
        inline XNearestNeighbourFoundAction(KernelDistributionManager* _dm, DataNeighbourhoodFinder* _dnf) : dm(_dm), dnf(_dnf) {
            dm->distCount = 0;
            dm->distSum = 0;
        }
        void operator()(std::size_t index, std::size_t neighbourIndex) {
            FloatType dist = dnf->distance(dm->uniqueX.template get<1>()[index].x, dm->uniqueX.template get<1>()[neighbourIndex].x);
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                ++dm->distCount;
                dm->distSum += dist;
            }
        }
    };

    struct ZNearestNeighbourFoundAction {
        KernelDistributionManager* dm;
        PredictionNeighbourhoodFinder* pnf;
        inline ZNearestNeighbourFoundAction(KernelDistributionManager* _dm, PredictionNeighbourhoodFinder* _pnf) : dm(_dm), pnf(_pnf) {
            dm->distCount = 0;
            dm->distSum = 0;
        }
        void operator()(std::size_t index, std::size_t neighbourIndex) {
            FloatType dist = pnf->distance(dm->uniqueZ.template get<1>()[index].z, dm->uniqueZ.template get<1>()[neighbourIndex].z);
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                ++dm->distCount;
                dm->distSum += dist;
            }
        }
    };

    /*
        In a first step: fill above structure
        2nd step: for each sample, follow neighbors in Z and for which X these appeared. Put all these X in the same list.
                  merge lists when an X is already in another list
        3rd step: if there are several lists, then these p(Z|X) never match on any sample => disjoint/parallel causal state search
    */
    struct XEntryToX {
        DataType& operator()(XEntry& xe) {return xe.x;}
        const DataType& operator()(const XEntry& xe) const {return xe.x;}
        typedef const DataType& result_type;
    };
    struct ZEntryToZ {
        PredictionType& operator()(ZEntry& ze) {return ze.z;}
        const PredictionType& operator()(const ZEntry& ze) const {return ze.z;}
        typedef const PredictionType& result_type;
    };

    KernelDistributionManager(DataSet& _dataset, DataKernel _dataKernel, PredictionKernel _predictionKernel, PredictionSampleContainer _samples, FloatType pthreshold = 0)
    : dataset(_dataset), dataKernel(_dataKernel), predictionKernel(_predictionKernel), samples(_samples),
      ndata(helpers::ContainerSizeWrapper<DataSet>::size(_dataset)),
      nsamples(helpers::ContainerSizeWrapper<PredictionSampleContainer>::size(_samples))
    {
        probaThreshold = pthreshold;
        // step 1: collect unique values of x and z in the data set
        for (typename KernelDistributionManager::DataSet::iterator it = dataset.begin(); it!= dataset.end(); ++it) {
//std::cout << it->ori << " " << it->x << " " << it->y << std::endl;
            DataType x = helpers::DataSetIteratorUnrefHelper<DataSet,TypeTraits>::getData(dataset,it);
            if (helpers::IsMissingHelper<DataSet,DataType>::isMissing(dataset,x)) continue;
            PredictionType z = helpers::DataSetIteratorUnrefHelper<DataSet,TypeTraits>::getPrediction(dataset,it);
            if (helpers::IsMissingHelper<DataSet,PredictionType>::isMissing(dataset,z)) continue;
            typename UniqueXContainer::iterator xit = uniqueX.find(x);
            if (xit==uniqueX.end()) xit = uniqueX.insert(XEntry(x)).first;
            else ++xit->count;
            typename UniqueZContainer::iterator zit = uniqueZ.find(z);
            if (zit==uniqueZ.end()) zit = uniqueZ.insert(ZEntry(z)).first;
            // get a random access iterator from the hashed unique iterator of z
            typename UniqueZContainer::template nth_index<1>::type::iterator razit = uniqueZ.template project<1>(zit);
            // find back index and store it in x
            // random access container: further inserts placed with higher indices, this index will remain valid
            int index = razit - uniqueZ.template get<1>().begin();
            // a find/insert would be more efficient, no std::pair creation for nothing
            int insertplace = xit->zidxcount.insert_unique_sorted(std::make_pair(index,0), helpers::LessOnFirst<std::pair<int,int> >());
            // whether the element existed or not we can increase the count (make_pair start at 0)
            ++xit->zidxcount[insertplace].second;
        }

//std::cout << "Nx: " << uniqueX.size() << std::endl;
//std::cout << "Nz: " << uniqueZ.size() << std::endl;

/*
timeval tv; gettimeofday(&tv,0);
unsigned long long thetime = tv.tv_sec * 1000000 + tv.tv_usec;
*/
        // Step 2: find neighbour points in X now the structure and indices are stable
        typedef helpers::ContainerAdaptor<typename UniqueXContainer::template nth_index<1>::type, DataType, XEntryToX> AdaptedX;
        DataNeighbourhoodFinder dnf;
        // pass container of reference points
        dnf.setup(AdaptedX(uniqueX.template get<1>()));

        // Auto-guess kernel size ?
        if (!dataKernel.valid()) {
            dnf.findNearest(AdaptedX(uniqueX.template get<1>()), XNearestNeighbourFoundAction(this,&dnf));
            // rule of thumb: take stdev = half dist avg
            FloatType ksize = distSum / distCount;
            ksize *= ksize;
//std::cout << "data kernel auto-size: " << ksize << std::endl;
            dataKernel.setSize(ksize);
        }

        // to find the neighbors for each point
        // this interface allows for parallel computations, including GPU
        dnf.find(
            // points at which the search should be performed
            AdaptedX(uniqueX.template get<1>()),
            // max distance for the search. Note: It's up to the user to provide neighbourhood finder and kernel with consistent distance/norm
            // default: both normal kernel and exhaustive finder use squared euclidian arguments
            dataKernel.neighbourhoodSize(),
            // action to perform for each search point, with the search point index and a container of indices of reference points as argument
            // the actions shall be thread-safe and may be called concurrently
            XNeighboursFoundAction(this)
        );

        // Step 3: find neighbour samples for each Z
        // Auto-guess kernel size ?
        typedef helpers::ContainerAdaptor<typename UniqueZContainer::template nth_index<1>::type, PredictionType, ZEntryToZ> AdaptedZ;
        if (!predictionKernel.valid()) {
            PredictionNeighbourhoodFinder pnf;
            pnf.setup(AdaptedZ(uniqueZ.template get<1>()));
            pnf.findNearest(AdaptedZ(uniqueZ.template get<1>()), ZNearestNeighbourFoundAction(this,&pnf));
            // rule of thumb: take stdev = dist avg
            FloatType ksize = distSum / distCount;
            ksize *= ksize;
//std::cout << "prediction kernel auto-size: " << ksize << std::endl;
            // rule of thumb: take variance = half dist avg
            predictionKernel.setSize(ksize);
        }
        // need a new finder, based on sample points, if auto-kern size was required
        PredictionNeighbourhoodFinder pnf;
        // pass container of reference points
        pnf.setup(samples);
        // find neighbours in samples for each z
        pnf.find(
            AdaptedZ(uniqueZ.template get<1>()),
            predictionKernel.neighbourhoodSize(),
            ZNeighboursFoundAction(this)
        );

/*gettimeofday(&tv,0);
thetime = tv.tv_sec* 1000000 + tv.tv_usec - thetime;
std::cout << "time spent in neighbour search: " << (float)thetime * 1e-6f << std::endl;
*/
        // Step 4: Build distributions at each unique X
        int xiend = uniqueX.template get<1>().size();
        helpers::SizedArray<Distribution>(xiend).swap(distributions);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int xi=0; xi<xiend; ++xi) {
            const XEntry& xe = uniqueX.template get<1>()[xi];
            // build temporary map for fast inserts, then copy to storage structures
            std::map<int,FloatType> mergedContribs;
            // merge in contributions at each neighbour
            // This includes this very x if the neighbour search is consistent
            // The advantage are:
            // - a uniform interface, compared to including directly xe
            // - the kernel(x,x) value
            for (int ni = 0; ni < (int)xe.influencialNeighbours.size(); ++ni) {
                const XEntry& ne = uniqueX.template get<1>()[xe.influencialNeighbours[ni].index];
                // repeat above code, but additionally scale by neighbour influence
                for(int zi = 0; zi<(int)ne.zidxcount.size(); ++zi) {
                    const ZEntry& ze = uniqueZ.template get<1>()[ne.zidxcount[zi].first];
                    for (int i=0; i<(int)ze.idxarray->size(); ++i) {
                        mergedContribs[(*ze.idxarray)[i]] += ne.count * ze.contributions[i] * ne.zidxcount[zi].second * xe.influencialNeighbours[ni].contribution;
                    }
                }
            }

            // remove values below threshold from map. step 1: get total
            FloatType initsumcontrib = 0;
            for (typename std::map<int,FloatType>::iterator it = mergedContribs.begin(); it != mergedContribs.end(); ++it) {
                initsumcontrib += it->second;
            }
            // step 2: remove insignificant elements
            for (typename std::map<int,FloatType>::iterator it = mergedContribs.begin(); it != mergedContribs.end();) {
                if (it->second < probaThreshold * initsumcontrib) {
                    typename std::map<int,FloatType>::iterator nextit = it; ++nextit;
                    mergedContribs.erase(it);
                    it = nextit;
                } else ++it;
            }

            // Done. copy map to distribution
            // first allocate objects
            int numentries = mergedContribs.size();
            boost::shared_ptr<helpers::SizedArray<int> > array = boost::shared_ptr<helpers::SizedArray<int> >(
                new helpers::SizedArray<int>(numentries)
            );
            distributions[xi]->contributions = boost::shared_array<FloatType>(new FloatType[numentries]);
            // then copy values
            distributions[xi]->sumcontribs = 0;
            int i = 0;
            for (typename std::map<int,FloatType>::iterator it = mergedContribs.begin(); it != mergedContribs.end(); ++it, ++i) {
                // shared array - or not, see below
                (*array)[i] = it->first;
                // contributions are unique to this distrib, store them directly
                distributions[xi]->contributions[i] = it->second;
                // their total is computed also on the fly
                distributions[xi]->sumcontribs += it->second;
//if (xi<20) std::cout << it->second << " ";
            }
//if (xi<20) std::cout << std::endl;
            // Access to allSampleIdx, multi_index no thread safe
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                // Now look for an existing array. Could have been possible without allocating the array, with a few iterator transforms and functors...
                typename AllSampleIdxContainer::iterator sit = allSampleIdx.find(*array);
                if (sit == allSampleIdx.end()) sit = allSampleIdx.insert(array).first;
                // indices are shared. This saves space, temporary array object is destroyed here if it already existed.
                distributions[xi]->idxarray = *sit;
            }
        }

    }

    std::size_t size() {return distributions.size();}
    typedef typename helpers::SizedArray<Distribution>::iterator iterator;
    iterator begin() { return distributions.begin(); }
    iterator end() { return distributions.end(); }
    Distribution& operator[](int index) {
        return distributions[index];
    }
    typedef Distribution value_type;

    typedef boost::unordered_map<DataType, Distribution> ExtraDists;
    ExtraDists extraDists;

    struct SingleDataContainer {
        const DataType& x;
        SingleDataContainer(const DataType& _x) : x(_x) {}
        typedef const DataType* iterator;
        inline iterator begin() const {return &x;}
        inline iterator end() const {return (&x)+1;}
        std::size_t size() const {return 1;}
        const DataType& operator[](int) const {return x;}
    };

    struct MakeExtraDistribution {
        KernelDistributionManager* dm;
        const DataType& x;
        Distribution& dist;
        MakeExtraDistribution(KernelDistributionManager* _dm, const DataType& _x, Distribution& _dist) : dm(_dm), x(_x), dist(_dist) {
        }
        template<typename Iterator>
        void operator()(std::size_t, const Iterator& begin, const Iterator& end) {
            // See comments in constructor of KernelDistributionManager
            std::map<int,FloatType> mergedContribs;
            // merge contributions at each neighbour
            for (Iterator it = begin; it != end; ++it) {
                const XEntry& ne = dm->uniqueX.template get<1>()[*it];
                // kernel contribution from that neighbour scales z values
                FloatType xscale = dm->dataKernel(ne.x,x);
                for(int zi = 0; zi<(int)ne.zidxcount.size(); ++zi) {
                    const ZEntry& ze = dm->uniqueZ.template get<1>()[ne.zidxcount[zi].first];
                    for (int i=0; i<(int)ze.idxarray->size(); ++i) {
                        mergedContribs[(*ze.idxarray)[i]] += ne.count * ze.contributions[i] * ne.zidxcount[zi].second * xscale;
                    }
                }
            }
            // Done. copy map to distribution
            int numentries = mergedContribs.size();
            dist->idxarray = boost::shared_ptr<helpers::SizedArray<int> >(new helpers::SizedArray<int>(numentries));
            dist->contributions = boost::shared_array<FloatType>(new FloatType[numentries]);
            dist->sumcontribs = 0;
            int i = 0;
            for (typename std::map<int,FloatType>::iterator it = mergedContribs.begin(); it != mergedContribs.end(); ++it, ++i) {
                (*dist->idxarray)[i] = it->first;
                dist->contributions[i] = it->second;
                dist->sumcontribs += it->second;
            }
        }
    };

    Distribution missingDist;
    Distribution& distribution(const DataType& x) {
        if (helpers::IsMissingHelper<DataSet,DataType>::isMissing(x)) return missingDist;
        typename UniqueXContainer::iterator it = uniqueX.find(x);
        if (it != uniqueX.end()) {
            // get a random access iterator from the hashed unique iterator of x. O(1) operation
            typename UniqueXContainer::template nth_index<1>::type::iterator raxit = uniqueX.template project<1>(it);
            // find back index and return the corresponding distribution
            return distributions[raxit - uniqueX.template get<1>().begin()];
        }
        // unknown x, must return an extra dist
        typename ExtraDists::iterator eit = extraDists.find(x);
        if (eit == extraDists.end()) {
            eit = extraDists.insert(std::make_pair(x,Distribution())).first;
            // find neighbour points in X for that new query point
            DataNeighbourhoodFinder dnf;
            typedef helpers::ContainerAdaptor<typename UniqueXContainer::template nth_index<1>::type, DataType, XEntryToX> AdaptedX;
            // pass container of reference points
            dnf.setup(AdaptedX(uniqueX.template get<1>()));
            dnf.find(
                SingleDataContainer(x),
                dataKernel.neighbourhoodSize(),
                // callback from neighbour finding routine. Directly creates the distribution
                MakeExtraDistribution(this,x,eit->second)
            );
        }
        return eit->second;
    }

    void* getInfo(Distribution& dist) {
        return dist->info;
    }
    void setInfo(Distribution& dist, void* info) {
        dist->info = info;
    }

    int getCount(Distribution& dist) {
        int idx = &dist - &distributions[0];
        if (idx >= 0 && idx < (int)distributions.size()) {
            // true idx, address match in the vector
            return uniqueX.template get<1>()[idx].count;
        }
        return 1; // extra dist
    }

    typedef typename Distribution::BhattacharyyaMatcher DistributionMatcher;

};



}

#endif
