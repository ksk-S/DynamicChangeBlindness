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

#ifndef DECISIONAL_STATES_SYMBOL_CONSTRAINTS_H
#define DECISIONAL_STATES_SYMBOL_CONSTRAINTS_H

#include "decisional_states_helpers.hpp"

#include <set>
#include <map> // ordered
#include <boost/unordered_map.hpp>

namespace decisional_states {

    
// TransitionFeeder does not support symbols => no constraint
template <class TransitionFeeder, class DistributionManager, class TypeTraits, int hasSymbol = 0>
struct SymbolConstraintsBase {
    typedef std::pair<int,int>* iterator;
    typedef std::pair<int,int>* const_iterator;
    iterator begin() const {return 0;}
    iterator end() const {return 0;}
    typedef typename TypeTraits::FloatType FloatType;
    SymbolConstraintsBase(TransitionFeeder&, DistributionManager&, FloatType) {}
    template<class IndexClassProvider>
    void getSplitConstraints(IndexClassProvider indexClassProvider) {}
    template<class IndexClassProvider>
    void getMergeConstraints(IndexClassProvider indexClassProvider) {}
};

// TransitionFeeder supports symbols, take them into account
template <class TransitionFeeder, class DistributionManager, class TypeTraits>
struct SymbolConstraintsBase<TransitionFeeder, DistributionManager,TypeTraits,1> {
    std::set< std::pair<int,int> > constraints;
    typedef typename std::set< std::pair<int,int> >::iterator iterator;
    typedef typename std::set< std::pair<int,int> >::const_iterator const_iterator;
    const_iterator begin() const {return constraints.begin();}
    const_iterator end() const {return constraints.end();}
    iterator begin() {return constraints.begin();}
    iterator end() {return constraints.end();}
    TransitionFeeder& transitionFeeder;
    DistributionManager& distributionManager;
    typedef typename TypeTraits::FloatType FloatType;
    FloatType tolerance;

    SymbolConstraintsBase(TransitionFeeder& _transitionFeeder, DistributionManager& _distributionManager, FloatType _tolerance)
    : transitionFeeder(_transitionFeeder), distributionManager(_distributionManager), tolerance(_tolerance) {}

    typedef typename TypeTraits::DataType DataType;
    typedef typename helpers::SymbolTypeWrapper<TransitionFeeder>::SymbolType SymbolType;
    typedef std::map<int, int> ClassCount;
    typedef boost::unordered_map<SymbolType, ClassCount> SymbolMap;
    typedef boost::unordered_map<int,SymbolMap> IndexedSymbolMap;

    struct TmpInfo {
        int index;
        void* info; // original info
        TmpInfo(int _index, void* _info) : index(_index), info(_info) {}
    };

    /// must provide pairs of (index, class to be merged)
    template<class IndexClassProvider>
    void getMergeConstraints(IndexClassProvider indexClassProvider) {
        constraints.clear();
        IndexedSymbolMap indexedSymbolMap;

        if (tolerance <= 0) return;

        if (transitionFeeder.begin()==transitionFeeder.end()) return;

        // temporarily squat the info field of the distributions
        // we have to provide constraints expressed in terms of the distributionManager indices
        // Logically, constraints on distributions that must be put in the same state
        int i = 0;
        for (typename DistributionManager::iterator distit = distributionManager.begin(); distit != distributionManager.end(); ++distit, ++i) {
            void* info = distributionManager.getInfo(*distit);
            distributionManager.setInfo(*distit, new TmpInfo(i,info));
        }

        for (typename TransitionFeeder::iterator tit = transitionFeeder.begin(); tit != transitionFeeder.end(); ++tit) {
            DataType prevx = transitionFeeder.getDataBeforeTransition(tit);
            DataType x = transitionFeeder.getDataAfterTransition(tit);
            SymbolType s = transitionFeeder.getSymbol(tit);
            int previndex = ((TmpInfo*)distributionManager.getInfo(distributionManager.distribution(prevx)))->index;
            int index = ((TmpInfo*)distributionManager.getInfo(distributionManager.distribution(x)))->index;
            int prevclass = indexClassProvider(previndex);
            int currentclass = indexClassProvider(index);
            SymbolMap& smap = indexedSymbolMap[prevclass];
            // count all transitions, ignore spurious ones below tolerance
            ++smap[s][currentclass];
        }
        for (typename IndexedSymbolMap::iterator it = indexedSymbolMap.begin(); it != indexedSymbolMap.end(); ++it) {
            SymbolMap& smap = it->second;
            // for each symbol
            for (typename SymbolMap::iterator sit = smap.begin(); sit != smap.end(); ++sit) {
                // counts
                int total = 0;
                for (typename ClassCount::iterator cit = sit->second.begin(); cit != sit->second.end(); ++cit) {
                    total += cit->second;
                }
                int mergeID = -1; bool hasMergeID = false;
                // counts above threshold should be merged. Only one class allowed above threshold
                for (typename ClassCount::iterator cit = sit->second.begin(); cit != sit->second.end(); ++cit) {
                    // merge second and after with first
                    if (cit->second >= tolerance * total) {
                        if (hasMergeID) {
                            constraints.insert(std::make_pair(std::min(mergeID,cit->first),std::max(mergeID,cit->first)));
                        } else {
                            hasMergeID = true;
                            mergeID = cit->first;
                        }
                    }
                }
            }
        }

        // restore original info fields
        for (typename DistributionManager::iterator distit = distributionManager.begin(); distit != distributionManager.end(); ++distit) {
            TmpInfo* tmpinfo = (TmpInfo*)distributionManager.getInfo(*distit);
            distributionManager.setInfo(*distit, tmpinfo->info);
            delete tmpinfo;
        }
    }


    /// must provide pairs of (index, class to be splitted from)
    template<class IndexClassProvider>
    void getSplitConstraints(IndexClassProvider indexClassProvider) {
        constraints.clear();
        IndexedSymbolMap indexedSymbolMap;

        if (tolerance <= 0) return;

        if (transitionFeeder.begin()==transitionFeeder.end()) return;

        int i = 0;
        for (typename DistributionManager::iterator distit = distributionManager.begin(); distit != distributionManager.end(); ++distit, ++i) {
            void* info = distributionManager.getInfo(*distit);
            // split phase: act on indices and split from class
            distributionManager.setInfo(*distit, new TmpInfo(i,info));
        }

        for (typename TransitionFeeder::iterator tit = transitionFeeder.begin(); tit != transitionFeeder.end(); ++tit) {
            DataType prevx = transitionFeeder.getDataBeforeTransition(tit);
            DataType x = transitionFeeder.getDataAfterTransition(tit);
            SymbolType s = transitionFeeder.getSymbol(tit);
            int previndex = ((TmpInfo*)distributionManager.getInfo(distributionManager.distribution(prevx)))->index;
            int index = ((TmpInfo*)distributionManager.getInfo(distributionManager.distribution(x)))->index;
            int prevclass = indexClassProvider(previndex);
            int currentclass = indexClassProvider(index);
            SymbolMap& smap = indexedSymbolMap[prevclass];
            // count all transitions, ignore spurious ones below tolerance
            ++smap[s][currentclass];
        }
        
        for (typename IndexedSymbolMap::iterator it = indexedSymbolMap.begin(); it != indexedSymbolMap.end(); ++it) {
            SymbolMap& smap = it->second;
            // for each symbol
            for (typename SymbolMap::iterator sit = smap.begin(); sit != smap.end(); ++sit) {
                // counts
                int total = 0;
                for (typename ClassCount::iterator cit = sit->second.begin(); cit != sit->second.end(); ++cit) {
                    total += cit->second;
                }
                // keep only counts above threshold
                bool hasIDtoKeep = false;
                for (typename ClassCount::iterator cit = sit->second.begin(); cit != sit->second.end();) {
                    if (cit->second < tolerance * total) {
                        typename ClassCount::iterator ncit = cit; ++ncit;
                        sit->second.erase(cit);
                        cit = ncit;
                    }
                    else {
                        if (hasIDtoKeep) {
                            cit->second = -1; // see below
                        } else {
                            hasIDtoKeep = true;
                            cit->second = 1;
                        }
                        ++cit;
                    }
                }
            }
        }
        
        // run over the data again, this time generate split constraints
        for (typename TransitionFeeder::iterator tit = transitionFeeder.begin(); tit != transitionFeeder.end(); ++tit) {
            DataType prevx = transitionFeeder.getDataBeforeTransition(tit);
            DataType x = transitionFeeder.getDataAfterTransition(tit);
            SymbolType s = transitionFeeder.getSymbol(tit);
            int previndex = ((TmpInfo*)distributionManager.getInfo(distributionManager.distribution(prevx)))->index;
            int index = ((TmpInfo*)distributionManager.getInfo(distributionManager.distribution(x)))->index;
            int prevclass = indexClassProvider(previndex);
            int currentclass = indexClassProvider(index);
            SymbolMap& smap = indexedSymbolMap[prevclass];
            typename ClassCount::iterator cit = smap[s].find(currentclass);
            if (cit != smap[s].end() && cit->second==-1) {
                constraints.insert(std::make_pair(previndex,prevclass));
            }
        }

        // restore original info fields
        for (typename DistributionManager::iterator distit = distributionManager.begin(); distit != distributionManager.end(); ++distit) {
            TmpInfo* tmpinfo = (TmpInfo*)distributionManager.getInfo(*distit);
            distributionManager.setInfo(*distit, tmpinfo->info);
            delete tmpinfo;
        }
    }

};

template <class TransitionFeeder, class DistributionManager, class TypeTraits>
struct SymbolConstraints : public SymbolConstraintsBase<TransitionFeeder,DistributionManager,TypeTraits, helpers::SymbolTypeWrapper<TransitionFeeder>::found> {
    typedef SymbolConstraintsBase<TransitionFeeder,DistributionManager,TypeTraits, helpers::SymbolTypeWrapper<TransitionFeeder>::found> Parent;
    typedef typename Parent::iterator iterator;
    typedef typename Parent::const_iterator const_iterator;
    typedef typename TypeTraits::FloatType FloatType;
    SymbolConstraints(TransitionFeeder& feeder, DistributionManager& distributionManager, FloatType tolerance) : Parent(feeder, distributionManager, tolerance) {}
};

}

#endif
