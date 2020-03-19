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

#ifndef DECISIONAL_STATES_H
#define DECISIONAL_STATES_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>

#include "decisional_states_helpers.hpp"
#include "decisional_states_discrete_distribution_manager.hpp"
//#include "decisional_states_single_link_clustering.hpp"
#include "decisional_states_aggregate_clustering.hpp"
#include "decisional_states_symbol_constraints.hpp"
#include "decisional_states_optimisers.hpp"

#include <vector>
#include <list>
#include <utility> // for std::pair

#include <assert.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
//#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/key_extractors.hpp>
#include <boost/shared_ptr.hpp>
//#include <boost/unordered_set.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

namespace decisional_states {

// Export this in the namespace so the user can change the default argument if needed.
template<typename FloatType>
struct ApproximateUtilityMatcher {
    FloatType epsilon;
    ApproximateUtilityMatcher(FloatType _epsilon = 1e-6f) : epsilon(_epsilon) {}
    bool operator()(FloatType a, FloatType b) {
        return helpers::fabs(a - b) <= epsilon;
    }
};

template<class _DistributionManager, class _DataSet, class _Optimiser = NullOptimiser, class TransitionFeeder = typename helpers::DefaultTypeTraits<_DataSet,_Optimiser>::TransitionFeeder,  class TypesTraits = helpers::DefaultTypeTraits<_DataSet,_Optimiser> >
class DecisionalStatesAnalyser {
public:
    typedef _DataSet DataSet;
    typedef _Optimiser Optimiser;
    typedef _DistributionManager DistributionManager;
    typedef typename helpers::DistributionWrapper<DistributionManager>::Distribution Distribution;
    //FloatType defaults to double, Prediction and Data would need a decltype to default to something
    typedef typename TypesTraits::FloatType FloatType;
    typedef typename TypesTraits::PredictionType PredictionType;
    typedef typename TypesTraits::DataType DataType;
    typedef typename TypesTraits::PredictionSet PredictionSet;
    // symbol type might be void
    typedef typename helpers::SymbolTypeWrapper<TransitionFeeder>::SymbolType SymbolType;
    
    /// user can provide either a custom class for the transitions
    DecisionalStatesAnalyser(DistributionManager& _manager, TransitionFeeder& _transitionFeeder, Optimiser& _optimiser) : distributionManager(_manager), transitionFeeder(&_transitionFeeder), ownedTransitionFeeder(false), optimiser(&_optimiser), ownedOptimiser(false) {}
    
    /// or simply the data set, we'll build a transition feeder around it (possibly a default one)
    DecisionalStatesAnalyser(DistributionManager& _manager, DataSet& _dataset, Optimiser& _optimiser) : distributionManager(_manager), transitionFeeder(new TransitionFeeder(_dataset)), ownedTransitionFeeder(true), optimiser(&_optimiser), ownedOptimiser(false) {}
    
    DecisionalStatesAnalyser(DistributionManager& _manager, TransitionFeeder& _transitionFeeder) : distributionManager(_manager), transitionFeeder(&_transitionFeeder), ownedTransitionFeeder(false), optimiser(new Optimiser()), ownedOptimiser(true) {}
    
    /// or simply the data set, we'll build a transition feeder around it (possibly a default one)
    DecisionalStatesAnalyser(DistributionManager& _manager, DataSet& _dataset) : distributionManager(_manager), transitionFeeder(new TransitionFeeder(_dataset)), ownedTransitionFeeder(true), optimiser(new Optimiser()), ownedOptimiser(true) {}
    
    /// Works only if the objet was not duplicated...
    /// TODO copy ops, shared ptr for feeder, etc
    ~DecisionalStatesAnalyser() {   
        if (ownedTransitionFeeder) delete transitionFeeder;
        if (ownedOptimiser) delete optimiser;
    }

    /// internal
    DistributionManager& distributionManager;
    // need to handle both an external ref and an owned object
    // => ptr + ownership flag
    TransitionFeeder* transitionFeeder;
    bool ownedTransitionFeeder;
    Optimiser* optimiser;
    bool ownedOptimiser;

    /** Data Structures:
        StateInfo: stored for each distribution 1 to 1.
        Causal states:
            - contain N distributions (supposedly all equivalent) and one aggregate dist
            - have a unique prediction set and utility value
        Decisional state: contains causal states, is intersection of IsoUtility and IsoPrediction
        IsoUtility state: contains N decisional states, partition
        IsoPrediction state: contains N decisional states, different partition
    */

/**
    Operational modes:
    - can cluster causal states
      - from data
      - by subdividing decisional/isoutility/isoprediction sets
    - can apply utility functor
      - to data directly
      - to causal states
    - cluster isoutility and/or isoprediction from functor result
    - as soon as both iso-X are available, decisional states automatically computed by intersecting
    - computeDecisionalStates does both auto at the same time.
*/

    struct DecisionalState;
    struct IsoUtilityState;
    struct IsoPredictionState;
    struct CausalState;

    struct StateInfo {                            // CS first         DS first
        // Distribution is derived from causal states (CS first) or from data (DS first)
        // the manager.setInfo points to resp. the common info in CS, or to individual info structs
        Distribution distribution;                // 2                1
        // Utility/predictions are always computed from the above distributions
        FloatType utility;                        // 3                2
        PredictionSet predictionSet;              // 3                2
        // Causal states handle their own distributions - either from data directly (CS) or using the DS as containers
        boost::shared_ptr<CausalState> causalState;                 // 1                5
        // Iso sets are always clustered from the above quantities
        boost::shared_ptr<IsoUtilityState> isoUtilityState;         // 4                3
        boost::shared_ptr<IsoPredictionState> isoPredictionState;   // 4                3
        // decisional states are always computed by intersecting the iso states
        boost::shared_ptr<DecisionalState> decisionalState;         // 5                4
        // Conclusion: causal states computed either first or last.
        // Thus we use causalStates.empty() as a flag if they were computed first or not

        // constructor: default objects and 0 pointers
        StateInfo() {} //: causalState(0), isoUtilityState(0), isoPredictionState(0), decisionalState(0) {}
    };

    // container: no key - can cluster together identical elements if needed
    std::vector<StateInfo> stateInfos;
    // Extra StateInfo objects for points that were met after the clustering
    // kept in a list so they are owned by this object and automatically destroyed
    // the list instead of vector ensures the pointers remain valid
    std::list<StateInfo> extraInfos;

    /* Graph structures
       Could use boost::graph ? Except we need to integrate our other structs
    */
    template<class State, typename Label = void>
    struct TransitionKey {
        Label label;
        // Anonymous union, thanks C++ designers for this funny beast
        union {
            // member state pointer, as in non-symbol version
            mutable State* state;
            // temporary map used while making epsilon-machine
            mutable std::map<State*,FloatType>* states;
        };
        typedef boost::multi_index::member<TransitionKey, Label, &TransitionKey::label> Key;
        Label getKey() {return label;}
        inline void initCountForState(State* _state) const {
            states = new std::map<State*,FloatType>();
            (*states)[_state] = 1;
        }
        inline void incrementCountForState(State* _state) const {
            // prerequisite: union in count mode
            ++(*states)[_state];
        }
        inline void retainMaxCountState() const {
            int maxCount = -1; State* maxCState = 0;
            for (typename std::map<State*,FloatType>::iterator it = states->begin(); it!=states->end(); ++it) {
                if (it->second>maxCount) {
                    maxCount = it->second;
                    maxCState = it->first;
                }
            }
            delete states;
            state = maxCState;
        }
        inline void rescaleStatesProba(int total) const {
            for (typename std::map<CausalState*,FloatType>::iterator sit = states->begin(); sit!=states->end(); ++sit) {
                sit->second /= total;
            }
        }
    };
    template<class State>
    struct TransitionKey<State,void> {
        State* state;
        State* states;
        typedef boost::multi_index::member<TransitionKey, State*, &TransitionKey::state> Key;
        inline State* getKey() {return state;}
        inline void incrementCountForState(State*) const {}
        inline void initCountForState(State*) const {}
        inline void retainMaxCountState() const {}
        inline void rescaleStatesProba(int total) const {}
        TransitionKey<State,void>() : state(0), states(0) {}
        TransitionKey<State,void>(State* key) : state(key), states(0) {}
    };
    template<class State, typename Label = void>
    struct Transition : public TransitionKey<State, Label> {
        typedef typename TransitionKey<State, Label>::Key Key;
        mutable FloatType probability;
        Transition() {}
        Transition(State* key) : TransitionKey<State, Label>(key) {}
    };

    typedef Transition<CausalState, typename helpers::SymbolTypeWrapper<TransitionFeeder>::SymbolType> CausalStateTransition;
    typedef Transition<IsoUtilityState> IsoUtilityStateTransition;
    typedef Transition<IsoPredictionState> IsoPredictionStateTransition;
    typedef Transition<DecisionalState> DecisionalStateTransition;
    // template typedef will be nice in C++0x
    typedef boost::multi_index_container<CausalStateTransition,
        boost::multi_index::indexed_by<boost::multi_index::hashed_unique<typename CausalStateTransition::Key> >
    > CausalStateTransitionSet;
    typedef boost::multi_index_container<IsoUtilityStateTransition,
        boost::multi_index::indexed_by<boost::multi_index::hashed_unique<typename IsoUtilityStateTransition::Key> >
    > IsoUtilityStateTransitionSet;
    typedef boost::multi_index_container<IsoPredictionStateTransition,
        boost::multi_index::indexed_by<boost::multi_index::hashed_unique<typename IsoPredictionStateTransition::Key> >
    > IsoPredictionStateTransitionSet;
    typedef boost::multi_index_container<DecisionalStateTransition,
        boost::multi_index::indexed_by<boost::multi_index::hashed_unique<typename DecisionalStateTransition::Key> >
    > DecisionalStateTransitionSet;

    // root state: no predecessor except self-loops, necessarily transient
    // recurrent state: within the recurrent part of the graph. Mutually exclusive with root
    enum {ROOT_STATE = 1, RECURRENT_STATE=2};
    struct State {
        int count;   // for stats
        void* user;  // space for a user property attached to the state
        int flags;   // one or more of the above flags (user extension possible)
        int index;
        inline State() : count(0), user(0), flags(0), index(-1) {}
        inline bool is_root() {return (flags&ROOT_STATE)!=0;}
        inline bool is_recurrent() {return (flags&RECURRENT_STATE)!=0;}
        boost::shared_ptr<std::vector<StateInfo*> > members;
    };
    
    // Each state contains its aggregate representative
    struct CausalState : public State {
        Distribution distribution;
        boost::shared_ptr<CausalStateTransitionSet> transitions;
    };
    struct IsoUtilityState : public State {
        FloatType utility;
        boost::shared_ptr<IsoUtilityStateTransitionSet> transitions;
    };
    struct IsoPredictionState : public State {
        PredictionSet predictionSet;
        boost::shared_ptr<IsoPredictionStateTransitionSet> transitions;
    };
    struct DecisionalState : public State {
        boost::shared_ptr<IsoUtilityState> isoUtilityState;
        boost::shared_ptr<IsoPredictionState> isoPredictionState;
        boost::shared_ptr<DecisionalStateTransitionSet> transitions;
    };

    // containers
    typedef std::vector<boost::shared_ptr<CausalState> > CausalStates;
    typedef std::vector<boost::shared_ptr<IsoUtilityState> > IsoUtilityStates;
    typedef std::vector<boost::shared_ptr<IsoPredictionState> > IsoPredictionStates;
    typedef std::vector<boost::shared_ptr<DecisionalState> > DecisionalStates;
    CausalStates causalStates;
    IsoUtilityStates isoUtilityStates;
    IsoPredictionStates isoPredictionStates;
    DecisionalStates decisionalStates;

    template <class UtilityFunctor>
    struct ExpectedUtilityFunctor {
        DistributionManager& distributionManager;
        Distribution distribution;
        UtilityFunctor& utility;

        ExpectedUtilityFunctor(DistributionManager& dm, Distribution dist, UtilityFunctor& util)
        : distributionManager(dm), distribution(dist), utility(util) {
        }

        FloatType operator()(const PredictionType& y) {
            return helpers::ExpectationHelper<FloatType, DistributionManager>::expectation(distributionManager,distribution, boost::bind<FloatType>(utility, boost::cref(y), _1));
        }

    };

    // Functors are stored for later use in function object wrappers.
    // They are used directly in the functions that cluster states so as to benefit from inlining when possible.
    // When an unknown point is met in a getXXXState function, the callback is then used for consistency with the clustering.
    typedef boost::function<FloatType (const PredictionType&, const PredictionType&)> UtilityFunctorCallback;
    UtilityFunctorCallback utilityFunctorCallback;
    boost::function<FloatType (const Distribution&, const Distribution&)> distributionMatcherCallback;
    boost::function<FloatType (const FloatType&, const FloatType&)> utilityMatcherCallback;
    boost::function<FloatType (const PredictionSet&, const PredictionSet&)> predictionSetMatcherCallback;


    // Apply the given utility functor to the data and derive the utility and prediction values
    template <class UtilityFunctor>
    void applyUtility(UtilityFunctor utilityFunctor) {
        assert(!utilityFunctorCallback);
        utilityFunctorCallback = utilityFunctor;
        // starting from data?
        if (causalStates.empty()) {
            // Yes: create info structs
            assert(stateInfos.empty());
            // size helper: reserve the right amount of memory
            stateInfos.resize(helpers::ContainerSizeWrapper<DistributionManager>::size(distributionManager));
            int i = 0;
            for (typename DistributionManager::iterator it = distributionManager.begin(); it != distributionManager.end(); ++it) {
                stateInfos[i].distribution = helpers::DistributionManagerIteratorHolder<DistributionManager>(distributionManager,it).deref();
                helpers::StateInfoHelper<DistributionManager>::setInfo(distributionManager,stateInfos[i].distribution,&stateInfos[i]);
                ++i;
            }
        }
        // the info vector corresponds to either data distributions or to the causal state entries
        // in either case use it as the reference for computing utility-related values
        int iend = stateInfos.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i=0; i<iend; i++) {
            stateInfos[i].utility = optimiser->optimise(ExpectedUtilityFunctor<UtilityFunctor>(distributionManager, stateInfos[i].distribution, utilityFunctor), stateInfos[i].predictionSet);
        }
    }

    typedef typename helpers::DistributionMatcherWrapper<DistributionManager,Distribution>::DistributionMatcher ClusterCausalStatesDefaultDistributionMatcher;
    typedef AggregateClustering<> ClusterCausalStatesDefaultClusteringAlgorithm;

// TODO: doc.
// The user must provide
// - MatchFunctor: takes as argument (Data,Data) or anything compatible
//   (like references to elements) and returns:
//   - match(a,b) > 0: the elements match
//   - match(a,b) <= 0: the elements do not match
//   This allows boolean predicates to work, and also provides a mean to encode distances
//   for a clustering algorithm (ex: return threshold - distance(a,b); )
//   It's also possible to ignore threshold: return distance(a,b)==0;  or even d==0?1:-d


    // chained invocation Java-like trick to get default arguments to work
    inline bool computeCausalStates() {
        return computeCausalStates(ClusterCausalStatesDefaultClusteringAlgorithm());
    }

    template<class ClusteringAlgorithm>
    inline bool computeCausalStates(ClusteringAlgorithm clusterer) {
        return computeCausalStates(clusterer,ClusterCausalStatesDefaultDistributionMatcher());
    }

    template<class ClusteringAlgorithm, class DistributionMatcher>
    inline bool computeCausalStates(ClusteringAlgorithm clusterer, DistributionMatcher matcher) {
        return computeCausalStates(clusterer,matcher,true);
    }

    template<class ClusteringAlgorithm, class DistributionMatcher>
    inline bool computeCausalStates(ClusteringAlgorithm clusterer, DistributionMatcher matcher, bool checkStatesConsistency) {
        return computeCausalStates(clusterer,matcher,checkStatesConsistency,FloatType(0.95f));
    }

    // cluster causal states
    // - by subdividing decisional / iso-x states if they are present
    // - from data
    // note: a tolerance <=0 disables checking for symbol constraints (non-deterministic automata)
    template<class ClusteringAlgorithm, class DistributionMatcher>
    bool computeCausalStates(ClusteringAlgorithm clusterer, DistributionMatcher matcher, bool checkStatesConsistency, FloatType tolerance) {
        distributionMatcherCallback = matcher;
        // Start from decisional states ?
        if (!decisionalStates.empty()) {
            return subpartitionCausalStatesFromOtherState(matcher,clusterer,decisionalStates,checkStatesConsistency);
        // either none or only one of the iso-x is present. when both are computed, so is decisional states => above case
        // hence the next two tests order does not matter
        } else if (!isoUtilityStates.empty()) {
            return subpartitionCausalStatesFromOtherState(matcher,clusterer,isoUtilityStates,checkStatesConsistency);
        } else if (!isoPredictionStates.empty()) {
            return subpartitionCausalStatesFromOtherState(matcher,clusterer,isoPredictionStates,checkStatesConsistency);
        }
        // from data
        assert(stateInfos.empty());
        assert(causalStates.empty());
        // TODO: detect operator[] and if it does not exist store iterators in a vector to simulate random access
        // build constraints for each unique x... represented as unique distributions
        clusterer.cluster(distributionManager, matcher, SymbolConstraints<TransitionFeeder,DistributionManager,TypesTraits>(*transitionFeeder, distributionManager,tolerance));
        // create info structs: one per causal state for further construction of other states
        stateInfos.resize(helpers::ContainerSizeWrapper<ClusteringAlgorithm>::size(clusterer));
        // create causal state container
        causalStates.resize(helpers::ContainerSizeWrapper<ClusteringAlgorithm>::size(clusterer));
        // Pointers to vector elements are now valid: fill info objects with states
        int i = 0;
        typedef helpers::ClusteringAlgorithmIteratorHolder<ClusteringAlgorithm> CAIteratorHolder;
        for (typename ClusteringAlgorithm::iterator it = clusterer.begin(); it != clusterer.end(); ++it) {
            typedef typename helpers::ClusterWrapper<ClusteringAlgorithm>::Cluster Cluster;
            typedef helpers::ClusterIteratorHolder<Cluster> CIteratorHolder;
            Cluster& cluster = CAIteratorHolder(clusterer,it).deref();
            causalStates[i] = boost::shared_ptr<CausalState>(new CausalState());
            causalStates[i]->index = i;
            causalStates[i]->distribution = helpers::AggregatorHelper<ClusteringAlgorithm,DistributionManager>::aggregate(clusterer,cluster,distributionManager);
            causalStates[i]->members.reset(new std::vector<StateInfo*>(1));
            (*causalStates[i]->members)[0] = &stateInfos[i];
            stateInfos[i].distribution = causalStates[i]->distribution;
            stateInfos[i].causalState = causalStates[i];
            for (typename Cluster::iterator cit = cluster.begin(); cit != cluster.end(); ++cit) {
                Distribution& dist = distributionManager[CIteratorHolder(cluster,cit).deref()];
                // All members of the cluster share the same info
                helpers::StateInfoHelper<DistributionManager>::setInfo(distributionManager,dist,&stateInfos[i]);
                // check consistency of the clustering? stop checking after the first failure
                // aggregate dist must match all member dists
                if (checkStatesConsistency) checkStatesConsistency = matcher(
                    causalStates[i]->distribution,
                    dist
                ) > 0;
            }
            ++i;
        }
        return checkStatesConsistency;
    }

    // internal helper function to avoid code duplication
    template<class DistributionMatcher, class ClusteringAlgorithm, class OtherState>
    bool subpartitionCausalStatesFromOtherState(DistributionMatcher& matcher, ClusteringAlgorithm& clusterer, std::vector<OtherState>& otherStates, bool checkStatesConsistency) {
        typedef helpers::ClusteringAlgorithmIteratorHolder<ClusteringAlgorithm> CAIteratorHolder;
        typedef typename helpers::ClusterWrapper<ClusteringAlgorithm>::Cluster Cluster;
        typedef helpers::ClusterIteratorHolder<Cluster> CIteratorHolder;
        // iterate over the higher-level states and subcluster
        for (typename std::vector<OtherState>::iterator oit = otherStates.begin(); oit!=otherStates.end(); ++oit) {
            ClusteringAlgorithm clust = clusterer;
            typedef helpers::PointerContainerMemberAdaptor<std::vector<StateInfo*>, Distribution, &StateInfo::distribution> AdaptedInfo;
            AdaptedInfo adaptedInfo(*(*oit)->members);
            clust.cluster(adaptedInfo, matcher);
            // clustering done. Retrieve clusters
            for (typename ClusteringAlgorithm::iterator it = clust.begin(); it != clust.end(); ++it) {
                Cluster& cluster = CAIteratorHolder(clust,it).deref();
                boost::shared_ptr<CausalState> cs = boost::shared_ptr<CausalState>(new CausalState());
                causalStates.push_back(cs);
                causalStates.back()->index = causalStates.size()-1;
                cs->distribution = helpers::AggregatorHelper<ClusteringAlgorithm,AdaptedInfo>::aggregate(clust,cluster,adaptedInfo);
                cs->members.reset(new std::vector<StateInfo*>(helpers::ContainerSizeWrapper<Cluster>::size(cluster)));
                int i = 0;
                for (typename Cluster::iterator cit = cluster.begin(); cit != cluster.end(); ++cit) {
                    StateInfo* info = (*(*oit)->members)[CIteratorHolder(cluster,cit).deref()];
                    (*cs->members)[i++] = info;
                    if (checkStatesConsistency) checkStatesConsistency = matcher(
                        cs->distribution,
                        info->distribution
                    ) > 0;
                }
            }
        }
        // the link from state info to causal states is maintained when the causal states vector is fully built
        for (typename std::vector<boost::shared_ptr<CausalState> >::iterator it = causalStates.begin(); it!=causalStates.end(); ++it) {
            for(typename std::vector<StateInfo*>::iterator sit = (*it)->members->begin(); sit != (*it)->members->end(); ++sit) {
                (*sit)->causalState = *it;
            }
        }
        return checkStatesConsistency;
    }

    typedef AggregateClustering<> ClusterIsoPredictionStatesDefaultClusteringAlgorithm;

    // chained invocation Java-like trick to get default arguments to work
    bool computeIsoPredictionStates() {
        return computeIsoPredictionStates(ClusterIsoPredictionStatesDefaultClusteringAlgorithm());
    }

    template<class ClusteringAlgorithm>
    bool computeIsoPredictionStates(ClusteringAlgorithm clusterer) {
        return computeIsoPredictionStates(clusterer,std::equal_to<PredictionSet>());
    }

    template<class ClusteringAlgorithm, class PredictionSetMatcher>
    bool computeIsoPredictionStates(ClusteringAlgorithm clusterer, PredictionSetMatcher matcher) {
        return computeIsoPredictionStates(clusterer,matcher,true);
    }

    template<class ClusteringAlgorithm, class PredictionSetMatcher>
    bool computeIsoPredictionStates(ClusteringAlgorithm clusterer, PredictionSetMatcher matcher, bool checkStatesConsistency) {
        // Always start from the info vector, whatever the way it was produced (data or causal states)
        // If it is empty, then there was no previous clustering or there was nothing to cluster
        if (stateInfos.empty()) return false;
        assert(!utilityFunctorCallback.empty());
        assert(isoPredictionStates.empty());
        predictionSetMatcherCallback = matcher;
        // Cluster the info states according to the prediction field
        typedef helpers::ContainerMemberAdaptor<std::vector<StateInfo>, PredictionSet, &StateInfo::predictionSet> AdaptedInfo;
        AdaptedInfo adaptedInfo(stateInfos);
        clusterer.cluster(adaptedInfo, matcher);
        isoPredictionStates.resize(helpers::ContainerSizeWrapper<ClusteringAlgorithm>::size(clusterer));
        int i = 0;
        typedef helpers::ClusteringAlgorithmIteratorHolder<ClusteringAlgorithm> CAIteratorHolder;
        for (typename ClusteringAlgorithm::iterator it = clusterer.begin(); it != clusterer.end(); ++it) {
            typedef helpers::ClusteringAlgorithmIteratorHolder<ClusteringAlgorithm> CAIteratorHolder;
            typedef typename helpers::ClusterWrapper<ClusteringAlgorithm>::Cluster Cluster;
            typedef helpers::ClusterIteratorHolder<Cluster> CIteratorHolder;
            Cluster& cluster = CAIteratorHolder(clusterer,it).deref();
            isoPredictionStates[i] = boost::shared_ptr<IsoPredictionState>(new IsoPredictionState());
            isoPredictionStates[i]->predictionSet = helpers::AggregatorHelper<ClusteringAlgorithm,AdaptedInfo>::aggregate(clusterer,cluster,adaptedInfo);
            isoPredictionStates[i]->members.reset(new std::vector<StateInfo*>(helpers::ContainerSizeWrapper<Cluster>::size(cluster)));
            int j = 0;
            for (typename Cluster::iterator cit = cluster.begin(); cit != cluster.end(); ++cit) {
                StateInfo* info = &stateInfos[CIteratorHolder(cluster,cit).deref()];
                (*isoPredictionStates[i]->members)[j++] = info;
                info->isoPredictionState = isoPredictionStates[i];
                if (checkStatesConsistency) checkStatesConsistency = matcher(
                    isoPredictionStates[i]->predictionSet,
                    info->predictionSet
                ) > 0;
            }
            ++i;
        }
        for (int i=0; i<(int)isoPredictionStates.size(); ++i) isoPredictionStates[i]->index=i;
        return checkStatesConsistency;
    }

    typedef AggregateClustering<> ClusterIsoUtilityStatesDefaultClusteringAlgorithm;

    // chained invocation Java-like trick to get default arguments to work
    bool computeIsoUtilityStates() {
        return computeIsoUtilityStates(ClusterIsoUtilityStatesDefaultClusteringAlgorithm());
    }

    template<class ClusteringAlgorithm>
    bool computeIsoUtilityStates(ClusteringAlgorithm clusterer) {
        return computeIsoUtilityStates(clusterer,ApproximateUtilityMatcher<FloatType>());
    }

    template<class ClusteringAlgorithm, class UtilityMatcher>
    bool computeIsoUtilityStates(ClusteringAlgorithm clusterer, UtilityMatcher matcher) {
        return computeIsoUtilityStates(clusterer,matcher,true);
    }

    template<class ClusteringAlgorithm, class UtilityMatcher>
    bool computeIsoUtilityStates(ClusteringAlgorithm clusterer, UtilityMatcher matcher, bool checkStatesConsistency) {
        // Always start from the info vector, whatever the way it was produced (data or causal states)
        // If it is empty, then there was no previous clustering or there was nothing to cluster
        if (stateInfos.empty()) return false;
        assert(!utilityFunctorCallback.empty());
        assert(isoUtilityStates.empty());
        utilityMatcherCallback = matcher;
        // Cluster the info states according to the prediction field
        typedef helpers::ContainerMemberAdaptor<std::vector<StateInfo>, FloatType, &StateInfo::utility> AdaptedInfo;
        AdaptedInfo adaptedInfo(stateInfos);
        clusterer.cluster(adaptedInfo, matcher);
        isoUtilityStates.resize(helpers::ContainerSizeWrapper<ClusteringAlgorithm>::size(clusterer));
        int i = 0;
        typedef helpers::ClusteringAlgorithmIteratorHolder<ClusteringAlgorithm> CAIteratorHolder;
        for (typename ClusteringAlgorithm::iterator it = clusterer.begin(); it != clusterer.end(); ++it) {
            typedef helpers::ClusteringAlgorithmIteratorHolder<ClusteringAlgorithm> CAIteratorHolder;
            typedef typename helpers::ClusterWrapper<ClusteringAlgorithm>::Cluster Cluster;
            typedef helpers::ClusterIteratorHolder<Cluster> CIteratorHolder;
            Cluster& cluster = CAIteratorHolder(clusterer,it).deref();
            isoUtilityStates[i] = boost::shared_ptr<IsoUtilityState>(new IsoUtilityState());
            isoUtilityStates[i]->utility = helpers::AggregatorHelper<ClusteringAlgorithm,AdaptedInfo>::aggregate(clusterer,cluster,adaptedInfo);
            isoUtilityStates[i]->members.reset(new std::vector<StateInfo*>(helpers::ContainerSizeWrapper<Cluster>::size(cluster)));
            int j = 0;
            for (typename Cluster::iterator cit = cluster.begin(); cit != cluster.end(); ++cit) {
                StateInfo* info = &stateInfos[CIteratorHolder(cluster,cit).deref()];
                (*isoUtilityStates[i]->members)[j++] = info;
                info->isoUtilityState = isoUtilityStates[i];
                if (checkStatesConsistency) checkStatesConsistency = matcher(
                    isoUtilityStates[i]->utility,
                    info->utility
                ) > 0;
            }
            ++i;
        }
        for (int i=0; i<(int)isoUtilityStates.size(); ++i) isoUtilityStates[i]->index=i;
        return checkStatesConsistency;
    }

    /// Called when both Iso-X states are available. Intersect them
    void computeDecisionalStates() {
        // If it is empty, then there was no previous clustering or there was nothing to cluster
        if (stateInfos.empty()) return;
        assert(!isoPredictionStates.empty());
        assert(!isoUtilityStates.empty());
        assert(decisionalStates.empty());

        typedef std::pair<boost::shared_ptr<IsoUtilityState>,boost::shared_ptr<IsoPredictionState> > DecisionalStateKey;
        typedef boost::shared_ptr<std::vector<StateInfo*> > Members;
        typedef boost::unordered_map<DecisionalStateKey, Members > DecisionalStateMap;
        DecisionalStateMap decisionalStateMap;

        // build the map: find unique combinations of iso-X states
        for(typename std::vector<StateInfo>::iterator it = stateInfos.begin(); it != stateInfos.end(); ++it) {
            Members& members = decisionalStateMap[std::make_pair(it->isoUtilityState,it->isoPredictionState)];
            // new state?
            if (members.get()==0) members = Members(new std::vector<StateInfo*>());
            members->push_back(&(*it));
        }

        decisionalStates.resize(decisionalStateMap.size());
        int i = 0;
        for(typename DecisionalStateMap::iterator it = decisionalStateMap.begin(); it != decisionalStateMap.end(); ++it) {
            decisionalStates[i] = boost::shared_ptr<DecisionalState>(new DecisionalState());
            decisionalStates[i]->isoUtilityState = it->first.first;
            decisionalStates[i]->isoPredictionState = it->first.second;
            decisionalStates[i]->members = it->second; // long life to shared pointers!
            decisionalStates[i]->index = i;
            // link back member info to the state
            for (typename std::vector<StateInfo*>::iterator mit = it->second->begin(); mit != it->second->end(); ++mit) {
                (*mit)->decisionalState = decisionalStates[i];
            }
            ++i;
        }
    }

    boost::shared_ptr<CausalState> getCausalState(DataType x, FloatType matchThreshold = 0, bool store = true) {
        if (!distributionMatcherCallback) return boost::shared_ptr<CausalState>();
        Distribution dist = distributionManager.distribution(x);
        return getCausalState(dist, matchThreshold, store, distributionMatcherCallback);
    }

    boost::shared_ptr<CausalState> getCausalState(Distribution dist, FloatType matchThreshold = 0, bool store = true) {
        if (!distributionMatcherCallback) return boost::shared_ptr<CausalState>();
        return getCausalState(dist, matchThreshold, store, distributionMatcherCallback);
    }

    template <class DistributionMatcher>
    boost::shared_ptr<CausalState> getCausalState(Distribution dist, FloatType matchThreshold, bool store, DistributionMatcher distributionMatcher) {
        // See the logic in getIsoUtilityState, some comments are not repeated here
        StateInfo* info = static_cast<StateInfo*>(helpers::StateInfoHelper<DistributionManager>::getInfo(distributionManager,dist));
        if (info && info->causalState) return info->causalState;
        // clustering done and no causal state => necessarily an "extra" info if there is an info
        if (store && !info) {
            extraInfos.push_back(StateInfo());
            info = &extraInfos.back();
            info->distribution = dist;
            helpers::StateInfoHelper<DistributionManager>::setInfo(distributionManager,dist,info);
        }
        boost::shared_ptr<CausalState> result;
        FloatType bestMatch = -std::numeric_limits<FloatType>::max();
        for (typename CausalStates::iterator it = causalStates.begin(); it != causalStates.end(); ++it) {
            FloatType match = distributionMatcher(dist,(*it)->distribution);
            if (match > bestMatch) {
                bestMatch = match;
                result = *it;
            }
        }
        if (result && bestMatch > matchThreshold) {
            if (info && store) info->causalState = result;
            return result;
        }
        result = boost::shared_ptr<CausalState>(new CausalState());
        result->distribution = dist;
        result->count = 1;
        if (info && store) {
            info->causalState = result;
            result->members->push_back(info);
        }
        return result;
    }

    boost::shared_ptr<IsoUtilityState> getIsoUtilityState(DataType x, FloatType matchThreshold = 0, bool store = true) {
        // No info, compute the utility for that new point
        if (!utilityFunctorCallback) return boost::shared_ptr<IsoUtilityState>();  // no utility function given, cannot compute state
        return getIsoUtilityState(x, matchThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    boost::shared_ptr<IsoUtilityState> getIsoUtilityState(DataType x, FloatType matchThreshold, bool store, UtilityFunctor u) {
        if (!utilityMatcherCallback) return boost::shared_ptr<IsoUtilityState>();
        return getIsoUtilityState(x, matchThreshold, store, u, utilityMatcherCallback);
    }

    template <class UtilityFunctor, class UtilityMatcher>
    boost::shared_ptr<IsoUtilityState> getIsoUtilityState(DataType x, FloatType matchThreshold, bool store, UtilityFunctor utilityFunctor, UtilityMatcher utilityMatcher) {
        Distribution dist = distributionManager.distribution(x);
        return getIsoUtilityState(dist, matchThreshold, store, utilityFunctor, utilityMatcher);
    }

    boost::shared_ptr<IsoUtilityState> getIsoUtilityState(Distribution dist, FloatType matchThreshold = 0, bool store = true) {
        // No info, compute the utility for that new point
        if (!utilityFunctorCallback) return boost::shared_ptr<IsoUtilityState>();  // no utility function given, cannot compute state
        return getIsoUtilityState(dist, matchThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    boost::shared_ptr<IsoUtilityState> getIsoUtilityState(Distribution dist, FloatType matchThreshold, bool store, UtilityFunctor u) {
        if (!utilityMatcherCallback) return boost::shared_ptr<IsoUtilityState>();
        return getIsoUtilityState(dist, matchThreshold, store, u, utilityMatcherCallback);
    }

    template <class UtilityFunctor, class UtilityMatcher>
    boost::shared_ptr<IsoUtilityState> getIsoUtilityState(Distribution dist, FloatType matchThreshold, bool store, UtilityFunctor utilityFunctor, UtilityMatcher utilityMatcher) {
        StateInfo* info = static_cast<StateInfo*>(helpers::StateInfoHelper<DistributionManager>::getInfo(distributionManager,dist));
        if (info && info->isoUtilityState) return info->isoUtilityState;
        PredictionSet predictionSet;
        FloatType utility;
        // if isoPredictionState exists then necessarily utility and predictionSet are valid
        // as they are both set at the same time (whether in the main object or as extra info)
        if (info && info->isoPredictionState) {
            utility = info->utility;
            predictionSet = info->predictionSet;
        // otherwise we need to compute them both, and possibly store them for later use
        } else {
            utility = optimiser->optimise(ExpectedUtilityFunctor<UtilityFunctor>(distributionManager, dist, utilityFunctor), predictionSet);
            if (store) {
                // this info object is special as it does not correspond to one of the original points used for clustering
                // keep it separately so as not to mess with the clustering algorithms
                extraInfos.push_back(StateInfo());
                // list => pointers to elements are not invalidated on insert
                info = &extraInfos.back();
                info->distribution = dist;
                info->utility = utility;
                info->predictionSet = predictionSet;
                // getInfo failed above so we can safely ask the distribution manager to store if possible without conflict
                helpers::StateInfoHelper<DistributionManager>::setInfo(distributionManager,dist,info);
            }
        }
        boost::shared_ptr<IsoUtilityState> result;
        FloatType bestMatch = -std::numeric_limits<FloatType>::max();
        // Use the matcher with existing clusters aggregate values to find a matching state, if any
        for (typename IsoUtilityStates::iterator it = isoUtilityStates.begin(); it != isoUtilityStates.end(); ++it) {
            FloatType match = utilityMatcher(utility,(*it)->utility);
            if (match > bestMatch) {
                bestMatch = match;
                result = *it;
            }
        }
        if (result && bestMatch > matchThreshold) {
            if (info && store) info->isoUtilityState = result;
            return result;
        }
        // No match. Create a temporary state so the user can access the computed utility for the given x point
        result = boost::shared_ptr<IsoUtilityState>(new IsoUtilityState());
        result->utility = utility;
        result->count = 1;
        if (info && store) {
            info->isoUtilityState = result;
            result->members->push_back(info); // so user can access distribution, etc.
        }
        return result;
    }

    boost::shared_ptr<IsoPredictionState> getIsoPredictionState(DataType x, FloatType matchThreshold = 0, bool store = true) {
        // No info, compute the utility for that new point
        if (!utilityFunctorCallback) return boost::shared_ptr<IsoPredictionState>();  // no utility function given, cannot compute state
        return getIsoPredictionState(x, matchThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    boost::shared_ptr<IsoPredictionState> getIsoPredictionState(DataType x, FloatType matchThreshold, bool store, UtilityFunctor u) {
        if (!predictionSetMatcherCallback) return boost::shared_ptr<IsoPredictionState>();
        return getIsoPredictionState(x, matchThreshold, store, u, predictionSetMatcherCallback);
    }

    template <class UtilityFunctor, class PredictionMatcher>
    boost::shared_ptr<IsoPredictionState> getIsoPredictionState(DataType x, FloatType matchThreshold, bool store, UtilityFunctor utilityFunctor, PredictionMatcher predictionMatcher) {
        Distribution dist = distributionManager.distribution(x);
        boost::shared_ptr<IsoPredictionState> ret = getIsoPredictionState(dist, matchThreshold, store, utilityFunctor, predictionSetMatcherCallback);
        return ret;
    }
    
    boost::shared_ptr<IsoPredictionState> getIsoPredictionState(Distribution dist, FloatType matchThreshold = 0, bool store = true) {
        // No info, compute the utility for that new point
        if (!utilityFunctorCallback) return boost::shared_ptr<IsoPredictionState>();  // no utility function given, cannot compute state
        return getIsoPredictionState(dist, matchThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    boost::shared_ptr<IsoPredictionState> getIsoPredictionState(Distribution dist, FloatType matchThreshold, bool store, UtilityFunctor u) {
        if (!predictionSetMatcherCallback) return boost::shared_ptr<IsoPredictionState>();
        return getIsoPredictionState(dist, matchThreshold, store, u, predictionSetMatcherCallback);
    }
    
    
    template <class UtilityFunctor, class PredictionMatcher>
    boost::shared_ptr<IsoPredictionState> getIsoPredictionState(Distribution dist, FloatType matchThreshold, bool store, UtilityFunctor utilityFunctor, PredictionMatcher predictionMatcher) {
        // See the logic in getIsoUtilityState, comments are not repeated here
        // the store parameter allows to avoid repeating computations
        StateInfo* info = static_cast<StateInfo*>(helpers::StateInfoHelper<DistributionManager>::getInfo(distributionManager,dist));
        if (info && info->isoPredictionState) return info->isoPredictionState;
        PredictionSet predictionSet;
        FloatType utility;
        if (info && info->isoUtilityState) {
            utility = info->utility;
            predictionSet = info->predictionSet;
        } else {
            utility = optimiser->optimise(ExpectedUtilityFunctor<UtilityFunctor>(distributionManager, dist, utilityFunctor), predictionSet);
            if (store) {
                extraInfos.push_back(StateInfo());
                info = &extraInfos.back();
                info->distribution = dist;
                info->utility = utility;
                info->predictionSet = predictionSet;
                helpers::StateInfoHelper<DistributionManager>::setInfo(distributionManager,dist,info);
            }
        }
        boost::shared_ptr<IsoPredictionState> result;
        FloatType bestMatch = -std::numeric_limits<FloatType>::max();
        for (typename IsoPredictionStates::iterator it = isoPredictionStates.begin(); it != isoPredictionStates.end(); ++it) {
            FloatType match = predictionMatcher(predictionSet,(*it)->predictionSet);
            if (match > bestMatch) {
                bestMatch = match;
                result = *it;
            }
        }
        if (result && bestMatch > matchThreshold) {
            if (info && store) info->isoPredictionState = result;
            return result;
        }
        // No match. Create a temporary state so the user can access the computed utility for the given x point
        result = boost::shared_ptr<IsoPredictionState>(new IsoPredictionState());
        result->predictionSet = predictionSet;
        result->count = 1;
        if (info && store) {
            info->isoPredictionState = result;
            result->members->push_back(info); // so user can access distribution, etc.
        }
        return result;
    }

    // take first match, see below. TODO: doc
    boost::shared_ptr<DecisionalState> getDecisionalState(DataType x, FloatType matchUtilityThreshold = 0, FloatType matchPredictionThreshold = 0, bool store = true) {
        // No info, compute the utility for that new point
        if (!utilityFunctorCallback) return boost::shared_ptr<DecisionalState>();  // no utility function given, cannot compute state
        return getDecisionalState(x, matchUtilityThreshold, matchPredictionThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    boost::shared_ptr<DecisionalState> getDecisionalState(DataType x, FloatType matchUtilityThreshold, FloatType matchPredictionThreshold, bool store, UtilityFunctor u) {
        if (!utilityMatcherCallback) return boost::shared_ptr<DecisionalState>();
        if (!predictionSetMatcherCallback) return boost::shared_ptr<DecisionalState>();
        return getDecisionalState(x, matchUtilityThreshold, matchPredictionThreshold, store, u, utilityMatcherCallback, predictionSetMatcherCallback);
    }

    template <class UtilityFunctor, class UtilityMatcher, class PredictionMatcher>
    boost::shared_ptr<DecisionalState> getDecisionalState(DataType x, FloatType matchUtilityThreshold, FloatType matchPredictionThreshold, bool store, UtilityFunctor u, UtilityMatcher utilityMatcher, PredictionMatcher predictionMatcher) {
        Distribution dist = distributionManager.distribution(x);
        return getDecisionalState(dist, matchUtilityThreshold, matchPredictionThreshold, store, u, utilityMatcher, predictionMatcher);
    }

    boost::shared_ptr<DecisionalState> getDecisionalState(Distribution dist, FloatType matchUtilityThreshold = 0, FloatType matchPredictionThreshold = 0, bool store = true) {
        // No info, compute the utility for that new point
        if (!utilityFunctorCallback) return boost::shared_ptr<DecisionalState>();  // no utility function given, cannot compute state
        return getDecisionalState(dist, matchUtilityThreshold, matchPredictionThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    boost::shared_ptr<DecisionalState> getDecisionalState(Distribution dist, FloatType matchUtilityThreshold, FloatType matchPredictionThreshold, bool store, UtilityFunctor u) {
        if (!utilityMatcherCallback) return boost::shared_ptr<DecisionalState>();
        if (!predictionSetMatcherCallback) return boost::shared_ptr<DecisionalState>();
        // This line triggers a spurious warning with some versions of gcc-4.4
        // The warning can be ignored, see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=40146
        return getDecisionalState(dist, matchUtilityThreshold, matchPredictionThreshold, store, u, utilityMatcherCallback, predictionSetMatcherCallback);
    }

    template <class UtilityFunctor, class UtilityMatcher, class PredictionMatcher>
    boost::shared_ptr<DecisionalState> getDecisionalState(Distribution dist, FloatType matchUtilityThreshold, FloatType matchPredictionThreshold, bool store, UtilityFunctor u, UtilityMatcher utilityMatcher, PredictionMatcher predictionMatcher) {
        // See the logic in getIsoUtilityState, comments are not repeated here
        // the store parameter allows to avoid repeating computations
        StateInfo* info = static_cast<StateInfo*>(helpers::StateInfoHelper<DistributionManager>::getInfo(distributionManager,dist));
        if (info && info->decisionalState) return info->decisionalState;
        PredictionSet predictionSet;
        FloatType utility;
        if (info && (info->isoUtilityState || info->isoPredictionState)) {
            utility = info->utility;
            predictionSet = info->predictionSet;
        } else {
            utility = optimiser->optimise(ExpectedUtilityFunctor<UtilityFunctor>(distributionManager, dist, u), predictionSet);
            if (store) {
                extraInfos.push_back(StateInfo());
                info = &extraInfos.back();
                info->distribution = dist;
                info->utility = utility;
                info->predictionSet = predictionSet;
                helpers::StateInfoHelper<DistributionManager>::setInfo(distributionManager,dist,info);
            }
        }
        // take first match, best match would require weighting between both thresholds on the pareto front
        for (typename DecisionalStates::iterator it = decisionalStates.begin(); it != decisionalStates.end(); ++it) {
            if ( ( (info && ((*it)->isoUtilityState == info->isoUtilityState)) || (utilityMatcher(utility,(*it)->isoUtilityState->utility) > matchUtilityThreshold))
              && ( (info && ((*it)->isoPredictionState == info->isoPredictionState)) || (predictionMatcher(predictionSet,(*it)->isoPredictionState->predictionSet))>matchPredictionThreshold)
            ) {
                if (info) { // extra info object was created
                    info->isoPredictionState = (*it)->isoPredictionState;
                    info->isoUtilityState = (*it)->isoUtilityState;
                    info->decisionalState = *it;
                }
                return *it;
            }
        }
        boost::shared_ptr<DecisionalState> result = boost::shared_ptr<DecisionalState>(new DecisionalState());
        result->isoUtilityState = boost::shared_ptr<IsoUtilityState>(new IsoUtilityState());
        result->isoUtilityState->utility = utility;
        result->isoUtilityState->count = 1;
        result->isoPredictionState = boost::shared_ptr<IsoPredictionState>(new IsoPredictionState());
        result->isoPredictionState->predictionSet = predictionSet;
        result->isoPredictionState->count = 1;
        result->count = 1;
        if (info && store) {
            info->isoUtilityState = result->isoUtilityState;
            info->isoPredictionState = result->isoPredictionState;
            info->decisionalState = result;
            result->members = boost::shared_ptr<std::vector<StateInfo*> >(new std::vector<StateInfo*>());
            result->members->push_back(info); // so user can access distribution, etc.
            result->isoUtilityState->members = boost::shared_ptr<std::vector<StateInfo*> >(new std::vector<StateInfo*>());
            result->isoUtilityState->members->push_back(info);
            result->isoPredictionState->members = boost::shared_ptr<std::vector<StateInfo*> >(new std::vector<StateInfo*>());
            result->isoPredictionState->members->push_back(info);
        }
        return result;
    }

    // statistics / frequency of each state. Total count is dataset.size()
    struct StateCounts {
        int mincount, maxcount;
        int totalcount;
        int recurrentcount;
        int nrecurrent;
        StateCounts() {reset();}
        void reset() {
            mincount=std::numeric_limits<int>::max();
            maxcount=-std::numeric_limits<int>::max();
            totalcount=0;
            recurrentcount=0;
            nrecurrent=0;
        }
        void update(int count, bool is_recurrent) {
            if (count<mincount) mincount = count;
            if (count>maxcount) maxcount = count;
            totalcount += count;
            if (is_recurrent) {
                recurrentcount += count;
                ++nrecurrent;
            }
        };
    };
    
    StateCounts causalStatesCounts, isoUtilityStatesCounts, isoPredictionStatesCounts, decisionalStatesCounts;
    
    template<class TFder, typename SType>
    struct SymbolKeyMaker {
        static CausalStateTransition getTransition(TFder& tfder, typename TFder::iterator it, CausalState* state) {
            CausalStateTransition ret;
            ret.label = tfder.getSymbol(it);
            ret.state = state;
            ret.probability = 1;
            return ret;
        }
        enum {HasSymbol = 1};
    };
    template<class TFder>
    struct SymbolKeyMaker<TFder, void> {
        static CausalStateTransition getTransition(TFder& tfder, typename TFder::iterator it, CausalState* state) {
            CausalStateTransition ret;
            ret.state = state;
            ret.probability = 1;
            return ret;
        }
        enum {HasSymbol = 0};
    };

    template<class State, class TransitionSet>
    void markRecurrentStates(std::vector<boost::shared_ptr<State> >& states, FloatType leakThreshold) {
        // compute root nodes - temporarily use the flags as indegree count
        for (int ci=0; ci<(int)states.size(); ++ci) states[ci]->flags = 0;
        bool newRootsFound = false;
        int passcount = 0;
        do {
            for (int ci=0; ci<(int)states.size(); ++ci) {
                if (states[ci]->flags<0) continue; // ignore previous roots, seek new ones
                for (typename TransitionSet::iterator it = states[ci]->transitions->begin(); it != states[ci]->transitions->end(); ++it) {
                    // do not count self-refs
                    if (it->state!=states[ci].get()) ++it->state->flags;
                    
                }
            }
            // roots are the nodes with no parent
            newRootsFound = false;
            for (int ci=0; ci<(int)states.size(); ++ci) {
                if (states[ci]->flags<0) continue; // ignore previous roots, seek new ones
                if (states[ci]->flags==0) {
                    if (passcount) states[ci]->flags = -1;
                    else states[ci]->flags = -2; // first loop = real roots
                    newRootsFound = true;
                }
                else states[ci]->flags = 0; // reset flags for next pass
            }
            ++passcount;
        } while (newRootsFound);
        for (int ci=0; ci<(int)states.size(); ++ci) {
            if (states[ci]->flags==-2) states[ci]->flags = ROOT_STATE;
            else if (states[ci]->flags>=0) states[ci]->flags = RECURRENT_STATE;
            else states[ci]->flags = 0;
        }
    }

    /// epsilon-machine when the data set provides symbols
    /// otherwise simply an unlabelled automaton
    void buildCausalStateGraph(FloatType leakThreshold = 0, bool deterministic = true) {
        // create transition objects, guaranteed to exist (possibly empty) after returning from this function
        for (typename CausalStates::iterator cit = causalStates.begin(); cit != causalStates.end(); ++cit) {
            (*cit)->transitions = boost::shared_ptr<CausalStateTransitionSet>(new CausalStateTransitionSet);
        }
        for (typename TransitionFeeder::iterator tit = transitionFeeder->begin(); tit != transitionFeeder->end(); ++tit) {
            DataType prevx = transitionFeeder->getDataBeforeTransition(tit);
            StateInfo* info = (StateInfo*)distributionManager.getInfo(distributionManager.distribution(prevx));
            if (!info) continue;
            CausalState* state = info->causalState.get();
            if (!state) continue;
            DataType x = transitionFeeder->getDataAfterTransition(tit);
            info = (StateInfo*)distributionManager.getInfo(distributionManager.distribution(x));
            if (!info) continue;
            CausalState* nextstate = info->causalState.get();
            if (!nextstate) continue;
            CausalStateTransition transition = SymbolKeyMaker<TransitionFeeder, SymbolType>::getTransition(*transitionFeeder,tit,nextstate);
            typename CausalStateTransitionSet::iterator it = state->transitions->find(transition.getKey());
            if (it!=state->transitions->end()) {
                it->probability += 1; // increment total for this symbol
                it->incrementCountForState(transition.state);
            }
            else {
                // proba set to one in SymbolKeyMaker
                it = state->transitions->insert(transition).first;
                it->initCountForState(transition.state);
            }
        }
        // convert counts to probabilities
        for (typename CausalStates::iterator cit = causalStates.begin(); cit != causalStates.end(); ++cit) {
            FloatType total = 0;
            for (typename CausalStateTransitionSet::iterator it = (*cit)->transitions->begin(); it != (*cit)->transitions->end(); ++it) {
                if (deterministic) it->retainMaxCountState();
                total += it->probability;  // probability that the symbol is taken, irrespectively of retained state
            }
            if (total==0) continue; // no transition
            for (typename CausalStateTransitionSet::iterator it = (*cit)->transitions->begin(); it != (*cit)->transitions->end(); ++it) {
                it->probability /= total;
                if (!deterministic) it->rescaleStatesProba(total);
            }
        }

        markRecurrentStates<CausalState, CausalStateTransitionSet>(causalStates, leakThreshold);

        // state counts for proba
        causalStatesCounts.reset();
        for (typename CausalStates::iterator it = causalStates.begin(); it != causalStates.end(); ++it) (*it)->count = 0;
        for (typename DistributionManager::iterator it = distributionManager.begin(); it != distributionManager.end(); ++it) {
            Distribution& dist = helpers::DistributionManagerIteratorHolder<DistributionManager>(distributionManager,it).deref();
            StateInfo* info = static_cast<StateInfo*>(helpers::StateInfoHelper<DistributionManager>::getInfo(distributionManager, dist));
            assert(info!=0);
            int count = distributionManager.getCount(dist);
            if (info->causalState) info->causalState->count += count;
        }
        for (typename CausalStates::iterator it = causalStates.begin(); it != causalStates.end(); ++it) causalStatesCounts.update((*it)->count, (*it)->is_recurrent());
    }

    // TODO: generic template for the 3 others to avoid code dup.
    // But needs member pointer to info struct, etc, too much trouble for now

    /// Unlabeled iso-prediction states graph
    void buildIsoPredictionStateGraph(FloatType leakThreshold = 0) {
        // create transition objects, guaranteed to exist (possibly empty) after returning from this function
        for (typename IsoPredictionStates::iterator cit = isoPredictionStates.begin(); cit != isoPredictionStates.end(); ++cit) {
            (*cit)->transitions = boost::shared_ptr<IsoPredictionStateTransitionSet>(new IsoPredictionStateTransitionSet);
        }
        for (typename TransitionFeeder::iterator tit = transitionFeeder->begin(); tit != transitionFeeder->end(); ++tit) {
            DataType prevx = transitionFeeder->getDataBeforeTransition(tit);
            StateInfo* info = (StateInfo*)distributionManager.getInfo(distributionManager.distribution(prevx));
            if (!info) continue;
            IsoPredictionState* state = info->isoPredictionState.get();
            if (!state) continue;
            DataType x = transitionFeeder->getDataAfterTransition(tit);
            info = (StateInfo*)distributionManager.getInfo(distributionManager.distribution(x));
            if (!info) continue;
            IsoPredictionState* nextstate = info->isoPredictionState.get();
            if (!nextstate) continue;
            // keys are here the next state, not symbols
            typename IsoPredictionStateTransitionSet::iterator it = state->transitions->find(nextstate);
            if (it!=state->transitions->end()) {
                it->probability += 1; // increment total for this state
            }
            else {
                it = state->transitions->insert(IsoPredictionStateTransition(nextstate)).first;
                it->probability = 1; // increment total for this state
            }
        }
        // convert counts to probabilities
        for (typename IsoPredictionStates::iterator cit = isoPredictionStates.begin(); cit != isoPredictionStates.end(); ++cit) {
            FloatType total = 0;
            for (typename IsoPredictionStateTransitionSet::iterator it = (*cit)->transitions->begin(); it != (*cit)->transitions->end(); ++it) {
                total += it->probability;
            }
            if (total==0) continue; // no transition
            for (typename IsoPredictionStateTransitionSet::iterator it = (*cit)->transitions->begin(); it != (*cit)->transitions->end(); ++it) {
                it->probability /= total;
            }
        }

        markRecurrentStates<IsoPredictionState, IsoPredictionStateTransitionSet>(isoPredictionStates, leakThreshold);
        
        // state counts for proba
        isoPredictionStatesCounts.reset();
        for (typename IsoPredictionStates::iterator it = isoPredictionStates.begin(); it != isoPredictionStates.end(); ++it) (*it)->count = 0;
        for (typename DistributionManager::iterator it = distributionManager.begin(); it != distributionManager.end(); ++it) {
            Distribution& dist = helpers::DistributionManagerIteratorHolder<DistributionManager>(distributionManager,it).deref();
            StateInfo* info = static_cast<StateInfo*>(helpers::StateInfoHelper<DistributionManager>::getInfo(distributionManager, dist));
            assert(info!=0);
            int count = distributionManager.getCount(dist);
            if (info->isoPredictionState) info->isoPredictionState->count += count;
        }
        for (typename IsoPredictionStates::iterator it = isoPredictionStates.begin(); it != isoPredictionStates.end(); ++it) isoPredictionStatesCounts.update((*it)->count, (*it)->is_recurrent());
    }

    /// Unlabelled iso-utility states graph
    void buildIsoUtilityStateGraph(FloatType leakThreshold = 0) {
        // create transition objects, guaranteed to exist (possibly empty) after returning from this function
        for (typename IsoUtilityStates::iterator cit = isoUtilityStates.begin(); cit != isoUtilityStates.end(); ++cit) {
            (*cit)->transitions = boost::shared_ptr<IsoUtilityStateTransitionSet>(new IsoUtilityStateTransitionSet);
        }
        for (typename TransitionFeeder::iterator tit = transitionFeeder->begin(); tit != transitionFeeder->end(); ++tit) {
            DataType prevx = transitionFeeder->getDataBeforeTransition(tit);
            StateInfo* info = (StateInfo*)distributionManager.getInfo(distributionManager.distribution(prevx));
            if (!info) continue;
            IsoUtilityState* state = info->isoUtilityState.get();
            if (!state) continue;
            DataType x = transitionFeeder->getDataAfterTransition(tit);
            info = (StateInfo*)distributionManager.getInfo(distributionManager.distribution(x));
            if (!info) continue;
            IsoUtilityState* nextstate = info->isoUtilityState.get();
            if (!nextstate) continue;
            // keys are here the next state, not symbols
            typename IsoUtilityStateTransitionSet::iterator it = state->transitions->find(nextstate);
            if (it!=state->transitions->end()) {
                it->probability += 1; // increment total for this state
            }
            else {
                it = state->transitions->insert(IsoUtilityStateTransition(nextstate)).first;
                it->probability = 1; // increment total for this state
            }
        }
        // convert counts to probabilities
        for (typename IsoUtilityStates::iterator cit = isoUtilityStates.begin(); cit != isoUtilityStates.end(); ++cit) {
            FloatType total = 0;
            for (typename IsoUtilityStateTransitionSet::iterator it = (*cit)->transitions->begin(); it != (*cit)->transitions->end(); ++it) {
                total += it->probability;
            }
            if (total==0) continue; // no transition
            for (typename IsoUtilityStateTransitionSet::iterator it = (*cit)->transitions->begin(); it != (*cit)->transitions->end(); ++it) {
                it->probability /= total;
            }
        }

        markRecurrentStates<IsoUtilityState, IsoUtilityStateTransitionSet>(isoUtilityStates, leakThreshold);
        
        // state counts for proba
        isoUtilityStatesCounts.reset();
        for (typename IsoUtilityStates::iterator it = isoUtilityStates.begin(); it != isoUtilityStates.end(); ++it) (*it)->count = 0;
        for (typename DistributionManager::iterator it = distributionManager.begin(); it != distributionManager.end(); ++it) {
            Distribution& dist = helpers::DistributionManagerIteratorHolder<DistributionManager>(distributionManager,it).deref();
            StateInfo* info = static_cast<StateInfo*>(helpers::StateInfoHelper<DistributionManager>::getInfo(distributionManager, dist));
            assert(info!=0);
            int count = distributionManager.getCount(dist);
            if (info->isoUtilityState) info->isoUtilityState->count += count;
        }
        for (typename IsoUtilityStates::iterator it = isoUtilityStates.begin(); it != isoUtilityStates.end(); ++it) isoUtilityStatesCounts.update((*it)->count, (*it)->is_recurrent());
    }

    /// Unlabelled decisional states graph
    void buildDecisionalStateGraph(FloatType leakThreshold = 0) {
        // create transition objects, guaranteed to exist (possibly empty) after returning from this function
        for (typename DecisionalStates::iterator cit = decisionalStates.begin(); cit != decisionalStates.end(); ++cit) {
            (*cit)->transitions = boost::shared_ptr<DecisionalStateTransitionSet>(new DecisionalStateTransitionSet);
        }
        for (typename TransitionFeeder::iterator tit = transitionFeeder->begin(); tit != transitionFeeder->end(); ++tit) {
            DataType prevx = transitionFeeder->getDataBeforeTransition(tit);
            StateInfo* info = (StateInfo*)distributionManager.getInfo(distributionManager.distribution(prevx));
            if (!info) continue;
            DecisionalState* state = info->decisionalState.get();
            if (!state) continue;
            DataType x = transitionFeeder->getDataAfterTransition(tit);
            info = (StateInfo*)distributionManager.getInfo(distributionManager.distribution(x));
            if (!info) continue;
            DecisionalState* nextstate = info->decisionalState.get();
            if (!nextstate) continue;
            // keys are here the next state, not symbols
            typename DecisionalStateTransitionSet::iterator it = state->transitions->find(nextstate);
            if (it!=state->transitions->end()) {
                it->probability += 1; // increment total for this state
            }
            else {
                it = state->transitions->insert(DecisionalStateTransition(nextstate)).first;
                it->probability = 1; // increment total for this state
            }
        }
        // convert counts to probabilities
        for (typename DecisionalStates::iterator cit = decisionalStates.begin(); cit != decisionalStates.end(); ++cit) {
            FloatType total = 0;
            for (typename DecisionalStateTransitionSet::iterator it = (*cit)->transitions->begin(); it != (*cit)->transitions->end(); ++it) {
                total += it->probability;
            }
            if (total==0) continue; // no transition
            for (typename DecisionalStateTransitionSet::iterator it = (*cit)->transitions->begin(); it != (*cit)->transitions->end(); ++it) {
                it->probability /= total;
            }
        }

        markRecurrentStates<DecisionalState, DecisionalStateTransitionSet>(decisionalStates, leakThreshold);
        
        // state counts for proba
        decisionalStatesCounts.reset();
        for (typename DecisionalStates::iterator it = decisionalStates.begin(); it != decisionalStates.end(); ++it) (*it)->count = 0;
        for (typename DistributionManager::iterator it = distributionManager.begin(); it != distributionManager.end(); ++it) {
            Distribution& dist = helpers::DistributionManagerIteratorHolder<DistributionManager>(distributionManager,it).deref();
            StateInfo* info = static_cast<StateInfo*>(helpers::StateInfoHelper<DistributionManager>::getInfo(distributionManager, dist));
            assert(info!=0);
            int count = distributionManager.getCount(dist);
            if (info->decisionalState) info->decisionalState->count += count;
        }
        for (typename DecisionalStates::iterator it = decisionalStates.begin(); it != decisionalStates.end(); ++it) decisionalStatesCounts.update((*it)->count, (*it)->is_recurrent());
    }

    template<typename LabelType, int haslabel>
    struct GetLabelHelper {
        static LabelType getLabel(typename CausalStateTransitionSet::iterator it) {
            return it->label;
        }
    };
    template<typename LabelType>
    struct GetLabelHelper<LabelType, 0> {
        static std::string getLabel(typename CausalStateTransitionSet::iterator it) {
            return "";
        }
    };

    // exporting in graphviz format
    void writeCausalStateGraph(std::ostream& dot, bool recurrent_only = true) {
        dot << "digraph G {" << std::endl;
        dot << "rankdir=LR;" << std::endl;
        dot << "node [shape = circle];" << std::endl;
        dot.precision(3);
        FloatType totalcount = recurrent_only ? causalStatesCounts.recurrentcount : causalStatesCounts.totalcount;
        for (typename CausalStates::iterator cit = causalStates.begin(); cit != causalStates.end(); ++cit) {
            if (recurrent_only && !(*cit)->is_recurrent()) continue;
            FloatType stateProba = (*cit)->count / totalcount;
            dot << (*cit)->index << " [label=\"p=" << stateProba << "\"]" << std::endl;
            for (typename CausalStateTransitionSet::iterator it = (*cit)->transitions->begin(); it != (*cit)->transitions->end(); ++it) {
                dot << (*cit)->index <<  " -> " << it->state->index << " [label=\"";
                if (SymbolKeyMaker<TransitionFeeder, SymbolType>::HasSymbol) {
                    dot << "s=" << GetLabelHelper<SymbolType, SymbolKeyMaker<TransitionFeeder, SymbolType>::HasSymbol>::getLabel(it) << ", ";
                }
                dot << "p=" << it->probability << "\"];" << std::endl;
            }
        }
        dot << "}" << std::endl;
    }

    // template is easier to do for the graphviz writer
    template<class UnlabeledStates, class UnlabeledStateTransitionSet>
    void writeUnlabeledStateGraph(UnlabeledStates& unlabeledStates, StateCounts& unlabeledStatesCounts, std::ostream& dot, bool recurrent_only) {
        dot << "digraph G {" << std::endl;
        dot << "rankdir=LR;" << std::endl;
        dot << "node [shape = circle];" << std::endl;
        dot.precision(3);
        FloatType totalcount = recurrent_only ? unlabeledStatesCounts.recurrentcount : unlabeledStatesCounts.totalcount;
        for (typename UnlabeledStates::iterator cit = unlabeledStates.begin(); cit != unlabeledStates.end(); ++cit) {
            if (recurrent_only && !(*cit)->is_recurrent()) continue;
            FloatType stateProba = (*cit)->count / totalcount;
            dot << (*cit)->index << " [label=\"p=" << stateProba << "\"]" << std::endl;
            for (typename UnlabeledStateTransitionSet::iterator it = (*cit)->transitions->begin(); it != (*cit)->transitions->end(); ++it) {
                dot << (*cit)->index <<  " -> " << it->state->index << " [label=\"";
                dot << "p=" << it->probability << "\"];" << std::endl;
            }
        }
        dot << "}" << std::endl;
    }
    
    // exporting in graphviz format
    void writeIsoPredictionStateGraph(std::ostream& dot, bool recurrent_only = true) {
        writeUnlabeledStateGraph<IsoPredictionStates, IsoPredictionStateTransitionSet>(isoPredictionStates, isoPredictionStatesCounts, dot, recurrent_only);
    }
    void writeIsoUtilityStateGraph(std::ostream& dot, bool recurrent_only = true) {
        writeUnlabeledStateGraph<IsoUtilityStates, IsoUtilityStateTransitionSet>(isoUtilityStates, isoUtilityStatesCounts, dot, recurrent_only);
    }
    void writeDecisionalStateGraph(std::ostream& dot, bool recurrent_only = true) {
        writeUnlabeledStateGraph<DecisionalStates, DecisionalStateTransitionSet>(decisionalStates, decisionalStatesCounts, dot, recurrent_only);
    }

    template<class States, class State>
    FloatType computeComplexity(States& states, StateCounts& statesCounts, bool recurrent_only) {
        FloatType C = 0;
        double totalCount = recurrent_only ? (double)statesCounts.recurrentcount : (double)statesCounts.totalcount;
        for (typename States::iterator cit = states.begin(); cit != states.end(); ++cit) {
            if (recurrent_only && !(*cit)->is_recurrent()) continue;
            State* currentState = cit->get();
            double stateProba = currentState->count / totalCount;
            if (stateProba) C -= stateProba * log2(stateProba);
        }
        return C;
    }

    FloatType statisticalComplexity(bool recurrent_only = true) {
        return computeComplexity<CausalStates, CausalState>(causalStates, causalStatesCounts, recurrent_only);
    }
    FloatType isoPredictionComplexity(bool recurrent_only = true) {
        return computeComplexity<IsoPredictionStates, IsoPredictionState>(isoPredictionStates, isoPredictionStatesCounts, recurrent_only);
    }
    FloatType isoUtilityComplexity(bool recurrent_only = true) {
        return computeComplexity<IsoUtilityStates, IsoUtilityState>(isoUtilityStates, isoUtilityStatesCounts, recurrent_only);
    }
    FloatType decisionalComplexity(bool recurrent_only = true) {
        return computeComplexity<DecisionalStates, DecisionalState>(decisionalStates, decisionalStatesCounts, recurrent_only);
    }

/// LOCAL Complexity fonctions. The actual code is very simple,
/// but the API covers all cases for generalizing when the data was not met before

    FloatType localStatisticalComplexity(DataType x, FloatType matchThreshold = 0, bool store = true) {
        if (helpers::IsMissingHelper<TransitionFeeder,DataType>::isMissing(*transitionFeeder,x)) return std::numeric_limits<FloatType>::quiet_NaN();
        if (!distributionMatcherCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        Distribution dist = distributionManager.distribution(x);
        return localStatisticalComplexity(dist, matchThreshold, store, distributionMatcherCallback);
    }

    FloatType localStatisticalComplexity(Distribution dist, FloatType matchThreshold = 0, bool store = true) {
        if (!distributionMatcherCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localStatisticalComplexity(dist, matchThreshold, store, distributionMatcherCallback);
    }

    template <class DistributionMatcher>
    FloatType localStatisticalComplexity(Distribution dist, FloatType matchThreshold, bool store, DistributionMatcher distributionMatcher) {
        boost::shared_ptr<CausalState> state = getCausalState(dist, matchThreshold, store, distributionMatcherCallback);
        return -helpers::log2((FloatType)state->count/causalStatesCounts.totalcount);
    }

    FloatType localIsoUtilityComplexity(DataType x, FloatType matchThreshold = 0, bool store = true) {
        // No info, compute the utility for that new point
        if (!utilityFunctorCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localIsoUtilityComplexity(x, matchThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    FloatType localIsoUtilityComplexity(DataType x, FloatType matchThreshold, bool store, UtilityFunctor u) {
        if (!utilityMatcherCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localIsoUtilityComplexity(x, matchThreshold, store, u, utilityMatcherCallback);
    }

    template <class UtilityFunctor, class UtilityMatcher>
    FloatType localIsoUtilityComplexity(DataType x, FloatType matchThreshold, bool store, UtilityFunctor utilityFunctor, UtilityMatcher utilityMatcher) {
        Distribution dist = distributionManager.distribution(x);
        return localIsoUtilityComplexity(dist, matchThreshold, store, utilityFunctor, utilityMatcher);
    }

    FloatType localIsoUtilityComplexity(Distribution dist, FloatType matchThreshold = 0, bool store = true) {
        // No info, compute the utility for that new point
        if (!utilityFunctorCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localIsoUtilityComplexity(dist, matchThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    FloatType localIsoUtilityComplexity(Distribution dist, FloatType matchThreshold, bool store, UtilityFunctor u) {
        if (!utilityMatcherCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localIsoUtilityComplexity(dist, matchThreshold, store, u, utilityMatcherCallback);
    }

    template <class UtilityFunctor, class UtilityMatcher>
    FloatType localIsoUtilityComplexity(Distribution dist, FloatType matchThreshold, bool store, UtilityFunctor utilityFunctor, UtilityMatcher utilityMatcher) {
        boost::shared_ptr<IsoUtilityState> state = getIsoUtilityState(dist, matchThreshold, store, utilityFunctor, utilityMatcher);
        return -helpers::log2((FloatType)state->count/isoUtilityStatesCounts.totalcount);
    }


    FloatType localIsoPredictionComplexity(DataType x, FloatType matchThreshold = 0, bool store = true) {
        // No info, compute the utility for that new point
        if (!utilityFunctorCallback) return std::numeric_limits<FloatType>::quiet_NaN(); // no utility function given, cannot compute state
        return localIsoPredictionComplexity(x, matchThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    FloatType localIsoPredictionComplexity(DataType x, FloatType matchThreshold, bool store, UtilityFunctor u) {
        if (!predictionSetMatcherCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localIsoPredictionComplexity(x, matchThreshold, store, u, predictionSetMatcherCallback);
    }

    template <class UtilityFunctor, class PredictionMatcher>
    FloatType localIsoPredictionComplexity(DataType x, FloatType matchThreshold, bool store, UtilityFunctor utilityFunctor, PredictionMatcher predictionMatcher) {
        Distribution dist = distributionManager.distribution(x);
        return localIsoPredictionComplexity(dist, matchThreshold, store, utilityFunctor, predictionSetMatcherCallback);
    }
    
    FloatType localIsoPredictionComplexity(Distribution dist, FloatType matchThreshold = 0, bool store = true) {
        if (!utilityFunctorCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localIsoPredictionComplexity(dist, matchThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    FloatType localIsoPredictionComplexity(Distribution dist, FloatType matchThreshold, bool store, UtilityFunctor u) {
        if (!predictionSetMatcherCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localIsoPredictionComplexity(dist, matchThreshold, store, u, predictionSetMatcherCallback);
    }
    
    
    template <class UtilityFunctor, class PredictionMatcher>
    FloatType localIsoPredictionComplexity(Distribution dist, FloatType matchThreshold, bool store, UtilityFunctor utilityFunctor, PredictionMatcher predictionMatcher) {
        boost::shared_ptr<IsoPredictionState> state = getIsoPredictionState(dist, matchThreshold, store, utilityFunctor, predictionMatcher);
        return -helpers::log2((FloatType)state->count/isoPredictionStatesCounts.totalcount);
    }


    FloatType localDecisionalComplexity(DataType x, FloatType matchUtilityThreshold = 0, FloatType matchPredictionThreshold = 0, bool store = true) {
        if (!utilityFunctorCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localDecisionalComplexity(x, matchUtilityThreshold, matchPredictionThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    FloatType localDecisionalComplexity(DataType x, FloatType matchUtilityThreshold, FloatType matchPredictionThreshold, bool store, UtilityFunctor u) {
        if (!utilityMatcherCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        if (!predictionSetMatcherCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localDecisionalComplexity(x, matchUtilityThreshold, matchPredictionThreshold, store, u, utilityMatcherCallback, predictionSetMatcherCallback);
    }

    template <class UtilityFunctor, class UtilityMatcher, class PredictionMatcher>
    FloatType localDecisionalComplexity(DataType x, FloatType matchUtilityThreshold, FloatType matchPredictionThreshold, bool store, UtilityFunctor u, UtilityMatcher utilityMatcher, PredictionMatcher predictionMatcher) {
        Distribution dist = distributionManager.distribution(x);
        return localDecisionalComplexity(dist, matchUtilityThreshold, matchPredictionThreshold, store, u, utilityMatcher, predictionMatcher);
    }

    FloatType localDecisionalComplexity(Distribution dist, FloatType matchUtilityThreshold = 0, FloatType matchPredictionThreshold = 0, bool store = true) {
        if (!utilityFunctorCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localDecisionalComplexity(dist, matchUtilityThreshold, matchPredictionThreshold, store, utilityFunctorCallback);
    }

    template <class UtilityFunctor>
    FloatType localDecisionalComplexity(Distribution dist, FloatType matchUtilityThreshold, FloatType matchPredictionThreshold, bool store, UtilityFunctor u) {
        if (!utilityMatcherCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        if (!predictionSetMatcherCallback) return std::numeric_limits<FloatType>::quiet_NaN();
        return localDecisionalComplexity(dist, matchUtilityThreshold, matchPredictionThreshold, store, u, utilityMatcherCallback, predictionSetMatcherCallback);
    }

    template <class UtilityFunctor, class UtilityMatcher, class PredictionMatcher>
    FloatType localDecisionalComplexity(Distribution dist, FloatType matchUtilityThreshold, FloatType matchPredictionThreshold, bool store, UtilityFunctor u, UtilityMatcher utilityMatcher, PredictionMatcher predictionMatcher) {
        boost::shared_ptr<DecisionalState> state = getDecisionalState(dist, matchUtilityThreshold, matchPredictionThreshold, store, u, utilityMatcher, predictionMatcher);
        return -helpers::log2((FloatType)state->count/decisionalStatesCounts.totalcount);
    }



};


}

#endif
