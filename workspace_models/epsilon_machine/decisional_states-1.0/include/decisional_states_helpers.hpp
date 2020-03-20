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

#ifndef DECISIONAL_STATES_HELPERS_H
#define DECISIONAL_STATES_HELPERS_H

#include "decisional_states_NearTree_neighbour_finder.hpp"
#include "decisional_states_helpers_sized_array.hpp"

#include <boost/unordered_map.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_class.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_floating_point.hpp>

namespace decisional_states {
namespace helpers {


DecisionalStates_Detail_checker(data)
DecisionalStates_Detail_checker(prediction)

// Generic template derefs the iterator with first/second
template <class DataSet, class TypeTraits, int derefMethod>
struct DataSetIteratorUnrefSelector {
    inline static typename TypeTraits::DataType& getData(DataSet&, typename DataSet::iterator it) {
        return it->first;
    };
    inline static typename TypeTraits::PredictionType& getPrediction(DataSet&, typename DataSet::iterator it) {
        return it->second;
    };
};

// Use data/prediction from index when available
template <class DataSet, class TypeTraits>
struct DataSetIteratorUnrefSelector<DataSet,TypeTraits,1> {
    inline static typename TypeTraits::DataType getData(DataSet& d, typename DataSet::iterator it) {
        return d.data(it);
    };
    inline static typename TypeTraits::PredictionType getPrediction(DataSet& d, typename DataSet::iterator it) {
        return d.prediction(it);
    };
};

template <class DataSet, class TypeTraits = DataSet>
struct DataSetIteratorUnrefHelper {
    typedef typename TypeTraits::DataType DataType;
    typedef typename TypeTraits::PredictionType PredictionType;
    typedef typename DataSet::iterator iterator;
    typedef DataSetIteratorUnrefSelector<DataSet, TypeTraits,
        // both must be present to get the specialization
        (Check_data<DataSet, DataType& (DataSet::*)(iterator)>::found
        && Check_prediction<DataSet, PredictionType& (DataSet::*)(iterator)>::found)
        || (Check_data<DataSet, DataType (DataSet::*)(iterator)>::found
        && Check_prediction<DataSet, PredictionType (DataSet::*)(iterator)>::found)
    > UnrefHelper;
    inline static DataType getData(DataSet& d, iterator it) {
        return UnrefHelper::getData(d,it);
    };
    inline static PredictionType getPrediction(DataSet& d, iterator it) {
        return UnrefHelper::getPrediction(d,it);
    };
};


#define DecisionalStates_Detail_iterator_holder(name,Name,Parent) \
DecisionalStates_Detail_checker(name) \
DecisionalStates_Detail_type_wrapper(Name, typename T::value_type) \
 \
template <class Parent, int derefMethod> \
struct Parent ## IteratorHolderBase { \
    typename Parent::iterator it; \
    Parent ## IteratorHolderBase(Parent&, typename Parent::iterator _it) : it(_it) {} \
    inline typename Name ## Wrapper<Parent>::Name& deref() { \
        return *it; \
    } \
    inline typename Name ## Wrapper<Parent>::Name const& deref() const { \
        return *it; \
    } \
}; \
 \
template <class Parent> \
struct Parent ## IteratorHolderBase<Parent, 1> { \
    Parent& p; \
    typename Parent::iterator it; \
    Parent ## IteratorHolderBase(Parent& _p, typename Parent::iterator _it) : p(_p), it(_it) {} \
    inline typename Name ## Wrapper<Parent>::Name& deref() { \
        return p.name(it); \
    } \
}; \
 \
template <class Parent> \
struct Parent ## IteratorHolder : public Parent ## IteratorHolderBase< \
    Parent, \
    Check_ ## name<Parent, typename Name ## Wrapper<Parent>::Name (Parent::*)(typename Parent::iterator)>::found  \
> \
{ \
    Parent ## IteratorHolder(Parent& _p, typename Parent::iterator _it) \
    : Parent ## IteratorHolderBase<Parent, Check_ ## name<Parent, typename Name ## Wrapper<Parent>::Name (Parent::*)(typename Parent::iterator)>::found>(_p,_it) \
    {} \
}; \

/*
template <class Parent, int derefMethod>
struct Parent ## IteratorUnrefSelector {
    inline static typename Name ## Wrapper<Parent>::Name get(Parent&, typename Parent::iterator it) {
        return *it;
    };
};
template <class Parent>
struct Parent ## IteratorUnrefSelector<Parent,1> {
    inline static typename Name ## Wrapper<Parent>::Name get(Parent& p, typename Parent::iterator it) {
        return p.name(it);
    };
    inline static const typename Name ## Wrapper<Parent>::Name get(const Parent& p, typename Parent::iterator it) {
        return p.name(it);
    };
};

template <class Parent>
struct Parent ## IteratorUnrefHelper {
    typedef Parent ## IteratorUnrefSelector<Parent,
        Check_ ## name<Parent, typename Name ## Wrapper<Parent>::Name (Parent::*)(typename Parent::iterator)>::found
    > UnrefHelper;
    inline static DataType getData(DataSet& d, iterator it) {
        return UnrefHelper::getData(d,it);
    };
    inline static PredictionType getPrediction(DataSet& d, iterator it) {
        return UnrefHelper::getPrediction(d,it);
    };
};
*/



////////////////////

// Handle clustering algorithm iterator to get cluster
DecisionalStates_Detail_iterator_holder(cluster,Cluster,ClusteringAlgorithm)
// Handle distribution manager
DecisionalStates_Detail_iterator_holder(distribution,Distribution,DistributionManager)
// Handle cluster iterator to get elements
DecisionalStates_Detail_iterator_holder(element,Element,Cluster)


// Refine the hack for cluster data aggregation
DecisionalStates_Detail_checker(aggregate)
    
template<typename Data, int isarith = boost::is_arithmetic<Data>::value >
struct AggregateSelectorArithmeticHelper {
    inline static Data nullData() { return 0; }
};

template<typename Data>
struct AggregateSelectorArithmeticHelper<Data, 0> {
    inline static Data nullData() { return Data(); }
};

// Default: aggregate=average using operator+= and /= for arithmetic types
template <class ClusteringAlgorithm, class Container, int aggregateProvider>
struct AggregateSelector {
    typedef typename ClusterWrapper<ClusteringAlgorithm>::Cluster Cluster;
    typedef typename Container::value_type Data;

    inline static Data aggregate(ClusteringAlgorithm& algo, Cluster& cluster, Container& container) {
        Data ret = AggregateSelectorArithmeticHelper<Data>::nullData();
        int count = 0;
        for (typename Cluster::iterator it = cluster.begin(); it != cluster.end(); ++it) {
            ret += container[ClusterIteratorHolder<Cluster>(cluster,it).deref()];
            ++count;
        }
        ret /= count;
        return ret;
    }
};

// Specialization when algo.aggregate(cluster) exists
template <class ClusteringAlgorithm, class Container>
struct AggregateSelector<ClusteringAlgorithm,Container,1> {
    typedef typename ClusterWrapper<ClusteringAlgorithm>::Cluster Cluster;
    typedef typename Container::value_type Data;
    inline static Data aggregate(ClusteringAlgorithm& algo, Cluster& cluster, Container& container) {
        return algo.aggregate(cluster,container);
    }
};

// Specialization when data.aggregate(other_data) exists
template <class ClusteringAlgorithm, class Container>
struct AggregateSelector<ClusteringAlgorithm,Container,2> {
    typedef typename ClusterWrapper<ClusteringAlgorithm>::Cluster Cluster;
    typedef typename Container::value_type Data;
    inline static Data aggregate(ClusteringAlgorithm& algo, Cluster& cluster, Container& container) {
        // start with default-constructed data:
        // - it is returned when the cluster is empty, so the user can check the default
        // - this avoids problems when data contains shared ptrs
        Data ret;
        for (typename Cluster::iterator it = cluster.begin(); it != cluster.end(); ++it) ret.aggregate(container[ClusterIteratorHolder<Cluster>(cluster,it).deref()]);
        return ret;
    }
};

template <class Data, class Enable = void>
struct Check_aggregate_Data_proxy {
    enum {found = 0};
};

template <class Data>
struct Check_aggregate_Data_proxy<Data, typename boost::enable_if<boost::is_class<Data> >::type> {
    enum {found = Check_aggregate<Data, void (Data::*)(Data)>::found
                 |Check_aggregate<Data, void (Data::*)(Data&)>::found
                 |Check_aggregate<Data, void (Data::*)(const Data&)>::found
         };
};

template <class ClusteringAlgorithm, class Container>
struct AggregatorHelper {
    typedef typename ClusterWrapper<ClusteringAlgorithm>::Cluster Cluster;
    typedef typename Container::value_type Data;
    inline static Data aggregate(ClusteringAlgorithm& algo, Cluster& cluster, Container& container) {
        // Should we allow arg as Value, Value& and const Value& and Value* and const Value* ?
        return AggregateSelector<ClusteringAlgorithm,Container,
            // a ? 1 : (b ? 2 : (c ? 3 : ...))
            (Check_aggregate<ClusteringAlgorithm, Data (ClusteringAlgorithm::*)(Cluster, Container)>::found
            |Check_aggregate<ClusteringAlgorithm, Data (ClusteringAlgorithm::*)(Cluster, Container&)>::found
            |Check_aggregate<ClusteringAlgorithm, Data (ClusteringAlgorithm::*)(Cluster, const Container&)>::found
            |Check_aggregate<ClusteringAlgorithm, Data (ClusteringAlgorithm::*)(Cluster&, Container)>::found
            |Check_aggregate<ClusteringAlgorithm, Data (ClusteringAlgorithm::*)(Cluster&, Container&)>::found
            |Check_aggregate<ClusteringAlgorithm, Data (ClusteringAlgorithm::*)(Cluster&, const Container&)>::found
            |Check_aggregate<ClusteringAlgorithm, Data (ClusteringAlgorithm::*)(const Cluster&, Container)>::found
            |Check_aggregate<ClusteringAlgorithm, Data (ClusteringAlgorithm::*)(const Cluster&, Container&)>::found
            |Check_aggregate<ClusteringAlgorithm, Data (ClusteringAlgorithm::*)(const Cluster&, const Container&)>::found
            ) ? 1 : (Check_aggregate_Data_proxy<Data>::found ? 2 : 0)
        >::aggregate(algo,cluster,container);
    }
};

template <class Data, int CanAggregate = Check_aggregate_Data_proxy<Data>::found>
struct DataAggregator {
    Data sum;
    Data avg;
    int count;
    // default-initialised member numeric data types have value 0 by C++ standard
    DataAggregator() : sum(), avg(), count(0) {}
    DataAggregator(const Data& data) : sum(data), avg(data), count(1) {}
    inline void aggregate(Data& data) {
        sum += data;
        avg = sum;
        avg /= ++count;
    }
    inline operator Data& () {return avg;}
};

template <class Data>
struct DataAggregator<Data,1> {
    Data aggregated;
    DataAggregator(const Data& data) : aggregated(data) {}
    inline void aggregate(Data& data) {
        aggregated.aggregate(data);
    }
    inline operator Data& () {return aggregated;}
};

// Detection for expectation operation on distributions
// Case 1: distribution.expectation(functor)
// Case 2: distrib_manager.expectation(dist,functor)
DecisionalStates_Detail_checker(expectation)

// TODO: generic is a "no found" case with explicit error message

template <class FloatType, class DistributionManager, int method>
struct ExpectationSelector {
    typedef typename DistributionWrapper<DistributionManager>::Distribution Distribution;
    template<class Functor>
    inline static FloatType expectation(DistributionManager&, Distribution& d, Functor f) {
        return d.expectation(f);
    }
};

template <class FloatType, class DistributionManager>
struct ExpectationSelector<FloatType, DistributionManager, 1> {
    typedef typename DistributionWrapper<DistributionManager>::Distribution Distribution;
    template<class Functor>
    inline static FloatType expectation(DistributionManager& dm, Distribution& d, Functor f) {
        return dm.expectation(d,f);
    }
};

template <class FloatType, class DistributionManager>
struct ExpectationHelper {
    typedef typename DistributionWrapper<DistributionManager>::Distribution Distribution;
    template<class Functor>
    inline static FloatType expectation(DistributionManager& dm, Distribution& d, Functor f) {
        return ExpectationSelector<FloatType, DistributionManager,
            Check_expectation<DistributionManager, FloatType (DistributionManager::*)(Distribution,Functor)>::found
            |Check_expectation<DistributionManager, FloatType (DistributionManager::*)(Distribution&,Functor)>::found
            |Check_expectation<DistributionManager, FloatType (DistributionManager::*)(const Distribution&,Functor)>::found
        >::expectation(dm,d,f);
    }
};

/*
Helper that can be used to store information needed by the analyser when
the distribution manager does not want to care.
It is recommended instead to store the information as part of the distribution
directly in the distribution manager
    void* getInfo(Distribution& dist) {
        return dist->info;
    }
    void setInfo(Distribution& dist, void* info) {
        dist->info = info;
    }
Alternatively, if the operation is unsupported (the get and set are not defined),
then the hash table is used without error as a fallback
*/

template <class Distribution>
struct MapStateInfoManager {
    typedef boost::unordered_map<Distribution,void*> StateDataMap;
    StateDataMap stateDataMap;
    void* getInfo(Distribution& dist) {
        typename StateDataMap::iterator it = stateDataMap.find(dist);
        if (it == stateDataMap.end()) return 0;
        return it->second;
    }
    void setInfo(Distribution& dist, void* info) {
        stateDataMap.insert(std::pair<Distribution,void*>(dist,info));
    }

};

DecisionalStates_Detail_checker(getInfo)
DecisionalStates_Detail_checker(setInfo)

// Generic template = unsuported case = use the hash table
template <class DistributionManager, int manager>
struct StateInfoSelector {
    typedef typename DistributionWrapper<DistributionManager>::Distribution Distribution;

    // Precondition: Distributions belong to one, and only one, manager
    // We can then put all distributions in the mapper irrespectively of manager
    static MapStateInfoManager<Distribution> mapper;

    inline static void* getInfo(DistributionManager&, Distribution& dist) {
        return mapper.getInfo(dist);
    };
    inline static void setInfo(DistributionManager&, Distribution& dist, void* info) {
        mapper.setInfo(dist,info);
    };
};

// Use the getter/setter when available
template <class DistributionManager>
struct StateInfoSelector<DistributionManager,1> {
    typedef typename DistributionWrapper<DistributionManager>::Distribution Distribution;

    inline static void* getInfo(DistributionManager& manager, Distribution& dist) {
        return manager.getInfo(dist);
    };
    inline static void setInfo(DistributionManager& manager, Distribution& dist, void* info) {
        manager.setInfo(dist,info);
    };
};

template <class DistributionManager>
struct StateInfoHelper {
    typedef typename DistributionWrapper<DistributionManager>::Distribution Distribution;
    typedef StateInfoSelector<DistributionManager,
        // both must be present to get the specialization
        (Check_getInfo<DistributionManager, void* (DistributionManager::*)(Distribution)>::found
        |Check_getInfo<DistributionManager, void* (DistributionManager::*)(Distribution&)>::found
        |Check_getInfo<DistributionManager, void* (DistributionManager::*)(const Distribution&)>::found
        ) &&
        (Check_setInfo<DistributionManager, void (DistributionManager::*)(Distribution, void* info)>::found
        |Check_setInfo<DistributionManager, void (DistributionManager::*)(Distribution&, void* info)>::found
        |Check_setInfo<DistributionManager, void (DistributionManager::*)(const Distribution&, void* info)>::found
        )
    > Selector;
    inline static void* getInfo(DistributionManager& manager, Distribution& dist) {
        return Selector::getInfo(manager,dist);
    };
    inline static void setInfo(DistributionManager& manager, Distribution& dist, void* info) {
        Selector::setInfo(manager,dist,info);
    };
};


// refine the wrapper to use std::equal_to with a given type
template <class BaseClass>
struct TypeCheck_DistributionMatcher {
    template <class Base> static long sfinaeOverload(typename Base::DistributionMatcher *);
    template <class Base> static char sfinaeOverload(...);
    enum { found = sizeof(sfinaeOverload<BaseClass>(0)) == sizeof(long) };
};
template<class T, class Dist, int x> struct TypeSelector_DistributionMatcher {
    typedef typename std::equal_to<Dist> DistributionMatcher;
};
template<class T, class Dist> struct TypeSelector_DistributionMatcher <T,Dist,1> {
    typedef typename T::DistributionMatcher DistributionMatcher;
};
template<class T, class Dist>
struct DistributionMatcherWrapper {
    typedef typename TypeSelector_DistributionMatcher<T, Dist, TypeCheck_DistributionMatcher<T>::found>::DistributionMatcher DistributionMatcher;
};


// Fallback to double when DataSet does not define FloatType
DecisionalStates_Detail_type_wrapper(FloatType,double)
// Fallback to the first and second type of std::pair when no data or prediction type is provided
DecisionalStates_Detail_type_wrapper(DataType, typename T::value_type::first_type)
DecisionalStates_Detail_type_wrapper(PredictionType, typename T::value_type::second_type)
// Fallback to near tree neighbours search
#define DecisionalStates_Detail_NTNFD NearTreeNeighbourhoodFinder<typename DataTypeWrapper<T>::DataType,typename FloatTypeWrapper<T>::FloatType>
#define DecisionalStates_Detail_NTNFP NearTreeNeighbourhoodFinder<typename PredictionTypeWrapper<T>::PredictionType,typename FloatTypeWrapper<T>::FloatType>
DecisionalStates_Detail_type_wrapper(DataNeighbourhoodFinder, DecisionalStates_Detail_NTNFD)
DecisionalStates_Detail_type_wrapper(PredictionNeighbourhoodFinder, DecisionalStates_Detail_NTNFP)



DecisionalStates_Detail_checker(isMissing)


// Helpers for missing data / NaN values
// delegate to DataSet::ismissing(data/prediction) if it exists
// if not, call isnan/isinf... on the data if it is a floating point type
//         integer and non-numeric data are considered valid in the absence
//         of user-specified ismissing
// Generic template = missing function
template <class DataSet, class DataOrPred, int manager>
struct IsMissingSelector {

    template<class T>
    inline static typename boost::enable_if<boost::is_floating_point<T>, bool>::type isMissingSFINAE(T f) {
//        std::cout << "float ismissing called" << std::endl;
        return !isfinite(f);
    }
    // In the event no global function is provided and this function is used,
    // the compiler shall optimize out dead code branches of the style "if (ismissing(x)) {}"
    template<class T>
    inline static typename boost::disable_if<boost::is_floating_point<T>, bool>::type isMissingSFINAE(T) {
//        std::cout << "dummy ismissing called" << std::endl;
        return false;
    }
    inline static bool isMissing(DataSet&, const DataOrPred& dop) {
        return isMissingSFINAE<DataOrPred>(dop);
    };
    inline static bool isMissing(const DataOrPred& dop) {
        return isMissingSFINAE<DataOrPred>(dop);
    };
};

// Template specialisation: the user-provided isMissing function exists only as member
template <class DataSet, class DataOrPred>
struct IsMissingSelector<DataSet,DataOrPred,1> {
    inline static bool isMissing(DataSet& ds, DataOrPred& dop) {
        return ds.isMissing(dop);
    };
    inline static bool isMissing(DataOrPred& dop) {
        return IsMissingSelector<DataSet,DataOrPred,0>::isMissing(dop);
    };
    inline static bool isMissing(DataSet& ds, const DataOrPred& dop) {
        return ds.isMissing(dop);
    };
    inline static bool isMissing(const DataOrPred& dop) {
        return IsMissingSelector<DataSet,DataOrPred,0>::isMissing(dop);
    };
};

// Template specialisation: the user-provided isMissing exists only as a static function
template <class DataSet, class DataOrPred>
struct IsMissingSelector<DataSet,DataOrPred,2> {
    inline static bool isMissing(DataSet&, DataOrPred& dop) {
        // then ignore the argument and call the function
        return DataSet::isMissing(dop);
    };
    inline static bool isMissing(DataOrPred& dop) {
        return DataSet::isMissing(dop);
    };
    inline static bool isMissing(DataSet&, const DataOrPred& dop) {
        // then ignore the argument and call the function
        return DataSet::isMissing(dop);
    };
    inline static bool isMissing(const DataOrPred& dop) {
        return DataSet::isMissing(dop);
    };
};
// both exist
template <class DataSet, class DataOrPred>
struct IsMissingSelector<DataSet,DataOrPred,3> {
    inline static bool isMissing(DataSet& ds, DataOrPred& dop) {
        return ds.isMissing(dop);
    };
    inline static bool isMissing(DataOrPred& dop) {
        return DataSet::isMissing(dop);
    };
    inline static bool isMissing(DataSet& ds, const DataOrPred& dop) {
        return ds.isMissing(dop);
    };
    inline static bool isMissing(const DataOrPred& dop) {
        return DataSet::isMissing(dop);
    };
};


template <class DataSet, class DataOrPred>
struct IsMissingHelper {

    typedef IsMissingSelector<DataSet,DataOrPred,
        // isMissing must be present one way or another to get the specialization
        (Check_isMissing<DataSet, bool (DataSet::*)(DataOrPred)>::found
        |Check_isMissing<DataSet, bool (DataSet::*)(DataOrPred&)>::found
        |Check_isMissing<DataSet, bool (DataSet::*)(const DataOrPred&)>::found
        ) | ((Check_isMissing<DataSet, bool (*)(DataOrPred)>::found
        |Check_isMissing<DataSet, bool (*)(DataOrPred&)>::found
        |Check_isMissing<DataSet, bool (*)(const DataOrPred&)>::found
        )*2)
    > Selector;

    inline static bool isMissing(DataSet& ds, DataOrPred& dop) {
        return Selector::isMissing(ds,dop);
    }
    inline static bool isMissing(DataOrPred& dop) {
        return Selector::isMissing(dop);
    };
    inline static bool isMissing(DataSet& ds, const DataOrPred& dop) {
        return Selector::isMissing(ds,dop);
    }
    inline static bool isMissing(const DataOrPred& dop) {
        return Selector::isMissing(dop);
    };
};


// Handle symbol constraints for causal states
DecisionalStates_Detail_type_wrapper(SymbolType,void)

template<class DataSet, typename TypesTraits, typename _SymbolType>
struct DefaultTransitionFeeder {
    typedef typename DataTypeWrapper<DataSet>::DataType DataType;
    typedef _SymbolType SymbolType;
    
    // specific iterator computes the symbol & data values only once
    // and skips missing transitions in the data set
    struct iterator {
        DataSet* dataset;
        DataType before, after;
        SymbolType symbol;
        typename DataSet::iterator datait;
        typename DataSet::iterator nextit;
        
        // makes an iterator at begin()
        iterator(DataSet* _dataset) : dataset(_dataset) {
            datait = dataset->begin(); // in case there is no ++ below
            nextit = dataset->begin();
            if (nextit != dataset->end()) {
                after = DataSetIteratorUnrefHelper<DataSet,TypesTraits>::getData(*dataset,nextit);
                ++(*this); // operator++ does the job
            }
        }
        // makes an end iterator
        // See the == operator
        iterator(DataSet* _dataset, int dummy) : dataset(_dataset) {
            datait = dataset->end();
            nextit = dataset->end();
        }

        // comparison: need to compare only dataset iterators
        // compare the nextit so we do not have to decrement dataset->end()
        bool operator==(const iterator& it) {
            return nextit == it.nextit;
        }
        bool operator!=(const iterator& it) {
            return nextit != it.nextit;
        }
        
        // increment: until there exists a symbol, skip erroneous transitions
        // like when a data set has several time series, there are no transitions
        // between the end of one series and the beginning of the next
        iterator& operator++() { // prefix ++
            while (nextit!=dataset->end()) {
                datait = nextit;
                before = after;
                ++nextit; if (nextit==dataset->end()) break;
                after = DataSetIteratorUnrefHelper<DataSet,TypesTraits>::getData(*dataset,nextit);
                if (dataset->getSymbol(before,after,symbol)) break; // stop when a symbol is found
            }
            return (*this);
        }
        
        iterator operator++(int) { // postfix ++
            iterator ret = *this;
            ++(*this);
            return ret;
        }
    };
    
    DataSet& dataset;
    iterator the_beginning;
    iterator the_end;
    
    DefaultTransitionFeeder(DataSet& _dataset) : dataset(_dataset), the_beginning(&_dataset), the_end(&_dataset,1) {
    }
    
    inline iterator begin() {return the_beginning;}
    inline iterator end() {return the_end;}
    
    // by construction the iterator skips missing symbols
    // so unlike the data set version we return directly the symbol here
    // as we know it is not missing for a valid iterator
    // This also uniformises the feeder API
    SymbolType& getSymbol(iterator& it) {
        return it.symbol;
    }
    
    DataType& getDataBeforeTransition(iterator& it) {
        return it.before;
    }
    
    DataType& getDataAfterTransition(iterator& it) {
        return it.after;
    }
    
    bool isMissing(DataType& x) {
        return helpers::IsMissingHelper<DataSet,DataType>::isMissing(dataset,x);
    }
};

// No symbol version
template<class DataSet, typename TypesTraits>
struct DefaultTransitionFeeder<DataSet, TypesTraits, void> {

    typedef typename DataTypeWrapper<DataSet>::DataType DataType;
    
    // specific iterator computes the data values only once
    // we do not have symbols to tell whether the transition is missing
    // The user may always provide her/his own transition feeder
    // for multi-series file for ex.
    struct iterator {
        DataSet* dataset;
        DataType before, after;
        typename DataSet::iterator datait;
        typename DataSet::iterator nextit;
        
        // makes an iterator at begin()
        iterator(DataSet* _dataset) : dataset(_dataset) {
            datait = dataset->begin(); // in case there is no ++ below
            nextit = dataset->begin();
            if (nextit != dataset->end()) {
                after = DataSetIteratorUnrefHelper<DataSet,TypesTraits>::getData(*dataset,nextit);
                ++(*this); // operator++ does the job
            }
        }
        // makes an end iterator
        // See the == operator
        iterator(DataSet* _dataset, int dummy) : dataset(_dataset) {
            datait = dataset->end();
            nextit = dataset->end();
        }

        // comparison: need to compare only dataset iterators
        // compare the nextit so we do not have to decrement dataset->end()
        bool operator==(const iterator& it) {
            return nextit == it.nextit;
        }
        bool operator!=(const iterator& it) {
            return nextit != it.nextit;
        }
        
        // increment: unlike the symbol version, take all transitions
        iterator& operator++() { // prefix ++
            if (nextit!=dataset->end()) {
                datait = nextit;
                before = after;
                ++nextit; if (nextit==dataset->end()) return (*this);
                after = DataSetIteratorUnrefHelper<DataSet,TypesTraits>::getData(*dataset,nextit);
            }
            return (*this);
        }
        
        iterator operator++(int) { // postfix ++
            iterator ret = *this;
            ++(*this);
            return ret;
        }
    };
    
    DataSet& dataset;
    iterator the_beginning;
    iterator the_end;
    
    DefaultTransitionFeeder(DataSet& _dataset) : dataset(_dataset), the_beginning(&_dataset), the_end(&_dataset,1) {
    }
    
    inline iterator begin() {return the_beginning;}
    inline iterator end() {return the_end;}
    
    DataType& getDataBeforeTransition(iterator& it) {
        return it.before;
    }
    
    DataType& getDataAfterTransition(iterator& it) {
        return it.after;
    }

    
    bool isMissing(DataType& x) {
        return helpers::IsMissingHelper<DataSet,DataType>::isMissing(dataset,x);
    }
};


template<class DataSet>
struct DefaultTypeTraitsDS {
    typedef typename FloatTypeWrapper<DataSet>::FloatType FloatType;
    typedef typename DataTypeWrapper<DataSet>::DataType DataType;
    typedef typename PredictionTypeWrapper<DataSet>::PredictionType PredictionType;
    typedef typename DataNeighbourhoodFinderWrapper<DataSet>::DataNeighbourhoodFinder DataNeighbourhoodFinder;
    typedef typename PredictionNeighbourhoodFinderWrapper<DataSet>::PredictionNeighbourhoodFinder PredictionNeighbourhoodFinder;
    typedef DefaultTransitionFeeder<DataSet, DefaultTypeTraitsDS<DataSet>, typename SymbolTypeWrapper<DataSet>::SymbolType> TransitionFeeder;
};

template<class DataSet, class Optimiser>
struct DefaultTypeTraits {
    // Fallback to double when DataSet does not define FloatType
    typedef typename FloatTypeWrapper<DataSet>::FloatType FloatType;
    typedef typename DataTypeWrapper<DataSet>::DataType DataType;
    typedef typename PredictionTypeWrapper<DataSet>::PredictionType PredictionType;
    typedef typename Optimiser::PredictionSet PredictionSet;
    typedef typename DataNeighbourhoodFinderWrapper<DataSet>::DataNeighbourhoodFinder DataNeighbourhoodFinder;
    typedef typename PredictionNeighbourhoodFinderWrapper<DataSet>::PredictionNeighbourhoodFinder PredictionNeighbourhoodFinder;
    typedef DefaultTransitionFeeder<DataSet, DefaultTypeTraitsDS<DataSet>, typename SymbolTypeWrapper<DataSet>::SymbolType> TransitionFeeder;
};

DecisionalStates_Detail_type_wrapper(SampleType, typename T::value_type)






}

}


#endif


