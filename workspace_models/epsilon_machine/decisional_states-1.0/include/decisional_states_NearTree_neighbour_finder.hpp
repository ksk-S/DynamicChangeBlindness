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

#ifndef DECISIONAL_STATES_NEARTREE_NEIGHBOUR_FINDER_H
#define DECISIONAL_STATES_NEARTREE_NEIGHBOUR_FINDER_H

#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "decisional_states_TNear.hpp"

#include "decisional_states_helpers_math.hpp"
#include "decisional_states_helpers_container.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>

namespace decisional_states {


namespace helpers {

template <class C, typename NumericType>
struct DistanceToL2 {
    template <class BaseClass, typename FunctionType>
    struct Check_member_distance_to {
        template <FunctionType> struct Finder;
        template <class Base> static long sfinaeOverload(Finder<&Base::distance_to> *);
        template <class Base> static char sfinaeOverload(...);
        enum { found = sizeof(sfinaeOverload<BaseClass>(0)) == sizeof(long) };
    };

    template <class BaseClass, typename FunctionType>
    struct Check_member_distance2_to {
        template <FunctionType> struct Finder;
        template <class Base> static long sfinaeOverload(Finder<&Base::distance2_to> *);
        template <class Base> static char sfinaeOverload(...);
        enum { found = sizeof(sfinaeOverload<BaseClass>(0)) == sizeof(long) };
    };

    template <typename X>
    struct Check_supports_member_func {
        template <class Base> static long sfinaeOverload(void(Base::*)(void));
        template <class Base> static char sfinaeOverload(...);
        enum { found = sizeof(sfinaeOverload<X>(0)) == sizeof(long) };
    };

    template <typename X, int is_class = 0>
    struct Selector_supports_member_func {
        enum { found = 0 };
    };
    template <typename X>
    struct Selector_supports_member_func<X,1> {
        enum { found = Check_member_distance_to<X, NumericType (X::*)(const X&) const>::found };
    };

    template <typename X, int is_class = 0>
    struct Selector_supports_member_func2 {
        enum { found = 0 };
    };
    template <typename X>
    struct Selector_supports_member_func2<X,1> {
        enum { found = Check_member_distance2_to<X, NumericType (X::*)(const X&) const>::found };
    };

    template <class BaseClass, int dist_method> struct Selector {
        template <typename T>
        inline static typename boost::disable_if<boost::is_arithmetic<T>, NumericType>::type distance_to_selected(const T& a, const T& b) {
            NumericType ret = 0;
            for (int i=0; i<(int)a.size(); ++i) ret += (a[i] - b[i])*(a[i] - b[i]);
            return helpers::sqrt(ret);
        }
        template <typename T>
        inline static typename boost::disable_if<boost::is_arithmetic<T>, NumericType>::type distance2_to_selected(const T& a, const T& b) {
            NumericType ret = 0;
            for (int i=0; i<(int)a.size(); ++i) ret += (a[i] - b[i])*(a[i] - b[i]);
            return ret;
        }
        template <typename T>
        inline static typename boost::enable_if<boost::is_arithmetic<T>, NumericType>::type distance_to_selected(const T& a, const T& b) {
            return helpers::fabs(a-b);
        }
        template <typename T>
        inline static typename boost::enable_if<boost::is_arithmetic<T>, NumericType>::type distance2_to_selected(const T& a, const T& b) {
            return (a-b)*(a-b);
        }
        inline static NumericType distance_to(const BaseClass& a,const BaseClass& b) {
            return distance_to_selected<BaseClass>(a,b);
        }
        inline static NumericType distance2_to(const BaseClass& a,const BaseClass& b) {
            return distance2_to_selected<BaseClass>(a,b);
        }
    };
    template <class BaseClass> struct Selector<BaseClass,1> {
        inline static NumericType distance_to(const BaseClass& a,const BaseClass& b) {
            return a.distance_to(b);
        }
        inline static NumericType distance2_to(const BaseClass& a,const BaseClass& b) {
            return a.distance2_to(b);
        }
    };

    inline static NumericType distance_to(const C& a, const C& b) {
        return Selector<C, Selector_supports_member_func<C, Check_supports_member_func<C>::found>::found >::distance_to(a,b);
    }

    inline static NumericType distance2_to(const C& a, const C& b) {
        return Selector<C, Selector_supports_member_func2<C, Check_supports_member_func<C>::found>::found >::distance2_to(a,b);
    }
};

struct IsFiniteHelper {
    template <typename T>
    inline static typename boost::disable_if<boost::is_arithmetic<T>, bool>::type isfinite(const T& a) {
        int iend = a.size();
        for (int i=0; i<iend; ++i) if (!std::isfinite(a[i])) return false;
        return true;
    }
    template <typename T>
    inline static typename boost::enable_if<boost::is_arithmetic<T>, bool>::type isfinite(T a) {
        return std::isfinite(a);
    }
};

}


template<class Element, typename FloatType = double, typename CNearTree_helpers::DistTypeMaker<Element,FloatType>::Distance glob_distance_to = &helpers::DistanceToL2<Element,FloatType>::distance_to, typename CNearTree_helpers::DistTypeMaker<Element,FloatType>::Distance glob_distance2_to = &helpers::DistanceToL2<Element,FloatType>::distance2_to>
struct NearTreeNeighbourhoodFinder {

    struct TreeElement {
        const Element* e; // element in original container
        int index;  // indice in original container
        TreeElement(const Element* _e, int _i) : e(_e), index(_i) {}
        TreeElement() : e(0), index(-1) {} // ensure crash if used, need to be copied from another object first
        inline FloatType distance_to(const TreeElement& other) const {
            return glob_distance_to(*e,*other.e);
        }
    };

    boost::shared_ptr<CNearTree<TreeElement,FloatType> > tree;

    /// API: prepare the internal structures with the given container of Elements
    template<class Container>
    void setup(const Container& c) {
        tree = boost::shared_ptr<CNearTree<TreeElement,FloatType> >(new CNearTree<TreeElement,FloatType>());
        int i = 0;
        for(typename Container::const_iterator it = c.begin(); it != c.end(); ++it) {
            if (helpers::IsFiniteHelper::isfinite(*it)) tree->Insert(TreeElement(&(*it), i++));
            else ++i;
        }
    }

    /// API: find all neighbors fullfilling the given predicate
    ///      calls back the given functor on elements that were found
    /// params:
    /// - container of elements at which to launch the search
    /// - query distance for being in neighbourhood. Here euclidian, squared, to match Gaussian Kernel implementation
    /// - action to perform on all elements that were found, passing them as a container of indices
    template<class Container, class Action>
    void find(const Container& points, FloatType sqDist, Action action) {
        // for each query point
        int nquery = points.size();
        FloatType radius = helpers::sqrt(sqDist);
#if defined(DECISIONAL_STATES_NEARTREE_PARALLEL_NEIGHBOUR_FINDER_H) && defined(_OPENMP)
#pragma omp parallel for
#endif
        for (int qidx = 0; qidx < nquery; ++qidx) {
            if (!helpers::IsFiniteHelper::isfinite(points[qidx])) continue;
            std::vector<TreeElement> found;
            tree->FindInSphere(radius, found, TreeElement(&points[qidx],qidx) );
            // launch user action on all indices that were found
            // API: query point index, begin iterator, end iterator
            helpers::ContainerMemberAdaptor<
                std::vector<TreeElement>,
                int,
                &TreeElement::index
            > idxContainer(found);
            action(qidx, idxContainer.begin(), idxContainer.end());
        }
    }

    template<class Container, class Action>
    void findNearest(const Container& points, Action action) {
        // for each query point
        int nquery = points.size();
        FloatType radius = std::numeric_limits<FloatType>::max();
#if defined(DECISIONAL_STATES_NEARTREE_PARALLEL_NEIGHBOUR_FINDER_H) && defined(_OPENMP)
#pragma omp parallel for
#endif
        for (int qidx = 0; qidx < nquery; ++qidx) {
            if (!helpers::IsFiniteHelper::isfinite(points[qidx])) continue;
            TreeElement nearest;
            bool found = tree->NearestNeighbor(radius, nearest, TreeElement(&points[qidx],qidx) );
            // launch user action on the nearest neighbour, if any
            if (found) action(qidx, nearest.index);
        }
    }

    static inline FloatType distance(const Element& e1, const Element& e2) {
        return glob_distance_to(e1,e2);
    }

    static inline FloatType distance2(const Element& e1, const Element& e2) {
        return glob_distance2_to(e1,e2);
    }

};

}

#endif
