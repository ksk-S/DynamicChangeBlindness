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

#ifndef DECISIONAL_STATES_EXHAUSTIVE_NEIGHBOUR_FINDER_H
#define DECISIONAL_STATES_EXHAUSTIVE_NEIGHBOUR_FINDER_H

#include "decisional_states_helpers_sized_array.hpp"

#include <vector>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>

#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace decisional_states {

template<class Element, typename FloatType = double>
struct ExhaustiveNeighbourhoodFinder {

    helpers::SizedArray<Element> elements;

    /// API: prepare the internal structures with the given container of Elements
    template<class Container>
    void setup(const Container& c) {
        elements = c;
    }

    /// API: find all neighbors fullfilling the given predicate
    ///      calls back the given functor on elements that were found
    /// params:
    /// - container of elements at which to launch the search
    /// - query distance for being in neighbourhood, euclidian, squared
    /// - action to perform on all elements that were found, passing them as a container of indices
    template<class Container, class Action>
    void find(const Container& points, FloatType sqDist, Action action) {
        // for each query point
        int nquery = points.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int qidx = 0; qidx < nquery; ++qidx) {
            const Element& queryPoint = points[qidx];
            // Brute-force exhaustive search
            std::vector<int> result;
            for (std::size_t i=0; i<elements.size(); ++i) {
                // if distSq returns NaN then result is not pushed
                if (distSq<Element>(elements[i],queryPoint)<=sqDist) {
                    result.push_back(i);
                }
            }
            // launch user action on all indices that were found
            // actions may be launched in parallel, user is warned
            action(qidx, result.begin(), result.end());
        }
    }


/// internal helpers
    template <typename T>
    typename boost::disable_if<boost::is_arithmetic<T>, FloatType>::type distSq(const T& a, const T&b) {
        std::size_t dim = a.size();
        assert(dim==b.size());
        FloatType ret = FloatType(0);
        for (std::size_t d=0; d<dim; ++d) ret += (a[d] - b[d]) * (a[d] - b[d]);
        return ret;
    }

    template <typename T>
    typename boost::enable_if<boost::is_arithmetic<T>, FloatType>::type distSq(const T& a, const T&b) {
        return (a-b)*(a-b);
    }


};


}
#endif
