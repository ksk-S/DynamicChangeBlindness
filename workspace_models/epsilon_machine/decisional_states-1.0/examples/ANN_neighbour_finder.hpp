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

#ifndef DECISIONAL_STATES_ANN_NEIGHBOUR_FINDER_H
#define DECISIONAL_STATES_ANN_NEIGHBOUR_FINDER_H

/// Implement the neighbour finder interface with the ANN backend

#include "ANN/ANN.h"

#include <boost/shared_ptr.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>

template<class Element>
struct ANN_NeighbourhoodFinder {

    boost::shared_ptr<ANNkd_tree> kdTree;
    ANNpointArray datapoints;

    ANN_NeighbourhoodFinder() : datapoints(0) {}
    ~ANN_NeighbourhoodFinder() {if (datapoints) annDeallocPts(datapoints); datapoints = 0;}

    /// API: prepare the internal structures with the given container of Elements
    template<class Container>
    void setup(const Container& c) {
        int npoints = c.size();
        if (npoints==0) return;
        int dim = dim_of<Element>(*c.begin());
        datapoints = annAllocPts(npoints, dim);
        // container of elements
        for (int i=0; i<npoints; ++i) {
            copyElementToPoint<Element>(datapoints[i], c[i]);
        }
        // finally create the kdTree
        kdTree = boost::shared_ptr<ANNkd_tree>(new ANNkd_tree(datapoints,npoints,dim));

    }

    /// API: find all neighbors fullfilling the given predicate
    ///      calls back the given functor on elements that were found
    /// params:
    /// - container of elements at which to launch the search
    /// - query distance for being in neighbourhood. Here euclidian, squared, to match Gaussian Kernel implementation
    /// - action to perform on all elements that were found, passing them as a container of indices
    template<class Container, class Action>
    void find(const Container& points, double sqDist, Action action) {
        // for each query point
        int nquery = points.size();
        // TODO: check thread-safety, then omp parallel for
        for (int qidx = 0; qidx < nquery; ++qidx) {
            // ANN manual for fixed-radius queries recommends:
            // 1. first perform a dummy search to get the number of points in radius
            ANNcoord q[kdTree->theDim()];
            copyElementToPoint<Element>(&q[0],points[qidx]);
            int numneighbours = kdTree->annkFRSearch(q, sqDist, 0, 0, 0);
            // 2. Allocate the right number of points and do the real search
            ANNidxArray neighboursIdx = new ANNidx[numneighbours];
            kdTree->annkFRSearch(q, sqDist, numneighbours, neighboursIdx, 0);
            // launch user action on all indices that were found
            // API: query point index, begin iterator, end iterator
            action(qidx, neighboursIdx, neighboursIdx+numneighbours);
            delete[] neighboursIdx;
        }
    }


/// internal helpers
    template <typename T>
    typename boost::disable_if<boost::is_arithmetic<T>, int>::type dim_of(const T& a) {
        return a.size();
    }

    template <typename T>
    typename boost::enable_if<boost::is_arithmetic<T>, int>::type dim_of(const T& a) {
        return 1;
    }


    template <typename T>
    typename boost::disable_if<boost::is_arithmetic<T> >::type copyElementToPoint(ANNpoint p, const T& e) {
        for (int i=0; i<(int)e.size(); ++i) p[i] = e[i];
    }

    template <typename T>
    typename boost::enable_if<boost::is_arithmetic<T> >::type copyElementToPoint(ANNpoint p, const T& e) {
        p[0] = e;
    }


};

#endif
