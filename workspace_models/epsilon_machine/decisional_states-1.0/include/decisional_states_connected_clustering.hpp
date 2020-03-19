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

#ifndef DECISIONAL_STATES_CONNECTED_CLUSTERING_H
#define DECISIONAL_STATES_CONNECTED_CLUSTERING_H

#include <vector>
#include <map>
#include <set>
#include <list>

#include <boost/shared_ptr.hpp>
#include <boost/functional/hash.hpp>
//#include "decisional_states_helpers.hpp"

#ifdef _OPENMP
#include <omp.h>
#include <string.h> // memcpy
#endif

namespace decisional_states {

// API
// The cluster method is passed a matcher and a random access container of unspecified but compatible type:
// - Container::value_type is defined
// - container[int] is defined and returns a reference of value_type
// - container.size() is valid
// - match(x,y) applies to variables x and y of this value_type
// The clustering algo must cluster together the indices of the elements in the container
// This allows:
// - type independance: any random access container can be passed, the int indices are universal.
//   thus the clustering algo can be created even if the container type to cluster is not yet known
// - the elements themselves are not duplicated in the clusters, just the indices


// Tags (traits) to specialize the behavior of the clustering algorithm below
struct NoAggregation {};
struct SumAggregation {};
struct AverageAggregation {};

// Performs a simple hierarchical clustering with the single-linkage method.
// It is guaranteed that data in different clusters do not match according to
// a user-provided predicate.
// Clusters are then maximal connected components, with the property that
// there is always a path of pairwise matching data (x,y), (y,z), (z,t)...
// even if (x,t) does not match.
// Note: This cluster method performs N(N-1)/2 match comparisons
//       and needs only O(N) memory
template<class Aggregation = NoAggregation>
struct ConnectedClustering {

    struct NoConstraints {
        typedef std::pair<int,int>* iterator;
        typedef std::pair<int,int>* const_iterator;
        iterator begin() const {return 0;}
        iterator end() const {return 0;}
        template<class IndexClassProvider>
        inline void getSplitConstraints(IndexClassProvider indexClassProvider) {}
        template<class IndexClassProvider>
        inline void getMergeConstraints(IndexClassProvider indexClassProvider) {}
    };

    typedef std::vector< int > intVector;
    struct Cluster : public boost::shared_ptr<intVector> {
        Cluster() {}
        Cluster(intVector* p) : boost::shared_ptr<intVector>(p) {}
        typedef intVector::value_type value_type;
        typedef intVector::iterator iterator;
        iterator begin() {return (*(boost::shared_ptr<intVector>*)this)->begin();}
        iterator end() {return (*(boost::shared_ptr<intVector>*)this)->end();}
        typedef intVector::const_iterator const_iterator;
        const_iterator begin() const {return (*(boost::shared_ptr<intVector>*)this)->begin();}
        const_iterator end() const {return (*(boost::shared_ptr<intVector>*)this)->end();}
        size_t size() const {return (*(boost::shared_ptr<intVector>*)this)->size();}
        int& operator[](int i) {return (*(boost::shared_ptr<intVector>*)this)->operator[](i);}
    };
    typedef std::vector< Cluster > ClusterVector;
    typedef typename ClusterVector::iterator iterator;
    boost::shared_ptr< ClusterVector > allclusters;
    iterator begin() {return allclusters->begin();}
    iterator end() {return allclusters->end();}
    typedef typename ClusterVector::value_type value_type;
    size_t size() const {return allclusters->size();}
    Cluster& operator[](int i) {return allclusters->operator[](i);}

    std::vector<int> clusterIndices; // gives the current cluster index for each of the N original elements
    struct IndexClassProvider {
        std::vector<int>& clusterIndices;
        IndexClassProvider(std::vector<int>& _clusterIndices) : clusterIndices(_clusterIndices) {}
        int operator()(int idx) {
            return clusterIndices[idx];
        }
    };


    ConnectedClustering() : allclusters(new ClusterVector()) {}

    template<class MatchPredicate, class Container>
    void cluster(Container& container, MatchPredicate matchPredicate) {
        cluster(container, matchPredicate, NoConstraints());
    }

    /**
      container: contains the elements that must be matched
      matchPredicate: acts on the container elements
      constraints (optional): pairs of indices in the container that shall be put together
    */
    template<class MatchPredicate, class Container, class ConstraintsSet>
    void cluster(Container& container, MatchPredicate matchPredicate, ConstraintsSet constraints) {

        int N = helpers::ContainerSizeWrapper<Container>::size(container);

        clusterIndices.resize(N);
        for (int i=0; i<N; ++i) clusterIndices[i] = i; // to begin with, each point in its own cluster

        std::map<int, std::set<int> > constraintsById;

        // merge constraints - one point per cluster = merge constraints defined on points here
        constraints.getMergeConstraints(IndexClassProvider(clusterIndices));
        for(typename ConstraintsSet::const_iterator it = constraints.begin(); it != constraints.end(); ++it) {
            constraintsById[it->first].insert(it->second);
            constraintsById[it->second].insert(it->first);
        }

        // remaining unmatched elements
        // once an element is affected to a cluster, it is removed from that list
        // each element in the cluster is matched against unaffected ones, so as to propagate connected components
        // the same way, elements are not rematched against elements from the same cluster
        // Once all elements in a cluster are matched, no remaining element can belong to that cluster
        // In each case anyway, a newly introduced element in any cluster is always
        // matched against all remaining ones => O(N(N-1)/2)
        
        // Removing one element:
        // - need to check whether it breaks its cluster: the point could be the only link between, say, 3 subgroups.
        // - => The only way is to recluster the old cluster completely (extreme case: all points are now separated)
        // - Other clusters are not impacted at all: removing a point does not create new links
        // - Similarly, points in the old cluster cannot be affected to other clusters: algo complexity is localized
        
        std::vector<int> unmatchedIndices(N);
        for(int i=0; i<N; ++i) unmatchedIndices[i] = i; // to begin with, the order is subject to change
        int Nunmatched = N;

        while (Nunmatched > 0) {
            // a new cluster for a subset of the remaining indices
            Cluster currentCluster(new intVector());
            currentCluster->push_back(unmatchedIndices[0]);
            unmatchedIndices[0] = unmatchedIndices[--Nunmatched];
            int currentClusterElement = 0;
            // repeat until all elements in the current cluster are done, and there remains elements in the list
            while (currentClusterElement < (int)currentCluster->size() && (Nunmatched > 0)) {
                typename Container::value_type currentValue = container[currentCluster[currentClusterElement]];
                const std::set<int>* idSetToMerge = 0;
                std::map<int, std::set<int> >::const_iterator cit = constraintsById.find(currentCluster[currentClusterElement]);
                if (cit != constraintsById.end()) idSetToMerge = &(cit->second);

                // match all remaining values to the current element value
#if defined(_OPENMP)
                int nthreads = omp_get_max_threads();
#pragma omp parallel
                {
                    int tnum = omp_get_thread_num();
                    int midx = tnum * Nunmatched / nthreads;
                    int end = (tnum+1) * Nunmatched / nthreads;
                    while (midx<end) {
                        int cidx = unmatchedIndices[midx];
                        if ((idSetToMerge!=0 && idSetToMerge->find(cidx)!=idSetToMerge->end()) || matchPredicate(currentValue,container[cidx]) > 0) {
                            unmatchedIndices[midx] = -1;
#pragma omp critical
                            currentCluster->push_back(cidx);
                        }
                        ++midx;
                    }
                }
                // now that all matchs are done, restore the global vector in the single-threaded section
                while(Nunmatched>0 && unmatchedIndices[Nunmatched-1] == -1) --Nunmatched;
                for (int midx=0; midx<Nunmatched;) {
                    if (unmatchedIndices[midx] == -1) {
                        unmatchedIndices[midx] = unmatchedIndices[--Nunmatched];
                        while(Nunmatched>0 && unmatchedIndices[Nunmatched-1] == -1) --Nunmatched;
                        continue;
                    }
                    ++midx;
                }
#else
                // single-threaded version can act directly on the global vector
                for (int midx=0; midx<Nunmatched;) {
                    int cidx = unmatchedIndices[midx];
                    if ((idSetToMerge!=0 && idSetToMerge->find(cidx)!=idSetToMerge->end()) || matchPredicate(currentValue,container[cidx]) > 0) {
                        currentCluster->push_back(cidx);
                        unmatchedIndices[midx] = unmatchedIndices[--Nunmatched];
                        continue;
                    }
                    ++midx;
                }
#endif

                // if nothing was added in this pass, the loop ends.
                // Otherwise, match the newly added element to remaining ones to propagate the connected component
                ++currentClusterElement;
            }

            int clusterIndex = allclusters->size(); // before the push
            // add this cluster to the global clusters
            allclusters->push_back(currentCluster);

            // maintain cluster indices for all elements in this cluster, now that the cluster index is known
            for (typename Cluster::iterator it = currentCluster.begin(); it != currentCluster.end(); ++it) {
                clusterIndices[*it] = clusterIndex;
            }
        }

        // Now clusters are merged, but some uncheked constraints may have appeared
        // ex: causal states + symbol point to unique other state
        //     possibly individually this property hold for each x in the causal state
        //     but now that they are merged there might exist x1->y1 and x2->y2 with the same symbol,
        //     and with x1 and x2 in the same state but not y1 and y2
        //     subdivide states so x1 and x2, that were incorrectly merged, do not appear together anymore
        // Perform cycles of split/merge until all constraints are respected
        // It's up to the user not to provide contradictions, or the loop will never end
        for (int convergenceCounter = 0; convergenceCounter < 3; ++convergenceCounter) {
            bool allConstraintsAreSatisfied = true;

            // we divide each state into (x that respect constraints, x that do not) = 2 sub-states
            // iteratively subdivide "x that do not" into consistent sub-states
            while(true) {
                constraints.getSplitConstraints(IndexClassProvider(clusterIndices));
                if (constraints.begin() == constraints.end()) break;

                std::map<int, std::set<int> > to_split;
                for(typename ConstraintsSet::const_iterator it = constraints.begin(); it != constraints.end(); ++it) {
                    int idx = it->first;
                    int clusterToSplit = it->second;
                    if (clusterToSplit == clusterIndices[idx]) to_split[clusterToSplit].insert(idx);
                }
                if (to_split.empty()) break;
                allConstraintsAreSatisfied = false;
                // split sets - can now run cluster by cluster thanks to the map aggregating constraints above
                for (std::map<int,std::set<int> >::iterator it = to_split.begin(); it != to_split.end(); ++it) {
                    // at this point there is at least one element to split: the entries in to_split are created above on access/insert
                    Cluster splittedCluster(new intVector());  // to hold the "x that do not respect constraint" above
                    int splittedClusterIndex = allclusters->size();
                    Cluster& currentCluster = (*allclusters)[it->first];
                    std::set<int>& splitted = it->second;
                    // transfer elements to split from the current to a new cluster
                    for (std::set<int>::iterator sit = splitted.begin(); sit != splitted.end(); ++sit) {
                        for (int cidx = 0; cidx < (int)currentCluster.size(); ++cidx) {
                            if (currentCluster[cidx]==*sit) {
                                currentCluster[cidx] = currentCluster->back();
                                currentCluster->pop_back();
                                splittedCluster->push_back(*sit);
                                clusterIndices[*sit] = splittedClusterIndex;
                            }
                        }
                    }
                    // add the new cluster
                    allclusters->push_back(splittedCluster);
                }
                // loop again to check for more split constraints, some extracted x might still not like each other
            }

            // Are all merge constraints respected ?
            std::vector<int> clusterTree(allclusters->size()); // original cluster => merge target
            for (int i=0; i<(int)clusterTree.size(); ++i) clusterTree[i] = i;
            constraints.getMergeConstraints(IndexClassProvider(clusterIndices));
            bool treeModified = false;
            for(typename ConstraintsSet::const_iterator it = constraints.begin(); it != constraints.end(); ++it) {
                int idx = it->first;
                int clusterToMerge = it->second;
                int tidx = std::min(clusterToMerge, clusterIndices[idx]);
                int midx = std::max(clusterToMerge, clusterIndices[idx]);
                // make the root node for midx point to the root node for tidx
                while (clusterTree[tidx]!=tidx) tidx = clusterTree[tidx];
                while (clusterTree[midx]!=midx) midx = clusterTree[midx];
                if (tidx==midx) continue; // merge already taken into account
                allConstraintsAreSatisfied = false;
                treeModified = true;
                int min_root = std::min(tidx, midx);
                int root_to_update = std::max(tidx, midx);
                clusterTree[root_to_update]=min_root;
            }

            // stop if all constraints are OK
            if (allConstraintsAreSatisfied) break;

            // otherwise, merge clusters that need it
            if (treeModified) {
                // Make each element point directly to its root
                for (int i=0; i<(int)clusterTree.size(); ++i) {
                    int root = clusterTree[i];
                    while (clusterTree[root]!=root) root = clusterTree[root];
                    clusterTree[i] = root;
                }
                // roots are preserved, all other clusters are merged and destroyed
                // thanks to clusters being shared_ptr, we can work on a second vector of clusters (itself shared_ptr, fast copy below)
                boost::shared_ptr< ClusterVector > newclusters = boost::shared_ptr< ClusterVector >(new ClusterVector());
                for (int i=0; i<(int)clusterTree.size(); ++i) {
                    // root ? => preserved
                    if (clusterTree[i] == i) newclusters->push_back((*allclusters)[i]);
                    // not a root => merge to the real root
                    else {
                        (*allclusters)[clusterTree[i]]->insert((*allclusters)[clusterTree[i]]->end(), (*allclusters)[i]->begin(), (*allclusters)[i]->end());
                        (*allclusters)[i]->clear();
                    }
                }
                // all clusters are now merged. restore the indices and main cluster vector
                allclusters = newclusters;
                for (int i = 0; i<(int)allclusters->size(); ++i) {
                    for (typename Cluster::iterator it = (*allclusters)[i]->begin(); it != (*allclusters)[i]->end(); ++it) {
                        clusterIndices[*it] = i;
                    }
                }
            }

        } // end of constraints validation loop

    }

};

template<>
struct ConnectedClustering<SumAggregation> : public ConnectedClustering<NoAggregation> {
    // propagate types
    typedef ConnectedClustering<NoAggregation>::Cluster Cluster;
    typedef ConnectedClustering<NoAggregation>::iterator iterator;
    typedef ConnectedClustering<NoAggregation>::value_type value_type;

    // Add in the aggregation method
    template<class Container>
    typename Container::value_type aggregate(Cluster& cluster, Container& container) {
        typename Cluster::iterator it = cluster.begin();
        if (it==cluster.end()) return typename Container::value_type();
        typename Container::value_type ret = container[*it];
        for (++it; it != cluster.end(); ++it) ret = ret + container[*it];
        return ret;
    }
};

template<>
struct ConnectedClustering<AverageAggregation> : public ConnectedClustering<SumAggregation> {
    // propagate types
    typedef ConnectedClustering<SumAggregation>::Cluster Cluster;
    typedef ConnectedClustering<SumAggregation>::iterator iterator;
    typedef ConnectedClustering<SumAggregation>::value_type value_type;

    // Specialize the aggregation method from sum to average
    template<class Container>
    typename Container::value_type aggregate(Cluster& cluster, Container& container) {
        if (cluster->empty()) return typename Container::value_type();
        return ConnectedClustering<SumAggregation>::aggregate(cluster,container) * (1.0 / cluster->size());
    }
};

}

#endif
