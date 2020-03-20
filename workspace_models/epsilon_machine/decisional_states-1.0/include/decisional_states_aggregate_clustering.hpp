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

#ifndef DECISIONAL_STATES_AGGREGATE_CLUSTERING_H
#define DECISIONAL_STATES_AGGREGATE_CLUSTERING_H

#include <vector>
#include <map>
#include <set>
#include <list>

#include <boost/shared_ptr.hpp>
#include <boost/functional/hash.hpp>
#include <boost/random.hpp>

namespace decisional_states {

struct FirstIndexProvider {
    FirstIndexProvider() {}
    int operator()(int rangemin, int rangemax) {return rangemin;}
};

struct RandomProvider {
    boost::mt19937 rng;
    RandomProvider(unsigned int _seed = 42) {rng.seed(_seed);}
    int operator()(int rangemin, int rangemax) {
        return rangemin + (rng() % (rangemax-rangemin+1));
    }
};

// clustering algorithm that aggregates the member of each cluster into one average value
// ex: clustering probability distribution, get average distribution
// matchs are performed between clusters with the refinement that clusters may be merged
// after a succesful insert
// Each point is only visited once, but the cluster merge is O(K) with K the final number
// of clusters => global algorithm is O(K.N)

template<class RandomProvider = decisional_states::RandomProvider>
struct AggregateClustering {

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
    
    int npass;
    RandomProvider randomProvider;

    // Default is set to 3-pass, use random algorithm default constructor
    AggregateClustering() : allclusters(new ClusterVector()), npass(3) {}
    AggregateClustering(int _npass) : allclusters(new ClusterVector()), npass(_npass>1?_npass:3) {}
    AggregateClustering(RandomProvider _randomProvider) : allclusters(new ClusterVector()), npass(3), randomProvider(_randomProvider) {}
    AggregateClustering(int _npass, RandomProvider _randomProvider) : allclusters(new ClusterVector()), npass(_npass>1?_npass:3), randomProvider(_randomProvider) {}

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

        std::vector<helpers::DataAggregator<typename Container::value_type> > cluster_aggregates;
        
        std::vector<int> permutation(N);
        for (int i=0; i<N; ++i) permutation[i]=i;

        for (int passcount = 0; passcount < npass; ++passcount) {
            
            // vector of aggregator objects
            // helpers that are :
            // - either a single data value as now when data.aggregate exists
            // - or two data (the sum and avg) + a count when data.aggregate does not exist
            //   then data.operator+=(data) and operator/=(count) are used
            // - can be converted back to data on request
            for (int i=0; i<N; ++i) {
                int ridx = randomProvider(i, N-1);
                std::swap(permutation[i], permutation[ridx]);
            }
            
            for (int ridx = 0; ridx<N; ++ridx) {
                int eltidx = permutation[ridx];
                std::vector<int> matchingcl;
                // O(K.N) algo
                for (int clidx = 0; clidx<(int)allclusters->size(); ++clidx) {
                    if (matchPredicate(cluster_aggregates[clidx],container[eltidx])>0)
                        matchingcl.push_back(clidx);
                }
                // empty ? => new cluster
                if (matchingcl.empty()) {
                    Cluster newCluster(new intVector(1,eltidx));
                    allclusters->push_back(newCluster);
                    // force a copy even when using shared semantics value_type
                    cluster_aggregates.push_back(typename Container::value_type());
                    cluster_aggregates.back().aggregate(container[eltidx]);
                    continue;
                }
                // merge all matching clusters
                // ordered loops above, we can safely remove last elements as we go
                for (int i=(int)matchingcl.size()-1; i>=1; --i) {
                    // merge cluster elements
                    (*allclusters)[matchingcl[0]]->insert((*allclusters)[matchingcl[0]]->end(), (*allclusters)[matchingcl[i]]->begin(), (*allclusters)[matchingcl[i]]->end());
                    // aggregate representative values
                    cluster_aggregates[matchingcl[0]].aggregate(cluster_aggregates[matchingcl[i]]);
                    // remove merged cluster
                    (*allclusters)[matchingcl[i]].swap(allclusters->back());
                    allclusters->pop_back();
                    cluster_aggregates[matchingcl[i]] = cluster_aggregates.back();
                    cluster_aggregates.pop_back();
                }
                (*allclusters)[matchingcl[0]]->push_back(eltidx);
                cluster_aggregates[matchingcl[0]].aggregate(container[eltidx]);
            }
            
            // remove empty clusters that were not matched in this pass
            for (int cliter = 0; cliter<(int)allclusters->size(); ) {
                if ((*allclusters)[cliter]->empty()) {
                    (*allclusters)[cliter].swap(allclusters->back());
                    allclusters->pop_back();
                    cluster_aggregates[cliter] = cluster_aggregates.back();
                    cluster_aggregates.pop_back();
                } else {
                    ++cliter;
                }
            }
            
            // Compensation for small clusters that were left around with no data
            // matching both the small cluster and the main one
            // O(K.(K-1)/2)
            for (int cliter = 0; cliter<(int)allclusters->size(); ++cliter) {
                std::vector<int> matchingcl;
                // match each cluster pair only once
                for (int clidx = cliter+1; clidx<(int)allclusters->size(); ++clidx) {
                    if (matchPredicate(cluster_aggregates[clidx],cluster_aggregates[cliter])>0)
                        matchingcl.push_back(clidx);
                }
                // no match? => next iter
                if (matchingcl.empty()) continue;
                // merge all matching clusters
                // ordered loops above, we can safely remove last elements as we go
                // entries in matchingcl are necessarily above cliter
                for (int i=(int)matchingcl.size()-1; i>=0; --i) {
                    // merge cluster elements
                    (*allclusters)[cliter]->insert((*allclusters)[cliter]->end(), (*allclusters)[matchingcl[i]]->begin(), (*allclusters)[matchingcl[i]]->end());
                    // aggregate representative values
                    cluster_aggregates[cliter].aggregate(cluster_aggregates[matchingcl[i]]);
                    // remove merged cluster
                    (*allclusters)[matchingcl[i]].swap(allclusters->back());
                    allclusters->pop_back();
                    cluster_aggregates[matchingcl[i]] = cluster_aggregates.back();
                    cluster_aggregates.pop_back();
                }
            }

            // prepare for next pass, if any.
            // Clears previous matching elements and redo the matching 
            // with existing clusters so as to be less dependent on the first-pass ordering
            if (passcount!=npass-1) for (int cliter = 0; cliter<(int)allclusters->size(); ++cliter) {
                (*allclusters)[cliter]->clear();
            }

        } // end of the pass loop

        clusterIndices.resize(N);

        for (int clusterIndex = 0; clusterIndex < (int)allclusters->size(); ++clusterIndex) {
            for (typename Cluster::iterator it = (*allclusters)[clusterIndex].begin(); it != (*allclusters)[clusterIndex].end(); ++it) {
                clusterIndices[*it] = clusterIndex;
            }
        }
        
        // Perform cycles of split/merge until all constraints are respected
        // It's up to the user not to provide contradictions, or the loop will never end
        for (int convergenceCounter = 0; convergenceCounter < 10; ++convergenceCounter) {
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

}

#endif
