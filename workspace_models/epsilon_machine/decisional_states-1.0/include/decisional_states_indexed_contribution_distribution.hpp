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

#ifndef DECISIONAL_STATES_INDEXED_CONTRIBUTION_DISTRIBUTION_H
#define DECISIONAL_STATES_INDEXED_CONTRIBUTION_DISTRIBUTION_H

#include "decisional_states_helpers_sized_array.hpp"
#include "decisional_states_helpers_math.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

namespace decisional_states {

template<typename FloatType>
struct IndexedContributionDistributionBase {
    IndexedContributionDistributionBase() : sumcontribs(0), info(0) {}

    // array of indices into samples
    // shared with other dist with the same support
    boost::shared_ptr<helpers::SizedArray<int> > idxarray;
    // array of p values at sample points, same size as idx list
    // shared so Distribution can be fast-copied
    boost::shared_array<FloatType> contributions;
    FloatType sumcontribs;
    void* info;
};


/// Distribution object use by kernel distribution managers
/// This class was made independant of the KDM type
template<typename FloatType>
struct IndexedContributionDistribution : public boost::shared_ptr<IndexedContributionDistributionBase<FloatType> > {

    IndexedContributionDistribution() : boost::shared_ptr<IndexedContributionDistributionBase<FloatType> >(new IndexedContributionDistributionBase<FloatType>()) {}

    // Allows clustering algorithm to fusion distributions
    // Aggregate the kernel contributions
    void aggregate(const IndexedContributionDistribution& d) {
        if ((*this)->idxarray.get()==0) {
            (*this)->idxarray = d->idxarray; // shallow copy OK here, the array is never modified after creation
            (*this)->contributions = d->contributions;
            (*this)->sumcontribs = d->sumcontribs;
            // info unique to original distribution
            return;
        }
        // keeping sumcontribs allow summing everything without bothering and thus aggregating in a row
        // So: aggregate( d1, aggregate(d2, d3)) works fine whatever the order of d1, d2, d3
        std::map<int, FloatType> mergedContribs;
        // merge all entries in the map
        for (int i=0; i<(int)(*this)->idxarray->size(); ++i) mergedContribs[(*(*this)->idxarray)[i]] += (*this)->contributions[i];
        for (int i=0; i<(int)d->idxarray->size(); ++i) mergedContribs[(*d->idxarray)[i]] += d->contributions[i];
        // copy the map back to our own arrays
        std::size_t numentries = mergedContribs.size();
        // Do NOT store this array in the global sampleidx container - aggregation of other dists would uselessly fill
        // that container with temporary objects that would hog memory until the whole thing gets deallocated
        (*this)->idxarray = boost::shared_ptr<helpers::SizedArray<int> >(new helpers::SizedArray<int>(numentries));
        (*this)->contributions = boost::shared_array<FloatType>(new FloatType[numentries]);
        int i = 0;
        for (typename std::map<int,FloatType>::iterator it = mergedContribs.begin(); it != mergedContribs.end(); ++it, ++i) {
            (*(*this)->idxarray)[i] = it->first;
            (*this)->contributions[i] = it->second;
        }
        (*this)->sumcontribs += d->sumcontribs;
    }




    // A few "distance" measures between distributions

    // Bhattacharyya Distance B(P||Q) = sqrt(1 - sum_i sqrt(p_i * q_i))
    // when all elements are equal, B(P||Q) = 0
    // matchlevel = 1-x
    // sqrt(1 - sum_i sqrt(p_i * q_i)) < x
    // 1 - x*x < sum_i sqrt(p_i * q_i)
    // But: 1 - (1-a)^2 = 1 - (1 - 2a + a^2) = 2a - a^2
    struct BhattacharyyaMatcher {
        // first field in struct
        FloatType matchlevel;

        BhattacharyyaMatcher(FloatType _matchlevel = 0.95) : matchlevel(_matchlevel * 2 - _matchlevel * _matchlevel) {}
        FloatType operator()(const IndexedContributionDistribution& a, const IndexedContributionDistribution& b) const {
            if (a->idxarray->empty() || b->idxarray->empty()) return -1; // no match
            std::size_t ai = 0, bi = 0;
            FloatType bhat = 0;
            while(true) {
                if ((*a->idxarray)[ai] < (*b->idxarray)[bi]) {
                    ++ai; if (ai>=a->idxarray->size()) break;
                    continue;
                }
                if ((*b->idxarray)[bi] < (*a->idxarray)[ai]) {
                    ++bi; if (bi>=b->idxarray->size()) break;
                    continue;
                }
                // both match, we can add non-null sqrt
                bhat += helpers::sqrt(a->contributions[ai] * b->contributions[bi]);
                ++ai; if (ai>=a->idxarray->size()) break;
                ++bi; if (bi>=b->idxarray->size()) break;
            }
            // if more elements remain, ignore and do not add 0
            // return >0 is match
            // b = sum_i sqrt((pza_i/sumca) * (pzb_i/sumcb))
            // b = sum_i sqrt(pza_i * pzb_i/(sumca*sumcb))
            // b = sum_i sqrt(pza_i * pzb_i)/sqrt(sumca*sumcb)
            // b = (sum_i sqrt(pza_i * pzb_i) )/sqrt(sumca*sumcb)
            // b > matchlevel => sum_i sqrt(pza_i * pzb_i) > matchlevel * sqrt(sumca*sumcb)
            return bhat - matchlevel * helpers::sqrt(a->sumcontribs * b->sumcontribs);
        }
    };

    // Known as "Variational Distance": V(P||Q) = sum_i |p_i - q_i|
    // when all elements are equal, V(P||Q) = 0
    struct VariationalMatcher {
        FloatType distmax;
        VariationalMatcher(FloatType _matchlevel = 0.95) : distmax(1 - _matchlevel) {}

        FloatType operator()(const IndexedContributionDistribution& a, const IndexedContributionDistribution& b) const {
            if (a->idxarray->empty() || b->idxarray->empty()) return -1; // no match
            FloatType v = 0;
            std::size_t ai = 0, bi = 0;
            while(true) {
                if ((*a->idxarray)[ai] < (*b->idxarray)[bi]) {
                    v += a->contributions[ai];
                    ++ai; if (ai>=a->idxarray->size()) break;
                    continue;
                }
                if ((*b->idxarray)[bi] < (*a->idxarray)[ai]) {
                    v += b->contributions[bi];
                    ++bi; if (bi>=b->idxarray->size()) break;
                    continue;
                }
                v += helpers::fabs(a->contributions[ai] - b->contributions[bi]);
                ++ai; if (ai>=a->idxarray->size()) break;
                ++bi; if (bi>=b->idxarray->size()) break;
            }
            // either a or b reached the end. add the other contributions
            while (ai<a->idxarray->size()) v += a->contributions[ai++];
            while (bi<b->idxarray->size()) v += b->contributions[bi++];
            // return >0 is match
            return distmax - v;
        }
    };

    // Known as "Harmonic mean Distance": M(P||Q) = 2 * sum_i p_i*q_i/(p_i+q_i)
    // when all elements are equal, M(P||Q) = 1
    struct HarmonicMatcher {
        // not a probability, just a match level
        FloatType matchlevel;
        HarmonicMatcher(FloatType _matchlevel = 0.95) : matchlevel(_matchlevel * 2) {}

        FloatType operator()(const IndexedContributionDistribution& a, const IndexedContributionDistribution& b) const {
            if (a->idxarray->empty() || b->idxarray->empty()) return -1; // no match
            std::size_t ai = 0, bi = 0;
            FloatType harmonicMean = 0;
            while(true) {
                if ((*a->idxarray)[ai] < (*b->idxarray)[bi]) {
                    ++ai; if (ai>=a->idxarray->size()) break;
                    continue;
                }
                if ((*b->idxarray)[bi] < (*a->idxarray)[ai]) {
                    ++bi; if (bi>=b->idxarray->size()) break;
                    continue;
                }
                // both match, we can add non-null sqrt
                harmonicMean += a->contributions[ai] * b->contributions[bi] / (a->contributions[ai] + b->contributions[bi]);
                ++ai; if (ai>=a->idxarray->size()) break;
                ++bi; if (bi>=b->idxarray->size()) break;
            }
            return harmonicMean - matchlevel;
        }
    };

    // Jensen-Shannon divergence (aka symmetrised Kullback-Leibler divergence)
    // JS(P||Q) = ( KL(P||(P+Q)/2) + KL(Q||(P+Q)/2) ) / 2
    // JS(P||Q) = 0 if both match
    struct JensenShannonMatcher {
        FloatType distmax;
        JensenShannonMatcher(FloatType _matchlevel = 0.95) : distmax(1 - _matchlevel) {}

        FloatType operator()(const IndexedContributionDistribution& a, const IndexedContributionDistribution& b) const {
            if (a->idxarray->empty() || b->idxarray->empty()) return -1; // no match
            FloatType res = 0;
            std::size_t ai = 0, bi = 0;
            while(true) {
                if ((*a->idxarray)[ai] < (*b->idxarray)[bi]) {
                    res += a->contributions[ai];
                    ++ai; if (ai>=a->idxarray->size()) break;
                    continue;
                }
                if ((*b->idxarray)[bi] < (*a->idxarray)[ai]) {
                    res += b->contributions[bi];
                    ++bi; if (bi>=b->idxarray->size()) break;
                    continue;
                }
                FloatType m = (a->contributions[ai] + b->contributions[bi]) * 0.5;
                res += a->contributions[ai] * helpers::log2(a->contributions[ai] / m) + b->contributions[bi] * helpers::log2(b->contributions[bi] / m);
                ++ai; if (ai>=a->idxarray->size()) break;
                ++bi; if (bi>=b->idxarray->size()) break;
            }
            // either a or b reached the end. add the other contributions
            while (ai<a->idxarray->size()) res += a->contributions[ai++];
            while (bi<b->idxarray->size()) res += b->contributions[bi++];
            return distmax - res * 0.5;
        }
    };

    struct ChiSquareMatcher {
        FloatType matchLevel;
        ChiSquareMatcher(FloatType _matchlevel = 0.95) : matchLevel(_matchlevel) {}

        FloatType operator()(const IndexedContributionDistribution& a, const IndexedContributionDistribution& b) const {
            if (a->idxarray->empty() || b->idxarray->empty()) return -1; // no match
            FloatType n1 = a->idxarray->size();
            FloatType n2 = b->idxarray->size();
            FloatType n1n2 = n1 * n2;
            FloatType n1overn2 = n1 / n2;
            FloatType n2overn1 = n2 / n1;
            FloatType chisq = 0;
            unsigned int ndegree = 0;
            std::size_t ai = 0, bi = 0;
            while(true) {
                ++ndegree;
                if ((*a->idxarray)[ai] < (*b->idxarray)[bi]) {
                    chisq += a->contributions[ai] * n2overn1;
                    ++ai; if (ai>=a->idxarray->size()) break;
                    continue;
                }
                if ((*b->idxarray)[bi] < (*a->idxarray)[ai]) {
                    chisq += b->contributions[bi] * n1overn2;
                    ++bi; if (bi>=b->idxarray->size()) break;
                    continue;
                }
                FloatType tmp = n2 * a->contributions[ai] - n1 * b->contributions[bi];
                chisq += tmp * tmp / (n1n2 * (a->contributions[ai] + b->contributions[bi]));
                ++ai; if (ai>=a->idxarray->size()) break;
                ++bi; if (bi>=b->idxarray->size()) break;
            }
            // either a or b reached the end. add the other contributions
            while (ai<a->idxarray->size()) {++ndegree; chisq += a->contributions[ai++] * n2overn1;}
            while (bi<b->idxarray->size()) {++ndegree; chisq += b->contributions[bi++] * n1overn2;}
            return helpers::ChiSquareMatch(ndegree, chisq, matchLevel);
        }
    };




};

}

#endif

