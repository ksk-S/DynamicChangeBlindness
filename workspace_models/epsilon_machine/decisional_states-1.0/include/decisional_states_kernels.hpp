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

#ifndef DECISIONAL_STATES_KERNELS_H
#define DECISIONAL_STATES_KERNELS_H

#include <limits>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/random.hpp>
#include <boost/shared_array.hpp>

#include "decisional_states_helpers.hpp"

// For M_LOG2E
#include <cmath>

/* Some useful kernels
*/


namespace decisional_states {


namespace helpers {
template<class T>
inline static T sqr(T x) {return x*x;}
}
/*
    Simple gaussian kernel
    - either one dimensional
    - or multivariate but with diagonal matrix where all elements are the same
    Good enough to not introduce a priori correlations, and much better performance than the full matrix version
    Not normalised, k(x,x) = 1

    Note: Full kernel with covariance matrix was first implemented with uBLAS. Even has Cholesky decomposition for sampling, etc.
          This version was overkill when feeding the kernel with a diagonal matrix as is often the case.
          It was not maintained and shall not compile with current code API.
          It is still available as a git branch.
*/
template<class FloatType = double>
struct SimpleGaussianKernel {
    // log2 of the constant that normalises the distribution
    // multivariate ex: normConst=1/( (2pi)^(D/2) * abs(Sigma)^(1/2))
    // log2Norm = - ( D*log2(2pi) + log2(abs(Sigma)) ) /2
    //FloatType log2Norm;
    FloatType inv2sigma;
    int dim;
    // variance of each individual dimension: the kernel is variance * Identity, with identity
    // having the correct dimension
    // this model with diag matrix of sigma^2 elements is equivalent to a product of
    // independant gaussians along each dimension, with variance_dim each.
    // exp(-0.5 (d_1^2 + d_2^2 + ... + d_n^2) / variance_dim)
    // = exp(-0.5 * d_1^2 / variance_dim -0.5 * d_2^2 / variance_dim + ... - -0.5 * d_n^2 / variance_dim) )
    // = exp(-0.5 * d_1^2 / variance_dim) * exp(-0.5 * d_2^2 / variance_dim) * ... up to n
    FloatType variance;
    FloatType cutoff;
    FloatType precision;
    boost::normal_distribution<FloatType> normal_dist;

    bool undetermined;

    bool valid() {
        return !undetermined;
    }
/*
    void setDeterminant(FloatType det) {
        if (det<=0) undetermined = true;
        else undetermined = false;
        inv2sigma = M_LOG2E * FloatType(-0.5)/ det;
        FloatType variance_dim = (FloatType)std::pow(det, 1.0 / dim);
        variance = variance_dim;
        normal_dist = boost::normal_distribution<FloatType>(FloatType(0),helpers::sqrt(variance_dim));
        cutoff = helpers::log2(precision) / inv2sigma;
    }
*/

    void setSize(FloatType variance_dim) {
        if (variance<=0 || !isfinite(variance)) undetermined = true;
        else undetermined = false;
        //inv2sigma = M_LOG2E * FloatType(-0.5)/ helpers::powN(variance_dim, dim);
        inv2sigma = M_LOG2E * FloatType(-0.5)/ variance_dim;
        variance = variance_dim;
        normal_dist = boost::normal_distribution<FloatType>(FloatType(0),helpers::sqrt(variance_dim));
        cutoff = helpers::log2(precision) / inv2sigma;
    }

    FloatType getSize() {
        return variance;
    }

    int getDimension() {
        return dim;
    }

    SimpleGaussianKernel() {
        undetermined = true;
    }

    // mono-dimensional
    SimpleGaussianKernel(FloatType variance, FloatType t = std::numeric_limits<FloatType>::min()) : dim(1) {
        if (variance<=0) undetermined = true;
        else undetermined = false;
        //log2Norm = helpers::log2(helpers::pi<FloatType>()*2*variance) * FloatType(-0.5);
        inv2sigma = M_LOG2E * FloatType(-0.5)/variance;
        this->variance = variance;
        normal_dist = boost::normal_distribution<FloatType>(FloatType(0),helpers::sqrt(variance));
        precision = t;
        // exp2(res*inv2sigma) < t
        // res*inv2sigma < log2(t)
        // res > log2(t) / inv2sigma  // positive numbers now
        // use t = min representable float number
        cutoff = helpers::log2(t) / inv2sigma;
    }
    // multivariate
    SimpleGaussianKernel(FloatType variance_dim, int _dim, FloatType t = std::numeric_limits<FloatType>::min()) : dim(_dim) {
        if (variance_dim<=0 || !isfinite(variance_dim)) undetermined = true;
        else undetermined = false;
        //log2Norm = variance_dim; for(int i=1; i<dim; ++i) log2Norm *= variance_dim; // determinant of diag matrix
        //log2Norm = (helpers::log2(helpers::pi<FloatType>()*2) * dim + helpers::log2(log2Norm) ) * FloatType(-0.5);
        
        // NO! this is (var * Id) and not (var^N * Id)
        // inv2sigma factorises var^-1 out of distance computation
        //inv2sigma = M_LOG2E * FloatType(-0.5)/ helpers::powN(variance_dim, _dim);

        inv2sigma = M_LOG2E * FloatType(-0.5)/variance_dim;
        variance = variance_dim;
        // normal dist along each dimension
        normal_dist = boost::normal_distribution<FloatType>(FloatType(0),helpers::sqrt(variance_dim));
        precision = t;
        cutoff = helpers::log2(t) / inv2sigma;
    }

    template <typename T>
    inline typename boost::enable_if<boost::is_arithmetic<T>, FloatType>::type distance2(T x, T y) {
        return (x-y) * (x-y);
    }
    template <typename T>
    inline typename boost::disable_if<boost::is_arithmetic<T>, FloatType>::type distance2(const T& x, const T& y) {
        FloatType res=0;
        for(int i=0; i<dim; ++i) res += (x[i] - y[i])*(x[i] - y[i]);
        return res;
    }


    // univariate kernel - for all arithmetic types. Pass by value then
    template <typename T>
    typename boost::enable_if<boost::is_arithmetic<T>, FloatType>::type operator()(T x, T y) {
        assert(dim==1);
        FloatType res = (x-y)*(x-y);
        if (res>cutoff || !isfinite(res)) return 0;
        return helpers::exp2(res*inv2sigma);// + log2Norm);
    }

    // multivariate kernel, for use with arrays and std::vectors...
    // TODO: enable_if has_operator[]
    // in the mean time allow the function for non-arithmetic types only
    template <typename T>
    typename boost::disable_if<boost::is_arithmetic<T>, FloatType>::type operator()(const T& x, const T& y) {
        FloatType res=0;
        for(int i=0; i<dim; ++i) res += (x[i] - y[i])*(x[i] - y[i]);
//std::cout << "d=" << res << ", cutoff=" << cutoff << ", d*inv2sigma=" << res*inv2sigma << ", exp...= " << helpers::exp2(res*inv2sigma) << std::endl;
        if (res>cutoff || !isfinite(res)) return 0;
        return helpers::exp2(res*inv2sigma);// + log2Norm);
    }

    template <typename T, class RNG>
    typename boost::disable_if<boost::is_arithmetic<T> >::type sample(const T& base_point, T& random_point_around_base, RNG& rng) {
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<FloatType> > normal_generator(rng, normal_dist);
        for(int i=0; i<dim; ++i) random_point_around_base[i] = normal_generator() + base_point[i];
    }

    template <typename T, class RNG>
    typename boost::enable_if<boost::is_arithmetic<T> >::type sample(const T base_point, T& random_point_around_base, RNG& rng) {
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<FloatType> > normal_generator(rng, normal_dist);
        random_point_around_base = normal_generator() + base_point;
    }

    /// Must be consistent with a neighbourhood finder
    /// Here we used Squared Euclidian distance
    inline FloatType neighbourhoodSize() {
        return cutoff;
    }

    /// Determines a kernel size from data
    /// The topic is controversial and still subject of research, there is currently no method for automatic setting of kernel size from data
    /// Yet we have several rules of thumb and AMISE estimates for one-dimension, but the problem is even worse for multivariate data...
    /// A cross-validation for searching the kernel that gives max Likelihood of the data solution was once implemented
    /// This may be fine in discrete cases (c.v. otherwise size=0 is a solution),
    /// but in continuous cases this does not give good results in regions of high density
    /// (where removing one sample and counting the c.v. on others does not help with size -> 0).
    /// So, back to the basics:
    /// - we need a size that depends on data: the more data in the same area, the smaller size.
    /// - we need a size that is not too sensitive to outliers, not too oversmoothed and undersmoothed
    /// - we need something reasonably fast to compute, whatever the dimension
    /// => taking the std. dev. as half the median nearest neighbour distance fills all these points.
    /// => taking the std. dev. as half the avg distance (default) is faster, and we can fall back to O(NlogN) median computation anyway if outliers are a problem.
    /// This is certainly not perfect, but still fills all the previous points.
    /// For real scenarios, I'd recommend a parameter search for the best kernel on a logarithmic scale around the default + refine on a regular grid at the best found scale
    template <class Container, class NeighborhoodFinder>
    FloatType setSizeFromSamples(const Container& container, NeighborhoodFinder& finder, bool useMedian = false) {
        /// Start with an approximate kernel size using the sample nearest-neighbor average distance
        finder.setup(container);
        NearestNeighbourFoundAction<Container, NeighborhoodFinder> nfa(container);
        finder.findNearest(container, nfa);
        // avg distance in dim N
        FloatType inisize;
        if (useMedian) inisize = nfa.getMed();
        else inisize = nfa.getAvg();
        inisize = inisize * inisize * 0.25; // variance
        setSize(inisize);
        return inisize;
    }

    template<class Container, class NearestNeighbourFinder>
    struct NearestNeighbourFoundAction {
        const Container& samples;
        boost::shared_array<FloatType> distances;

        NearestNeighbourFoundAction(const Container& c) : samples(c), distances(new FloatType[c.size()]) {
            // flag for missing values where there was no action callback
            for (int i=0; i<(int)c.size(); ++i) distances[i] = -1;
        }

        void operator()(std::size_t index, std::size_t neighbourIndex) {
            distances[index] = NearestNeighbourFinder::distance(samples[index], samples[neighbourIndex]);
        }

        FloatType getAvg() {
            FloatType sumDist = 0;
            int validcount = 0;
            for (int i=0; i<(int)samples.size(); ++i) if (isfinite(distances[i]) && distances[i]>=0) {
                sumDist += distances[i];
                ++validcount;
            }
            return sumDist / (FloatType)validcount; // no count, undefined
        }
        FloatType getMed() {
            // std::sort + unsorted NaN is an issue with current C++ specs apparently
            // thus, remove the NaN first as O(N) before doing the O(NlogN) sort
            int last = samples.size()-1;
            for (int i=0; i<=last;) {
                if (distances[i]<0 || !isfinite(distances[i])) {
                    std::swap(distances[i],distances[last]);
                    --last;
                    continue;
                }
                ++i;
            }
            if (last<0) return distances[0]; // not finite
            std::sort(&distances[0], &distances[last]);
            if ((last&1)==0) return distances[last/2];
            return (distances[last/2] + distances[(last+1)/2]) *0.5f;
        }
    };
};

}

#endif
