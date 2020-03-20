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

#ifndef DECISIONAL_STATES_HELPERS_MATH_H
#define DECISIONAL_STATES_HELPERS_MATH_H

#include <math.h>
#include <stdlib.h>

namespace decisional_states {
namespace helpers {

// In helpers namespace, we can reuse the base function name
#define DecisionalStates_Detail_math_fun_adaptor(fun,name) \
template<typename FloatType> \
inline FloatType name(FloatType x) { \
    return ::fun(x); \
} \
template<> \
inline float name<float>(float x) { \
    return ::fun ## f(x); \
} \
template<> \
inline long double name<long double>(long double x) { \
    return ::fun ## l(x); \
}

DecisionalStates_Detail_math_fun_adaptor(log2,log2std)
DecisionalStates_Detail_math_fun_adaptor(exp2,exp2std)
DecisionalStates_Detail_math_fun_adaptor(sqrt,sqrt)
DecisionalStates_Detail_math_fun_adaptor(fabs,fabs)
DecisionalStates_Detail_math_fun_adaptor(floor,floor)

template<> inline int fabs<int>(int x) { return ::abs(x); }
template<> inline long fabs<long>(long x) { return ::labs(x); }
template<> inline long long fabs<long long>(long long x) { return ::llabs(x); }

// fast exponentiation by positive integer
template<typename FloatType>
FloatType powN(FloatType x, int N) {
/*
    if (N==0) return 1;
    if (N==1) return x;
    if (N&1 == 0) return powN(x*x, N/2);
    return powN(x*x, N/2) * x;
*/
// actually, with modern architectures, the recursive call cost + if tests are very inefficient for small N!
// whereas loops in L1 instruction cache, + hardware-optimised multiplications on FPU on very few cycles...
// ... where is the boundary at which the recursive version is more efficient? N quite large I guess...
FloatType ret = 1;
for (int i=0; i<N; ++i) ret *= x;
return ret;
}


template<typename FloatType>
inline FloatType exp2(FloatType x) {
#ifndef NO_X86_ASSEMBLY
    FloatType ret;
    // from gnu libm, stripped down to cover only relevant cases here
    asm(
    "fld     %%st(0)\n"
    "frndint\n"                    // int(x), current rounding mode doesn't matter
    "fsubr   %%st(0),%%st(1)\n"    // fract(x), between -1.0 and 1.0, OK for f2xm1
    "fxch\n"
    "f2xm1\n"                      // 2^fract(x) - 1
    "fld1\n"
    "faddp\n"                      // 2^fract(x)
    "fscale\n"                     // 2^x
    "fstp    %%st(1)\n"
    : "=t" (ret)
    : "0" (x)
    : "st(1)"
    );
    return ret;
#else
    return exp2std(x);
#endif
}

template<typename FloatType>
inline FloatType log2(FloatType x) {
#ifndef NO_X86_ASSEMBLY
    FloatType ret;
    // custom implementation
    asm(
    "fld1\n"
    "fxch\n"
    "fyl2x"
    : "=t" (ret)
    : "0" (x)
    : "st(1)"
    );
    return ret;
#else
    return log2std(x);
#endif
}

template<typename FloatType>
inline FloatType pi() {
#ifdef M_PI
    return M_PI;
#else
    return 3.14159265358979323846;
#endif
}

template<>
inline long double pi<long double>() {
#ifdef M_PIl
    return M_PIl;
#else
    return 3.1415926535897932384626433832795029L;
#endif
}

bool ChiSquareMatch(unsigned int ndegree, float chisq, float matchLevel) {
    // See below, igamma(a,0) = gamma(a), then normalize, nig = 1.0f
    if (chisq<1e-6f) return true;

    float x = 0.5f*chisq;
    float negxdivln2;       // set below
    float l2gamma, nig;     // log-gamma and normalized incomplete gamma

    // Custom from-scratch normalized incomplete gamma Q(0.5f*ndegree, 0.5f*chisq) implementation
    // that is faster than usual power-series approx for the special cases of this program.
    // Iterate using the formulas:
    // - Q(a,x) = igamma(a,x) / gamma(a);
    // - twice_a = ndegree is integer.
    // - igamma(1/2, x) = sqrt(pi).erfc(sqrt(x)), gamma(1/2) = sqrt(pi)
    // - igamma(1) = e^-x, gamma(1) = 1
    // - igamma(a, x) = (a-1)*igamma(a-1,x) + x^(a-1).e^-x  (for a>=3/2)
    // - gamma(a) = (a-1)*gamma(a-1) (for a>=3/2)
    unsigned int twice_a = ndegree & 1;
    if (twice_a) {
        // twice a is odd, start with gamma(1/2)
        nig = ::erfcf(::sqrtf(x));
        if (ndegree==1) return nig > matchLevel;
        l2gamma = 0.825748064739; // = log2(sqrt(pi))
        // loop constant
        negxdivln2 = x / (-0.69314718056f);
    } else {
        // start with gamma(1)
        negxdivln2 = x / (-0.69314718056f);
        nig = helpers::exp2(negxdivln2);
        if (ndegree==2) return nig > matchLevel;
        l2gamma = 0.0f;
    }
    float halflog2x = 0.5f * helpers::log2<float>(x); // constant through the loop

    // now iterate using the above formulas
    for (twice_a += 2; twice_a <= ndegree; twice_a += 2) {
        l2gamma += helpers::log2<float>(twice_a) - 1.0f;
        float sumTerm = helpers::exp2(twice_a * halflog2x + negxdivln2 - l2gamma);
        nig += sumTerm;
        // sum has converged
        if (sumTerm < 1e-6) return nig > matchLevel;
        // early break, since we're always adding >0 numbers
        if (nig > matchLevel) return true;
    }
    return false; // since true would have been detected by the loop
}


}
}

#endif
