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

#ifndef DECISIONAL_STATES_HELPERS_DETECT_H
#define DECISIONAL_STATES_HELPERS_DETECT_H

namespace decisional_states {
namespace helpers {

// Technique for detecting the existence of member function, inspired from
//   http://www.codeproject.com/KB/architecture/Detector.aspx
// Except that here we detect only what we need, resulting in simpler code

#define DecisionalStates_Detail_checker(name) \
template <class BaseClass, typename FunctionType> \
struct Check_ ## name { \
    template <FunctionType> struct Finder; \
    template <class Base> static long sfinaeOverload(Finder<&Base::name> *); \
    template <class Base> static char sfinaeOverload(...); \
    enum { found = sizeof(sfinaeOverload<BaseClass>(0)) == sizeof(long) }; \
};

// if "name" type exist in class T
// then define wrapper name = T::name
// otherwise wrapper name = fallback
#define DecisionalStates_Detail_type_wrapper(name,fallback) \
template <class BaseClass> \
struct TypeCheck_ ## name { \
    template <class Base> static long sfinaeOverload(typename Base::name *); \
    template <class Base> static char sfinaeOverload(...); \
    enum { found = sizeof(sfinaeOverload<BaseClass>(0)) == sizeof(long) }; \
}; \
template<class T, int x> struct TypeSelector_ ## name { \
    typedef fallback name; \
}; \
template<class T> struct TypeSelector_ ## name <T,1> { \
    typedef typename T::name name; \
}; \
template<class T> \
struct name ## Wrapper { \
    enum { found = TypeCheck_ ## name<T>::found }; \
    typedef typename TypeSelector_ ## name<T, found>::name name; \
};

}

}

#endif


