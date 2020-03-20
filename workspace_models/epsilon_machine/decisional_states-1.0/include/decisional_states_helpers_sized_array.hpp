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

#ifndef DECISIONAL_STATES_HELPERS_SIZED_ARRAY_H
#define DECISIONAL_STATES_HELPERS_SIZED_ARRAY_H

#include <boost/functional/hash.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>


namespace decisional_states {
namespace helpers {

// Use vectors/lists for inserts, then this for storage
template<class Element>
struct SizedArray {
    Element* array;
    std::size_t array_size;
    // Special check for self-affectation
    SizedArray& operator=(const SizedArray& c) {
        if (c.array==array) return *this; // self-equal
        array_size = c.array_size;
        delete[] array;
        array = new Element[array_size];
        for (int i=0; i<(int)array_size; ++i) array[i] = c.array[i];
        return *this;
    }
    template<class Container>
    SizedArray& operator=(const Container& c) {
        array_size = c.size();
        delete[] array;
        array = new Element[array_size];
        int i = 0;
        for (typename Container::const_iterator it = c.begin(); it != c.end(); ++it) {
            array[i++] = *it;
        }
        return *this;
    }
    // also works for copy constructor
    // except that we redirect arithmetic types to the size_t constructor
    template<class Container>
    SizedArray(const Container& c, typename boost::disable_if<boost::is_arithmetic<Container> >::type* dummy = 0) : array(0), array_size(0) {
        *this = c;
    }
    ~SizedArray() {delete[] array;}
    SizedArray(std::size_t s) : array(new Element[s]), array_size(s) {}
    SizedArray() : array(0), array_size(0) {}
    void reset(std::size_t s) {
        delete[] array;
        array_size = s;
        array = new Element[array_size];
    }
    typedef Element value_type;
    typedef const Element* const_iterator;
    typedef Element* iterator;
    iterator begin() {return array;}
    const_iterator begin() const {return array;}
    iterator end() {return array+array_size;}
    const_iterator end() const {return array+array_size;}
    std::size_t size() const {return array_size;}
    bool empty() const {return array_size==0;}
    template<class Container>
    bool operator==(const Container& c) const {
        if (array_size!=c.size()) return false;
        int i=0;
        for (typename Container::const_iterator it = c.begin(); it != c.end(); ++it) {
            if (array[i++] != *it) return false;
        }
        return true;
    }
    Element& operator[](std::size_t idx) {
        return array[idx];
    }
    const Element& operator[](std::size_t idx) const {
        return array[idx];
    }
    // precondition: the array is sorted
    // effect: element x is inserted at its place, no change if it already exists
    // returns: the index of x.
    // postcondition: the array is sorted and contains x
    // Note: unlike STL containers we do not return pair<iterator,bool>. But one can always compare size before and after the operation.
    template<class Comparator>
    int insert_unique_sorted(const Element& x, Comparator comp) {
        // dichotomic search
        int a = 0, b = array_size;
        while (a!=b) {
            int c = (a + b) / 2;
            // a!=b, c is always valid (possibly a)
            if (comp(array[c],x)) {a = c+1; continue;}
            if (comp(x,array[c])) {b = c; continue;}
            // array[c]==x, or rather, it is neither < nor >
            return c;
        }
        // not found, but b gives the insertion point
        Element* new_array = new Element[++array_size];
        // for really small arrays, faster than function call. Inefficient for large arrays
        for (int i=0; i<b; ++i) new_array[i] = array[i];
        new_array[b] = x;
        for (int i=b+1; i<(int)array_size; ++i) new_array[i] = array[i-1];
        delete[] array;
        array = new_array;
        return b;
    }
    // default comparator: std::less
    int insert_unique_sorted(const Element& x) {
        return insert_unique_sorted(x, std::less<Element>());
    }
    // TODO: operator< if needed, swap
    void swap(SizedArray& other) {
        Element* tmparray;
        std::size_t tmparray_size;
        tmparray = array; array = other.array; other.array = tmparray;
        tmparray_size = array_size; array_size = other.array_size; other.array_size = tmparray_size;
    }
};

template<class Element>
std::size_t hash_value(const SizedArray<Element>& array) {
    size_t seed = 0;
    for(std::size_t i=0; i<array.size(); ++i) boost::hash_combine(seed, array[i]);
    return seed;
}


// utility for main program file
template<class HasFirst>
struct LessOnFirst {
    inline bool operator()(const HasFirst& a,const HasFirst& b) const {return a.first < b.first;}
};
template<class HasFirst>
struct GreaterOnFirst {
    inline bool operator()(const HasFirst& a,const HasFirst& b) const {return a.first > b.first;}
};

template<class Range, class Element>
struct SizedArrayRangeCompare {
    inline bool operator()(const SizedArray<Element>& x, const Range& y) const {return x == y;}
    inline bool operator()(const Range& x, const SizedArray<Element>& y) const {return y == x;}
};
template<class Range, class Element>
struct SizedArrayRangeHash {
    inline std::size_t operator()(const SizedArray<Element>& x) const {return hash_value(x);}
    inline std::size_t operator()(const Range& x) const {
        size_t seed = 0;
        for(typename Range::iterator it = x.begin(); it != x.end(); ++it) boost::hash_combine(seed, *it);
        return seed;
    }
};

}
}

#endif
