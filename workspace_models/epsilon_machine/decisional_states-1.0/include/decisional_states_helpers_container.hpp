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

#ifndef DECISIONAL_STATES_HELPERS_CONTAINER_H
#define DECISIONAL_STATES_HELPERS_CONTAINER_H

#include "decisional_states_helpers_detect.hpp"

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/pointee.hpp>

namespace decisional_states {
namespace helpers {

DecisionalStates_Detail_checker(size)

template <class Container, int has_size>
struct ContainerSizeSelector {
    // no size() member: use iterator, O(N) operation
    inline static std::size_t size(Container& c) {
        std::size_t ret = 0;
        for (typename Container::iterator it = c.begin(); it != c.end(); ++it) ++ret;
        return ret;
    }
    inline static std::size_t size(const Container& c) {
        std::size_t ret = 0;
        for (typename Container::const_iterator it = c.begin(); it != c.end(); ++it) ++ret;
        return ret;
    }
};
template <class Container>
struct ContainerSizeSelector<Container,1> {
    inline static std::size_t size(Container& c) {
        return c.size();
    }
    inline static std::size_t size(const Container& c) {
        return c.size();
    }
};
template <class Container>
struct ContainerSizeWrapper {
    inline static std::size_t size(Container& c) {
        return ContainerSizeSelector<Container, Check_size<Container, std::size_t (Container::*)()>::found || Check_size<Container, std::size_t (Container::*)() const>::found>::size(c);
    }
    inline static std::size_t size(const Container& c) {
        return ContainerSizeSelector<Container, Check_size<Container, std::size_t (Container::*)() const>::found>::size(c);
    }
};


// C++_0x lambda calculus with unnamed functions or Java would be helpful too
template<class T, class R, R T::*member>
struct ConvertToMember {
    typedef R& result_type;
    R& operator()(T& t) const {return t.*member;}
    const R& operator()(const T& t) const {return t.*member;}
};
template<class T, class R, R T::*member>
struct ConstConvertToMember {
    typedef const R& result_type;
    const R& operator()(const T& t) const {return t.*member;}
};
// for use with multi_index containers that force non-key constness
template<class T, class R, R T::*member>
struct ConvertToMemberDiscardConst {
    typedef R& result_type;
    R& operator()(const T& t) const {return const_cast<R&>(t.*member);}
};

template<class T, class R> struct MemberPointerType {
    typedef R T::*type;
};

template<class Container, class R, typename MemberPointerType<typename Container::value_type,R>::type member>
struct ContainerMemberAdaptor {
    typedef ConvertToMember<typename Container::value_type,R,member> Converter;
    typedef ConstConvertToMember<const typename Container::value_type,R,member> ConstConverter;
    Container& container;
    ContainerMemberAdaptor(Container& _container) : container(_container) {}
    // must have operator[]. TODO: simulate with O(N) loop
    R& operator[](int i) { return container[i].*member; }
    const R& operator[](int i) const { return container[i].*member; }
    inline std::size_t size() {return ContainerSizeWrapper<Container>::size(container);}
    inline std::size_t size() const {return ContainerSizeWrapper<Container>::size(container);}
    typedef boost::transform_iterator<Converter, typename Container::iterator> iterator;
    typedef boost::transform_iterator<ConstConverter, typename Container::const_iterator> const_iterator;
    const_iterator begin() const {return boost::make_transform_iterator(container.begin(), ConstConverter());}
    const_iterator end() const {return boost::make_transform_iterator(container.end(), ConstConverter());}
    iterator begin() {return boost::make_transform_iterator(container.begin(), Converter());}
    iterator end() {return boost::make_transform_iterator(container.end(), Converter());}
    typedef R value_type;
    typedef R& reference;
    typedef const R& const_reference;
    typedef R* pointer;
    typedef typename boost::iterator_difference<iterator>::type difference_type;
};

template<class T, class R, typename MemberPointerType<typename boost::pointee<T>::type, R>::type member>
struct ConvertToPointedMember {
    typedef R& result_type;
    R& operator()(T& t) const {return t->*member;}
    const R& operator()(const T& t) const {return t->*member;}
};
template<class T, class R, typename MemberPointerType<typename boost::pointee<T>::type, R>::type member>
struct ConstConvertToPointedMember {
    typedef const R& result_type;
    const R& operator()(const T& t) const {return t->*member;}
};

template<class Container, class R, typename MemberPointerType<typename boost::pointee<typename Container::value_type>::type, R>::type member>
struct PointerContainerMemberAdaptor {
    typedef ConvertToPointedMember<typename Container::value_type,R,member> Converter;
    typedef ConstConvertToPointedMember<typename Container::value_type,R,member> ConstConverter;
    Container& container;
    PointerContainerMemberAdaptor(Container& _container) : container(_container) {}
    // must have operator[]. TODO: simulate with O(N) loop
    R& operator[](int i) { return container[i]->*member; }
    const R& operator[](int i) const { return container[i]->*member; }
    inline std::size_t size() {return ContainerSizeWrapper<Container>::size(container);}
    inline std::size_t size() const {return ContainerSizeWrapper<Container>::size(container);}
    typedef typename boost::transform_iterator<Converter, typename Container::iterator> iterator;
    typedef typename boost::transform_iterator<ConstConverter, typename Container::const_iterator> const_iterator;
    const_iterator begin() const {return boost::make_transform_iterator(container.begin(), ConstConverter());}
    const_iterator end() const {return boost::make_transform_iterator(container.end(), ConstConverter());}
    iterator begin() {return boost::make_transform_iterator(container.begin(), Converter());}
    iterator end() {return boost::make_transform_iterator(container.end(), Converter());}
    typedef R value_type;
    typedef R& reference;
    typedef const R& const_reference;
    typedef R* pointer;
    typedef typename boost::iterator_difference<iterator>::type difference_type;
};


template<class Container, class R, class Converter>
struct ContainerAdaptor : public Converter {
    Container& container;
    ContainerAdaptor(Container& _container) : container(_container) {}
    R& operator[](int i) { return (*this)(container[i]); }
    const R& operator[](int i) const { return (*this)(container[i]); }
    inline std::size_t size() {return ContainerSizeWrapper<Container>::size(container);}
    inline std::size_t size() const {return ContainerSizeWrapper<Container>::size(container);}
    typedef boost::transform_iterator<Converter, typename Container::iterator> iterator;
    typedef boost::transform_iterator<Converter, typename Container::const_iterator> const_iterator;
    const_iterator begin() const {return boost::make_transform_iterator(container.begin(), (*this));}
    const_iterator end() const {return boost::make_transform_iterator(container.end(), (*this));}
    iterator begin() {return boost::make_transform_iterator(container.begin(), (*this));}
    iterator end() {return boost::make_transform_iterator(container.end(), (*this));}
    typedef R value_type;
    typedef R& reference;
    typedef const R& const_reference;
    typedef R* pointer;
};

template<class Content>
struct SingleElementContainerWrapper {
    typedef Content value_type;
    typedef Content* iterator;
    Content* content;
    iterator begin() const {return content;}
    iterator end() const {return content + 1;}
    int size() const {return 1;}
    value_type& operator[](int idx) const {return *content;}
    SingleElementContainerWrapper(Content& c) : content(&c) {}
    SingleElementContainerWrapper(const Content& c) : content(const_cast<Content*>(&c)) {}
};

template<class Content>
struct SingleElementContainerCopy {
    typedef Content value_type;
    typedef Content* iterator;
    Content content;
    iterator begin() {return &content;}
    iterator end() {return (&content) + 1;}
    int size() {return 1;}
    value_type& operator[](int idx) {return content;}
    SingleElementContainerCopy(const Content& c) : content(c) {}
};

}
}

#endif


