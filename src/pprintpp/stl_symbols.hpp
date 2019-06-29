/*
 * MIT License
 *
 * Copyright (c) 2017 Jacek Galowicz
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#pragma once

#ifdef PPRINTPP_AVOID_STL
#include <stddef.h>
#else
#include <cstddef>
#include <type_traits>
#endif


#ifndef PPRINTPP_AVOID_STL
#endif
#include "typelist.hpp"

namespace pprintpp {


#ifndef PPRINTPP_AVOID_STL

using nullptr_t = std::nullptr_t;

template <typename T>
using remove_cv_t = typename std::remove_cv<T>::type;

template <typename A, typename B>
using is_same = std::is_same<A, B>;

template <bool C, typename A, typename B>
using conditional = std::conditional<C, A, B>;

template <typename T>
using remove_ptr = std::remove_pointer<T>;

template <typename T>
using is_int_type = std::is_integral<T>;

#else

using nullptr_t = decltype(nullptr);

template<typename T>
struct remove_c { using type = T; };
template<typename T>
struct remove_c<const T> { using type = T; };

template<typename T>
struct remove_v { using type = T; };
template<typename T>
struct remove_v<volatile T> { using type = T; };

template <typename T>
using remove_cv_t = typename remove_v<
        typename remove_c<T>::type
    >::type;

template <typename A, typename B>
struct is_same { static constexpr bool value {false}; };

template <typename A>
struct is_same<A, A> { static constexpr bool value {true}; };

template <bool C, typename A, typename B>
struct conditional { using type = A; };

template <typename A, typename B>
struct conditional<false, A, B> { using type = B; };

template <typename T>
struct remove_ptr { using type = T; };

template <typename T>
struct remove_ptr<T*> { using type = T; };

template <typename T>
struct is_int_type {
    using ints = typelist::make_t<char, signed char, int, long, long long, unsigned char, unsigned, unsigned long, unsigned long long>;
    static constexpr bool value {typelist::contains<ints, T>::value};
};

#endif

}
