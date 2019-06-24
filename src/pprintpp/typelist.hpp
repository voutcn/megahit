/*
 * MIT License
 *
 * Copyright (c) 2016 Jacek Galowicz
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

namespace typelist {

struct null_t {
    using head = null_t;
    using tail = null_t;
};

template <class T, class U>
struct tl
{
    using head = T;
    using tail = U;
};

template <class TList>
using head_t = typename TList::head;

template <class TList>
using tail_t = typename TList::tail;


template <class ... Ts> struct make;

template <class T, class ... REST>
struct make<T, REST...> { using type = tl<T, typename make<REST...>::type>; };
template <>
struct make<> { using type = null_t; };

template <class ... Ts>
using make_t = typename make<Ts...>::type;


template <class TList, class T>
struct append;
template <>
struct append<null_t, null_t> { using type = null_t; };
template <class T>
struct append<null_t, T> { using type = make_t<T>; };
template <class Head, class T>
struct append<null_t, tl<Head, T>> { using type = tl<Head, T>; };
template <class Head, class Tail, class T>
struct append<tl<Head, Tail>, T>
{ using type = tl<Head, typename append<Tail, T>::type>; };

template <class TList, class T>
using append_t = typename append<TList, T>::type;

template <class TList, class T>
struct contains;

template <class T>
struct contains<null_t, T> { static constexpr bool value {false}; };
template <class T, class L>
struct contains<tl<T, L>, T> { static constexpr bool value {true}; };
template <class S, class T, class L>
struct contains<tl<S, L>, T> : contains<L, T> {};


template <class TList, class T>
struct remove;

template <class T>
struct remove<null_t, T> { using type = null_t; };
template <class T, class L>
struct remove<tl<T, L>, T> { using type = typename remove<L, T>::type; };
template <class S, class T, class L>
struct remove<tl<S, L>, T> { using type = tl<S, typename remove<L, T>::type>; };

template <class TL, class T>
using remove_t = typename remove<TL, T>::type;


template <class TList, class T, class TS>
struct substitute;

template <class T, class TS>
struct substitute<null_t, T, TS> { using type = null_t; };
template <class T, class L, class TS>
struct substitute<tl<T, L>, T, TS> { using type = tl<TS, typename substitute<L, T, TS>::type>; };
template <class S, class T, class L, class TS>
struct substitute<tl<S, L>, T, TS> { using type = tl<S, typename substitute<L, T, TS>::type>; };

template <class TL, class T, class TS>
using substitute_t = typename substitute<TL, T, TS>::type;

} // namespace tl
