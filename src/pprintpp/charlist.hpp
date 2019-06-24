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


#include "typelist.hpp"

namespace charlist {

using namespace typelist;

template <char val> struct char_t { static const constexpr char value {val}; };

template <char c, char ... chars>
struct char_tl;

template <char c, char ... chars>
struct char_tl {
    using type = typelist::tl<char_t<c>, typename char_tl<chars...>::type>;
};
template <char c>
struct char_tl<c> {
    using type = typelist::tl<char_t<c>, typelist::null_t>;
};

template <char ... chars>
using char_tl_t = typename char_tl<chars...>::type;

template <class Str, size_t Pos, char C>
struct string_list;

template <class Str, size_t Pos, char C>
struct string_list {
    using next_piece = typename string_list<
                            Str,
                            Pos + 1,
                            Str::str()[Pos + 1]
                        >::type;
    using type = typelist::tl<char_t<C>, next_piece>;
};

template <class Str, size_t Pos>
struct string_list<Str, Pos, '\0'> {
    using type = typelist::null_t;
};

template <class Str>
using string_list_t = typename string_list<Str, 0, Str::str()[0]>::type;

template <typename typelist, char ... chars>
struct tl_to_varlist;

template <char c, typename restlist, char ... chars>
struct tl_to_varlist<typelist::tl<char_t<c>, restlist>, chars...>
    : public tl_to_varlist<restlist, chars..., c>
{ };

template <>
struct tl_to_varlist<typelist::null_t> {
    static const char * const str() { return ""; }
};
template <char ... chars>
struct tl_to_varlist<typelist::null_t, chars...> {
    using list = char_tl<chars...>;

    static const char * const str() {
        static constexpr const char string[] = {chars..., '\0'};
        return string;
    }
};

}
