/*
 *  MEGAHIT
 *  Copyright (C) 2014 The University of Hong Kong
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file functional.h
 * @brief Useful Functionals.
 * @author Yu Peng
 * @version 1.0.0
 * @date 2011-08-24
 */

#ifndef __BASIC_FUNCTIONALS_H_

#define __BASIC_FUNCTIONALS_H_

#include <functional>


/**
 * @brief It is a basic functional to return the input itself.
 *
 * @tparam T
 */
template <typename T>
struct Identity
{
    const T &operator ()(const T &value) const
    { return value; }
};

/**
 * @brief It is a basic functional to select the first element of a stl pair.
 *
 * @tparam Pair
 */
template <typename Pair>
struct Select1st
{
    typedef typename Pair::first_type value_type;
    const value_type &operator ()(const Pair &pair) const
    { return pair.first; }
};

/**
 * @brief It is a basic functional to get the key value from a key-value pair.
 *
 * @tparam Key
 * @tparam Value
 */
template <typename Key, typename Value>
struct GetKey
{
    const Key &operator ()(const Value &value) const
    { return value.key(); }
};

#endif

