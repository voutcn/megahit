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
 * @file hash.h
 * @brief Hash Functionals.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-24
 */

#ifndef __BASIC_HASH_H_

#define __BASIC_HASH_H_

#include <stdint.h>

#include <functional>

template <typename T>
struct Hash: public std::unary_function<T, uint64_t>
{
    uint64_t operator ()(const T &value) const
    { return value.hash(); }
};

template <>
struct Hash<int8_t>: public std::unary_function<int8_t, uint64_t>
{
    uint64_t operator ()(int8_t value) const
    { return (uint64_t(value)* 1299709 + 104729) % 323780508946331ULL; }
};

template <>
struct Hash<uint8_t>: public std::unary_function<uint8_t, uint64_t>
{
    uint64_t operator ()(uint8_t value) const
    { return (uint64_t(value)* 1299709 + 104729) % 323780508946331ULL; }
};

template <>
struct Hash<int16_t>: public std::unary_function<int16_t, uint64_t>
{
    uint64_t operator ()(int16_t value) const
    { return (uint64_t(value)* 1299709 + 104729) % 323780508946331ULL; }
};

template <>
struct Hash<uint16_t>: public std::unary_function<uint16_t, uint64_t>
{
    uint64_t operator ()(uint16_t value) const
    { return (uint64_t(value)* 1299709 + 104729) % 323780508946331ULL; }
};

template <>
struct Hash<int32_t>: public std::unary_function<int32_t, uint64_t>
{
    uint64_t operator ()(int32_t value) const
    { return (uint64_t(value)* 1299709 + 104729) % 323780508946331ULL; }
};

template <>
struct Hash<uint32_t>: public std::unary_function<uint32_t, uint64_t>
{
    uint64_t operator ()(uint32_t value) const
    { return (uint64_t(value)* 1299709 + 104729) % 323780508946331ULL; }
};

template <>
struct Hash<int64_t>: public std::unary_function<int64_t, uint64_t>
{
    uint64_t operator ()(int64_t value) const
    { return (uint64_t(value)* 1299709 + 104729) % 323780508946331ULL; }
};

template <>
struct Hash<uint64_t>: public std::unary_function<uint64_t, uint64_t>
{
    uint64_t operator ()(uint64_t value) const
    { return (uint64_t(value)* 1299709 + 104729) % 323780508946331ULL; }
};

#endif

