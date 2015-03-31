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

