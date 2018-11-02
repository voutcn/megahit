/**
 * @file histgram.h
 * @brief Histgram Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-13
 * @modified by Dinghua Li
 */

#ifndef __BASIC_HISTGRAM_H__

#define __BASIC_HISTGRAM_H__

#include <stdint.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <map>


/**
 * @brief It is a class which contains a set of number for drawing histgram.
 *
 * @tparam T
 */
template <class T>
class Histgram {
  public:
    typedef T value_type;

    Histgram(): size_(0) {
        lock_ = 0;
    }
    Histgram(const Histgram<T> &hist)
        : map_(hist.map_), size_(hist.size_) {
        lock_ = 0;
    }
    ~Histgram() { }

    const Histgram<T> &operator =(const Histgram<T> &hist) {
        map_ = hist.map_;
        size_ = hist.size_;
        return *this;
    }

    void insert(value_type value, size_t count = 1) {
        while (__sync_lock_test_and_set(&lock_, 1)) while (lock_);

        map_[value] += count;
        size_ += count;
        __sync_lock_release(&lock_);
    }

    size_t count(value_type value) const {
        typename std::map<value_type, size_t>::const_iterator iter = map_.find(value);

        if (iter == map_.end())
            return 0;
        else
            return iter->second;
    }

    size_t count(value_type from, value_type to) const {
        if (from >= to)
            return 0;

        size_t sum = 0;

        for (typename std::map<value_type, size_t>::const_iterator iter = map_.lower_bound(from);
                iter != map_.end() && iter->first < to; ++iter)
            sum += iter->second;

        return sum;
    }

    value_type minimum() const {
        return empty() ? 0 : map_.begin()->first;
    }
    value_type maximum() const {
        return empty() ? 0 : map_.rbegin()->first;
    }

    double mean() const {
        if (size() == 0) {
            return 0;
        }

        return sum() / size();
    }

    double sum() const {
        double sum_ = 0;

        for (typename std::map<value_type, size_t>::const_iterator iter = map_.begin(); iter != map_.end(); ++iter)
            sum_ += 1.0 * iter->first * iter->second;

        return sum_;
    }

    double variance() const {
        double sum = 0;
        double square_sum = 0;
        size_t n = size();

        for (typename std::map<value_type, size_t>::const_iterator iter = map_.begin(); iter != map_.end(); ++iter) {
            sum += 1.0 * iter->first * iter->second;
            square_sum += 1.0 * iter->first * iter->first * iter->second;
        }

        return square_sum / n - (sum / n) * (sum / n);
    }

    double sd() const {
        return std::sqrt(variance());
    }

    value_type median() const {
        size_t half = size() / 2;
        size_t sum = 0;

        for (typename std::map<value_type, size_t>::const_iterator iter = map_.begin(); iter != map_.end(); ++iter) {
            sum += iter->second;

            if (sum > half)
                return iter->first;
        }

        return 0;
    }

    value_type percentile(double p) const {
        size_t half = size() * p;
        size_t sum = 0;

        for (typename std::map<value_type, size_t>::const_iterator iter = map_.begin(); iter != map_.end(); ++iter) {
            sum += iter->second;

            if (sum > half)
                return iter->first;
        }

        return 0;
    }

    double up_mean(double p) const {
        size_t count = size() * p;
        size_t sum = 0;
        double total = 0;

        for (typename std::map<value_type, size_t>::const_reverse_iterator iter = map_.rbegin(); iter != map_.rend(); ++iter) {
            sum += iter->second;
            total += iter->second * iter->first;

            if (sum > count)
                break;
        }

        return (sum != 0) ? total / sum : 0;
    }

    value_type Nx(double x) {
        double total = 0;

        for (typename std::map<value_type, size_t>::const_reverse_iterator iter = map_.rbegin(); iter != map_.rend(); ++iter) {
            total += iter->second * iter->first;

            if (total >= x)
                return iter->first;
        }

        return 0;
    }

    /** Return the first local minimum or zero if a minimum is not
     * found. */
    value_type FirstLocalMinimum() const
    {
        const size_t kSmoothing = 4;
        auto minimum = map_.begin();
        size_t count = 0;
        for (auto it = map_.begin(); it != map_.end(); ++it) {
            if (it->second <= minimum->second) {
                minimum = it;
                count = 0;
            } else if (++count >= kSmoothing)
                break;
        }
        if (minimum->first == maximum())
            return 0;
        return minimum->first;
    }

    size_t Trim(double fraction) {
        size_t trim_size = size_t(size() * fraction / 2 + 0.5);
        std::deque<value_type> trim_values;

        size_t sum = 0;

        for (typename std::map<value_type, size_t>::iterator iter = map_.begin(); iter != map_.end(); ++iter) {
            if (sum + iter->second <= trim_size) {
                sum += iter->second;
                trim_values.push_back(iter->first);
            }
            else
                break;
        }

        size_ -= sum;

        sum = 0;

        for (typename std::map<value_type, size_t>::reverse_iterator iter = map_.rbegin(); iter != map_.rend(); ++iter) {
            if (sum + iter->second <= trim_size) {
                sum += iter->second;
                trim_values.push_back(iter->first);
            }
            else
                break;
        }

        size_ -= sum;

        size_t trimmed = 0;

        for (size_t i = 0; i < trim_values.size(); ++i) {
            trimmed += map_[trim_values[i]];
            map_.erase(trim_values[i]);
        }

        return trimmed;
    }

    size_t TrimLow(value_type threshold) {
        std::deque<value_type> trim_values;

        size_t sum = 0;

        for (typename std::map<value_type, size_t>::iterator iter = map_.begin(); iter != map_.end(); ++iter) {
            if (iter->first < threshold) {
                sum += iter->second;
                trim_values.push_back(iter->first);
            }
            else
                break;
        }

        size_ -= sum;

        size_t trimmed = 0;

        for (size_t i = 0; i < trim_values.size(); ++i) {
            trimmed += map_[trim_values[i]];
            map_.erase(trim_values[i]);
        }

        return trimmed;
    }

    void swap(Histgram<value_type> &histgram) {
        std::swap(size_, histgram.size_);
        map_.swap(histgram.map_);
    }

    bool empty() const {
        return map_.empty();
    }

    size_t size() const {
        return size_;
    }

    void clear() {
        map_.clear();
        size_ = 0;
        lock_ = 0;
    }

  private:
    std::map<value_type, size_t> map_;
    size_t size_;
    volatile int lock_;
};

namespace std {
template <typename T> inline void swap(Histgram<T> &histgram1, Histgram<T> &histgram2) {
    histgram1.swap(histgram2);
}
}

#endif

