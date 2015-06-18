/**
 * @file vertex_status.h
 * @brief VertexStatus Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-05
 */

#ifndef __GRAPH_VERTEX_STATUS_H_

#define __GRAPH_VERTEX_STATUS_H_

#include <stdint.h>
#include <algorithm>

/**
 * @brief It is a class for storing the status of a vertex. It provides many 
 * useful functions to access the status of a vertex.
 */
class VertexStatus
{
public:
    VertexStatus(): status_(0) {}
    VertexStatus(const VertexStatus &vertex_status): status_(vertex_status.status_) {}

    const VertexStatus &operator =(const VertexStatus &vertex_status)
    { status_ = vertex_status.status_; return *this; }

    void SetUsedFlag() { SetFlag(kVertexStatusFlagUsed); }
    void ResetUsedFlag() { ResetFlag(kVertexStatusFlagUsed); }
    bool IsUsed() const { return GetFlag(kVertexStatusFlagUsed); }

    void SetDeadFlag() { SetFlag(kVertexStatusFlagDead); }
    void ResetDeadFlag() { ResetFlag(kVertexStatusFlagDead); }
    bool IsDead() const { return GetFlag(kVertexStatusFlagDead); }

    int GetLockID()
    {
        if (status_ & kVertexStatusFlagLock)
            return status_ & kVertexStatusMaskLock;
        return -1;
    }

    bool Lock(int id)
    {
        uint16_t old_status = status_;
        if (old_status & kVertexStatusFlagLock)
            return false;

        status_ = (old_status & ~kVertexStatusMaskLock) | kVertexStatusFlagLock | id;
        return true;
    }

    bool LockPreempt(int id)
    {
        uint16_t old_status = status_;
        int old_id = -1;
        if (old_status & kVertexStatusFlagLock)
            old_id = old_status & kVertexStatusMaskLock;

        if (old_id >= id)
            return false;

        status_ = (old_status & ~kVertexStatusMaskLock) | kVertexStatusFlagLock | id;
        return true;
    }

    void swap(VertexStatus &vertex_status)
    { if (this != &vertex_status) std::swap(status_, vertex_status.status_); }

    void clear() { status_ = 0; }

    static const uint16_t kVertexStatusFlagDead = 0x8000U;
    static const uint16_t kVertexStatusFlagUsed = 0x4000U;
    static const uint16_t kVertexStatusFlagLock = 0x2000U;
    static const uint16_t kVertexStatusMaskLock = 0x1FFFU;

private:
    bool GetFlag(uint16_t flag) const { return status_ & flag; }
    void SetFlag(uint16_t flag) { status_ |= flag; }
    void ResetFlag(uint16_t flag) { status_ &= ~flag; }

    uint16_t status_;
};

namespace std
{
template <> inline void swap(VertexStatus &x, VertexStatus &y) { x.swap(y); }
}

#endif

