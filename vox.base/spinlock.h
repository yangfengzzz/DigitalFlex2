//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <atomic>

namespace vox {
class SpinLock {
public:
    void lock() {
        while (m_lock.test_and_set(std::memory_order_acquire)) {
        }
    }

    void unlock() { m_lock.clear(std::memory_order_release); }

    SpinLock() = default;
    SpinLock(SpinLock const &other){};
    SpinLock &operator=(SpinLock const &other) { return *this; }

private:
    std::atomic_flag m_lock = ATOMIC_FLAG_INIT;
};
}  // namespace vox
