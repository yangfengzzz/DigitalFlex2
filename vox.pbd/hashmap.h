//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <cstdlib>
#include <map>

#include "vox.base/common.h"

namespace vox::utility {
template <class KeyType>
inline unsigned int hashFunction(const KeyType &key) {
    return 0u;
}

template <class KeyType, class ValueType>
class Hashmap {
public:
    typedef typename std::map<unsigned int, ValueType> KeyValueMap;

private:
    KeyValueMap **m_hashMap;
    unsigned int m_bucketCount;
    unsigned int m_moduloValue;

protected:
    FORCE_INLINE void init() {
        m_hashMap = new KeyValueMap *[m_bucketCount];
        for (unsigned int i = 0; i < m_bucketCount; i++) {
            m_hashMap[i] = nullptr;
        }
    }

    FORCE_INLINE void cleanup() {
        if (m_hashMap) {
            for (unsigned int i = 0; i < m_bucketCount; i++) {
                if (m_hashMap[i] != nullptr) {
                    m_hashMap[i]->clear();
                    delete m_hashMap[i];
                }
            }
            delete[] m_hashMap;
            m_hashMap = nullptr;
        }
    }

public:
    FORCE_INLINE explicit Hashmap(const unsigned int bucketCount) {
        // Use a bucket count of 2^n => faster modulo
        unsigned int val = bucketCount;
        unsigned int powerOfTwo = 1u;
        while (powerOfTwo < val) powerOfTwo <<= 1;
        m_bucketCount = powerOfTwo;
        m_moduloValue = m_bucketCount - 1u;
        init();
    }

    ~Hashmap() { cleanup(); }

    FORCE_INLINE void clear() {
        cleanup();
        init();
    }

    FORCE_INLINE KeyValueMap *getKeyValueMap(const unsigned int index) { return m_hashMap[index]; }

    FORCE_INLINE void reset() {
        for (unsigned int i = 0; i < m_bucketCount; i++) {
            if (m_hashMap[i] != nullptr) {
                m_hashMap[i]->clear();
            }
        }
    }

    /** Return the bucket count.
     */
    [[nodiscard]] FORCE_INLINE unsigned int bucket_count() const { return m_bucketCount; }

    /** Find element.
     */
    FORCE_INLINE ValueType *find(const KeyType &key) {
        const unsigned int hashValue = hashFunction<KeyType>(key);
        const unsigned int mapIndex = hashValue & m_moduloValue;
        if (m_hashMap[mapIndex] != nullptr) {
            typename KeyValueMap::iterator &iter = (*m_hashMap[mapIndex]).find(hashValue);
            if (iter != (*m_hashMap[mapIndex]).end()) return &iter->second;
        }
        return nullptr;
    }

    /** Insert element.
     */
    FORCE_INLINE void insert(const KeyType &key, const ValueType &value) {
        const unsigned int hashValue = hashFunction<KeyType>(key);
        const unsigned int mapIndex = hashValue & m_moduloValue;
        if (m_hashMap[mapIndex] == nullptr) {
            m_hashMap[mapIndex] = new KeyValueMap();
        }
        (*m_hashMap[mapIndex])[hashValue] = value;
    }

    /** Remove the given element and return true, if the element was found.
     */
    FORCE_INLINE void remove(const KeyType &key) {
        const unsigned int hashValue = hashFunction<KeyType>(key);
        const unsigned int mapIndex = hashValue & m_moduloValue;
        if (m_hashMap[mapIndex] != nullptr) {
            m_hashMap[mapIndex]->erase(hashValue);
            if (m_hashMap[mapIndex]->size() == 0) {
                delete m_hashMap[mapIndex];
                m_hashMap[mapIndex] = nullptr;
            }
        }
    }

    FORCE_INLINE ValueType &operator[](const KeyType &key) {
        const int hashValue = hashFunction<KeyType>(key);
        const unsigned int mapIndex = hashValue & m_moduloValue;
        if (m_hashMap[mapIndex] == nullptr) {
            m_hashMap[mapIndex] = new KeyValueMap();
        }
        return (*m_hashMap[mapIndex])[hashValue];
    }

    FORCE_INLINE const ValueType &operator[](const KeyType &key) const {
        const unsigned int hashValue = hashFunction<KeyType>(key, m_bucketCount);
        const unsigned int mapIndex = hashValue & m_moduloValue;
        if (m_hashMap[mapIndex] == nullptr) {
            m_hashMap[mapIndex] = new KeyValueMap();
        }
        return (*m_hashMap[mapIndex])[hashValue];
    }

    FORCE_INLINE const ValueType *query(const KeyType &key) const {
        const unsigned int hashValue = hashFunction<KeyType>(key);
        const unsigned int mapIndex = hashValue & m_moduloValue;
        if (m_hashMap[mapIndex] == nullptr) {
            return nullptr;
        }
        typename KeyValueMap::iterator &it = m_hashMap[mapIndex]->find(hashValue);
        if (it != m_hashMap[mapIndex]->end()) return &it->second;
        return nullptr;
    }

    FORCE_INLINE ValueType *query(const KeyType &key) {
        const unsigned int hashValue = hashFunction<KeyType>(key);
        const unsigned int mapIndex = hashValue & m_moduloValue;
        if (m_hashMap[mapIndex] == nullptr) {
            return nullptr;
        }
        const typename KeyValueMap::iterator &it = m_hashMap[mapIndex]->find(hashValue);
        if (it != m_hashMap[mapIndex]->end()) return &it->second;
        return nullptr;
    }
};
}  // namespace vox::utility
