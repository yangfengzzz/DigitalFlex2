//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <atomic>
#include <vector>

#include "vox.base/config.h"

namespace vox {
//MARK: - PointID
struct PointID {
    unsigned int point_set_id;
    unsigned int point_id;

    bool operator==(PointID const &other) const {
        return point_id == other.point_id && point_set_id == other.point_set_id;
    }
};

//MARK: - HashKey
struct HashKey {
    HashKey() = default;
    HashKey(int i, int j, int k) { this->k[0] = i, this->k[1] = j, this->k[2] = k; }

    HashKey &operator=(HashKey const &other) {
        k[0] = other.k[0];
        k[1] = other.k[1];
        k[2] = other.k[2];
        return *this;
    }

    bool operator==(HashKey const &other) const {
        return k[0] == other.k[0] && k[1] == other.k[1] && k[2] == other.k[2];
    }

    bool operator!=(HashKey const &other) const { return !(*this == other); }

    int k[3]{};
};

//MARK: - HashEntry
struct HashEntry {
    HashEntry() : n_searching_points(0u) { indices.reserve(INITIAL_NUMBER_OF_INDICES); }

    HashEntry(PointID const &id) : n_searching_points(0u) { add(id); }

    void add(PointID const &id) { indices.push_back(id); }

    void erase(PointID const &id) {
        auto it = std::find(indices.begin(), indices.end(), id);
        if (it != indices.end()) indices.erase(it);
    }

    [[nodiscard]] unsigned int n_indices() const { return static_cast<unsigned int>(indices.size()); }

    std::vector<PointID> indices;
    unsigned int n_searching_points;
};

struct SpatialHasher {
    std::size_t operator()(HashKey const &k) const {
        return static_cast<size_t>(static_cast<int64_t>(73856093) * static_cast<int64_t>(k.k[0]) ^
                                   static_cast<int64_t>(19349663) * static_cast<int64_t>(k.k[1]) ^
                                   static_cast<int64_t>(83492791) * static_cast<int64_t>(k.k[2]));
    }
};

//MARK: - ActivationTable
class ActivationTable {
private:
    std::vector<std::vector<unsigned char>> m_table;

public:
    bool operator==(ActivationTable const &other) const { return m_table == other.m_table; }

    bool operator!=(ActivationTable const &other) const { return !(m_table == other.m_table); }

    /** Add point set. If search_neighbors is true, neighbors in all other point
     * sets are searched. If find_neighbors is true, the new point set is
     * activated in the neighborhood search of all other point sets.
     */
    void add_point_set(bool search_neighbors = true, bool find_neighbors = true) {
        // add column to each row
        auto size = m_table.size();
        for (auto i = 0u; i < size; i++) {
            m_table[i].resize(size + 1);
            m_table[i][size] = static_cast<unsigned char>(find_neighbors);
        }

        // add new row
        m_table.resize(size + 1);
        m_table[size].resize(size + 1);
        for (auto i = 0u; i < size + 1; i++) m_table[size][i] = static_cast<unsigned char>(search_neighbors);
    }

    /** Activate/Deactivate that neighbors in point set index2 are found when
     * searching for neighbors of point set index1.
     */
    void set_active(unsigned int index1, unsigned int index2, bool active) {
        m_table[index1][index2] = static_cast<unsigned char>(active);
    }

    /** Activate/Deactivate all point set pairs containing the given index. If
     * search_neighbors is true, neighbors in all other point sets are searched.
     * If find_neighbors is true, the new point set is activated in the
     * neighborhood search of all other point sets.
     */
    void set_active(unsigned int index, bool search_neighbors = true, bool find_neighbors = true) {
        auto size = m_table.size();
        for (auto i = 0u; i < size; i++) {
            m_table[i][index] = static_cast<unsigned char>(find_neighbors);
            m_table[index][i] = static_cast<unsigned char>(search_neighbors);
        }
        m_table[index][index] = static_cast<unsigned char>(search_neighbors && find_neighbors);
    }

    /** Activate/Deactivate all point set pairs.
     */
    void set_active(bool active) {
        auto size = m_table.size();
        for (auto i = 0u; i < size; i++)
            for (auto j = 0u; j < size; j++) m_table[i][j] = static_cast<unsigned char>(active);
    }

    [[nodiscard]] bool is_active(unsigned int index1, unsigned int index2) const {
        return m_table[index1][index2] != 0;
    }

    [[nodiscard]] bool is_searching_neighbors(unsigned int const index) const {
        for (unsigned char i : m_table[index]) {
            if (i) {
                return true;
            }
        }
        return false;
    }
};
}  // namespace vox
