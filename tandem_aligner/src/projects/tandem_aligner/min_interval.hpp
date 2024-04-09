//
// Created by Andrey Bzikadze on 05/18/22.
//

#pragma once

#include <utility>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <experimental/iterator>
#include <experimental/filesystem>
#include "lcp_interval.hpp"
#include "common/dir_utils.hpp"

namespace tandem_aligner {

/**
 * @brief 具有最小长度的区间以及该区间再两条序列中的起始坐标 
 *  
 */
class MinInterval {
    int len{0};
    std::vector<int> fst_coords;
    std::vector<int> snd_coords;

 public:
    MinInterval() = default;
    explicit MinInterval(const int len) : len{len} {}

    void PushBackFst(const int start) { fst_coords.push_back(start); }
    void PushBackSnd(const int start) { snd_coords.push_back(start); }

    [[nodiscard]] int GetLen() const { return len; }
    [[nodiscard]] const std::vector<int> &GetFstCoords() const { return fst_coords; }
    [[nodiscard]] const std::vector<int> &GetSndCoords() const { return snd_coords; }
};

std::ostream &operator<<(std::ostream &os, const MinInterval &interval);

/**
 * @brief 最大不相交区间集合
 *
 * 该类表示一个最大不相交区间集合，其中每个区间都有一个频率，频率可以是fst_freq或snd_freq。
 *
 */
class MaxDisjointIntervalCollection {
    int fst_freq{1}, snd_freq{1}; // 分别表示第一序列和第二序列的频率。
    std::unordered_map<int, MinInterval> intervals; // 键是区间的分类标识，值是对应的 MinInterval 对象

 public:
    MaxDisjointIntervalCollection(const int fst_freq, const int snd_freq) :
        fst_freq{fst_freq}, snd_freq{snd_freq} {}

    /* 重载运算符，用于通过分类标识访问或添加 MinInterval 对象 */
    MinInterval &operator[](const int &clas) { return intervals[clas]; }

    decltype(intervals)::iterator begin() { return intervals.begin(); }
    decltype(intervals)::iterator end() { return intervals.end(); }
    [[nodiscard]] decltype(intervals)::const_iterator begin() const {
        return intervals.begin();
    }
    [[nodiscard]] decltype(intervals)::const_iterator end() const {
        return intervals.end();
    }
    [[nodiscard]] decltype(intervals)::const_iterator cbegin() const {
        return intervals.cbegin();
    }
    [[nodiscard]] decltype(intervals)::const_iterator cend() const {
        return intervals.cend();
    }

    template<class... Args>
    std::pair<decltype(intervals)::iterator, bool> Emplace(Args &&... args) {
        return intervals.emplace(args...);// emplace 方法在容器内部直接构造新元素
    }

    int GetFstFreq() const { return fst_freq; }
    int GetSndFreq() const { return snd_freq; }
};

std::ostream &operator<<(std::ostream &os, const MaxDisjointIntervalCollection &col);

using MaxDisjointIntervalCollections = std::vector<MaxDisjointIntervalCollection>;

std::ostream &operator<<(std::ostream &os, const MaxDisjointIntervalCollections &cols);


/**
 * @brief 最大不相交区间查找器
 * 
 * 
*/
class MaxDisjointIntervalFinder {
    const int max_freq{1};
    const int min_freq{1};
    bool force_highfreq_search{false};
    bool exprt{true};
    std::experimental::filesystem::path outdir;

    /// @brief 表示第一个(第二个)seq中的fst_freq(snd_freq)次数的最短和最长前缀
    struct MinMaxRarePrefixArray {
        // represents shortest and longest prefixes present fst_freq (snd_freq) times in the first (second) seq
        struct MinMaxRarePrefix {
            int clas{-1};
            int min_len{std::numeric_limits<int>::max()}; // length of shortest (fst_freq,snd_freq)-prefix
            int max_len{std::numeric_limits<int>::max()}; // length of longest (fst_freq,snd_freq)-prefix

            [[nodiscard]] bool IsInit() const {
                return clas!=-1 and min_len!=std::numeric_limits<int>::max();
            }
        };

        const int fst_freq{0}, snd_freq{0};
        std::vector<MinMaxRarePrefix> mrp; // shortest rare prefix 最短稀有前缀

        MinMaxRarePrefixArray(const int fst_freq, const int snd_freq,
                      const int length) :
            fst_freq{fst_freq}, snd_freq{snd_freq}, mrp(length) {}

        // 重载下标运算符，用于访问 mrp[i]
        MinMaxRarePrefix &operator[](int i) { return mrp[i]; }
        const MinMaxRarePrefix &operator[](int i) const { return mrp[i]; }

        [[nodiscard]] int Size() const { return mrp.size(); }
    };

    [[nodiscard]] std::vector<MinMaxRarePrefixArray> GetMRP(const LCP &lcp,
                                                    const int fst_len) const;

    struct ClassCoord {
        int lcp_class{-1};
        int coord{-1};
    };

    [[nodiscard]] MaxDisjointIntervalCollections
    PrefixesToIntervals(const std::vector<MinMaxRarePrefixArray> &mrp_vec,
                        const int fst_len) const;

    void OutputMRP(const std::vector<MinMaxRarePrefixArray> &srp_vec,
                   std::experimental::filesystem::path outpath) const;

 public:
    MaxDisjointIntervalFinder(int max_freq,
                      bool force_highfreq_search,
                      bool exprt,
                      std::experimental::filesystem::path outdir) :
                      max_freq{max_freq}, exprt{exprt},
                      force_highfreq_search{force_highfreq_search},
                      outdir{std::move(outdir)} {
        ensure_dir_existance(this->outdir);
    }

    [[nodiscard]] MaxDisjointIntervalCollections Find(const suffix_array::LCP<std::string> &lcp,
                                              const int fst_len) const;
};

} // namespace MinInterval
