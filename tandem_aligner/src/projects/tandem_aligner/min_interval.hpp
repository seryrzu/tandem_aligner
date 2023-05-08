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
#include "common/logging.hpp"

namespace tandem_aligner {

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

class MinIntervalCollection {
    int fst_freq{1}, snd_freq{1};
    std::unordered_map<int, MinInterval> intervals;

 public:
    MinIntervalCollection(const int fst_freq, const int snd_freq) :
        fst_freq{fst_freq}, snd_freq{snd_freq} {}

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
        return intervals.emplace(args...);
    }

    int GetFstFreq() const { return fst_freq; }
    int GetSndFreq() const { return snd_freq; }
};

std::ostream &operator<<(std::ostream &os, const MinIntervalCollection &col);

using MinIntervalCollections = std::vector<MinIntervalCollection>;

std::ostream &operator<<(std::ostream &os, const MinIntervalCollections &cols);

class MinIntervalFinder {
    const int max_freq{1};
    const int min_freq{1};
    bool force_highfreq_search{false};
    bool exprt{true};
    std::experimental::filesystem::path outdir;
    bool max_unique{false};
    logging::Logger &logger;

    struct MinRarePrefix {
        struct ClassMinLen {
            int clas{-1};
            int len{std::numeric_limits<int>::max()};
            int maxlen{std::numeric_limits<int>::max()};

            [[nodiscard]] bool IsInit() const {
                return clas!=-1 and len!=std::numeric_limits<int>::max();
            }
        };

        const int fst_freq{0}, snd_freq{0};
        std::vector<ClassMinLen> mrp; // shortest rare prefix

        MinRarePrefix(const int fst_freq, const int snd_freq,
                      const int length) :
            fst_freq{fst_freq}, snd_freq{snd_freq}, mrp(length) {}

        ClassMinLen &operator[](int i) { return mrp[i]; }
        const ClassMinLen &operator[](int i) const { return mrp[i]; }

        [[nodiscard]] int Size() const { return mrp.size(); }
    };

    [[nodiscard]] std::vector<MinRarePrefix> GetMRP(const LCP &lcp,
                                                    const int fst_len) const;

    struct ClassCoord {
        int lcp_class{-1};
        int coord{-1};
    };

    [[nodiscard]] MinIntervalCollections
    PrefixesToIntervals(const std::vector<MinRarePrefix> &mrp_vec,
                        const int fst_len) const;

    void OutputMRP(const std::vector<MinRarePrefix> &srp_vec,
                   std::experimental::filesystem::path outpath) const;

 public:
    MinIntervalFinder(int max_freq,
                      bool force_highfreq_search,
                      bool exprt,
                      std::experimental::filesystem::path outdir,
                      bool max_unique,
                      logging::Logger &logger) :
                      max_freq{max_freq}, exprt{exprt},
                      force_highfreq_search{force_highfreq_search},
                      outdir{std::move(outdir)},
                      max_unique{max_unique},
                      logger{logger} {
        ensure_dir_existance(this->outdir);
    }

    [[nodiscard]] MinIntervalCollections Find(const suffix_array::LCP<std::string> &lcp,
                                              const int fst_len) const;
};

} // namespace MinInterval
