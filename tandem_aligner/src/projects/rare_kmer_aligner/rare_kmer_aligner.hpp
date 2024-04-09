//
// Created by Andrey Bzikadze on 2/7/21.
//

#pragma once

#include <optional>
#include <set>
#include <unordered_set>
#include <unordered_map>

#include <common/logging.hpp>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>

#include "rolling_hash/rolling_hash.hpp"

namespace rare_kmer_aligner {

    namespace details {

        /**
         * @brief 获取稀有 k-mer
         *
         * 根据给定的 Contig 对象、k 值、最大稀有频率和哈希器对象，获取稀有 k-mer 的集合。
         *
         * @param contig Contig 对象
         * @param k k 值
         * @param max_rare_freq 最大稀有频率
         * @param hasher 哈希器对象
         *
         * @return 稀有 k-mer 的集合
         */
        template <typename htype>
        std::set<htype>
        // std::set<std::string>
        get_rare_kmers(const Contig & contig, const size_t k, const size_t max_rare_freq,
                       const RollingHash<htype> & hasher) {
            std::unordered_map<htype, size_t> counter;
            // KWH<htype> kwh(hasher, contig.seq, 0);
            // while(true) {
            //     counter[kwh.hash()] += 1;
            //     if (!kwh.hasNext()) {
            //         break;
            //     }
            //     kwh = kwh.next();
            // }
            for (size_t i = 0; i < contig.seq.size() - k + 1; ++i) {
                const auto kmer = contig.seq.Subseq(i, i+k);
                counter[hasher.hash(kmer, 0)] += 1;
            }

            std::set<htype> rare_kmers;
            for (const auto & [hash, cnt] : counter) {
                if (cnt <= max_rare_freq) {
                    rare_kmers.insert(hash);
                }
            }
            return rare_kmers;
        }

        template <typename htype>
        struct PosHash {
            size_t pos {0};
            htype hash {0};
        };

        template<typename htype>
        std::ostream & operator << (std::ostream & os, const PosHash<htype> & pos_hash) {
            os << pos_hash.pos;
            return os;
        }


        template <typename htype>
        using PosHashVector = std::vector<PosHash<htype>>;

        /**
         * @brief 过滤 kmers
         *
         * 根据给定的稀有 kmers 集合，从给定的 Contig 中筛选出稀有的 kmers，并返回其位置和哈希值。
         *
         * @param contig Contig 对象
         * @param rare_kmers 稀有 kmers 集合
         * @param k k-mer 的长度
         * @param hasher RollingHash 对象
         *
         * @return 包含稀有 kmers 位置和哈希值的 PosHashVector 对象
         */
        template <typename htype>
        PosHashVector<htype> filter_kmers(const Contig & contig,
                                   const std::unordered_set<htype> & rare_kmers,
                                   const size_t k,
                                   const RollingHash<htype> & hasher) {
            PosHashVector<htype> v;
            for (size_t i = 0; i < contig.seq.size() - k + 1; ++i) {
                const auto kmer = contig.seq.Subseq(i, i+k);
                const htype hash = hasher.hash(kmer, 0);
                if (rare_kmers.contains(hash)) {
                    v.push_back({i, hash});
                }
            }
            return v;
        }

        template <typename htype>
        struct AlignedPosHash {
            std::optional<PosHash<htype>> left;
            std::optional<PosHash<htype>> right;

            bool is_match() const {
                return left.has_value() and right.has_value();
            }

            bool is_del() const {
                return not is_match() and right.has_value();
            }

            /**
             * 检查是否 right 不包含值（即 right 为空）而 left 包含值。
             * 这种情况通常表示一个“插入”操作，即 left 指向的位置在比较序列中存在，但在原始序列中不存在。
            */
            bool is_ins() const {
                return not is_match() and left.has_value();
            }
        };

        template <typename htype>
        using KmerAlignment = std::vector<AlignedPosHash<htype>>; // kmer对齐表示的数据结构

        template<typename htype>
        KmerAlignment<htype> glob_align(const PosHashVector<htype> & ph1, const PosHashVector<htype> & ph2,
                                        logging::Logger & logger,
                                        const int64_t match = 1,
                                        const int64_t mismatch = -1000000,
                                        const int64_t indel = 0) {
            size_t n = ph1.size() + 1;
            size_t m = ph2.size() + 1;
            std::vector<std::vector<int64_t>> score;
            for (size_t i = 0; i < n; ++i) {
                score.emplace_back(m, 0);
            }

            for (size_t i = 1; i < n; ++i) {
                score[i][0] = score[i-1][0] + indel;
            }

            for (size_t j = 1; j < m; ++j) {
                score[0][j] = score[0][j-1] + indel;
            }

            for (size_t i = 1; i < n; ++i) {
                if (i % 1000 == 0) {
                    logger.info() << "Alignment iteration " << i << " out of " << n << std::endl;
                }
                for (size_t j = 1; j < m; ++j) {
                    const htype & h1 = ph1[i-1].hash;
                    const htype & h2 = ph2[j-1].hash;
                    int64_t add = h1 == h2 ? match : mismatch;
                    score[i][j] = std::max({score[i-1][j-1] + add,
                                            score[i-1][j] + indel,
                                            score[i][j-1] + indel});
                }
            }
            logger.info() << "Alignment score: " << score[n-1][m-1] << std::endl;

            KmerAlignment<htype> alignment;
            size_t i = n-1;
            size_t j = m-1;
            while (i > 0 and j > 0) {
                const htype & h1 = ph1[i-1].hash;
                const htype & h2 = ph2[j-1].hash;
                int64_t add = h1 == h2 ? match : mismatch;
                if (score[i][j] == score[i-1][j-1] + add) {
                    alignment.push_back({ph1[i-1], ph2[j-1]});
                    --i;
                    --j;
                } else if (score[i][j] == score[i-1][j] + indel) {
                    alignment.push_back({ph1[i-1], std::nullopt});
                    --i;
                } else {
                    VERIFY(score[i][j] == score[i][j-1] + indel);
                    alignment.push_back({std::nullopt, ph2[j-1]});
                    --j;
                }
            }
            while (j > 0) {
                alignment.push_back({std::nullopt, ph2[j-1]});
                --j;
            }
            while (i > 0) {
                alignment.push_back({ph1[i-1], std::nullopt});
                --i;
            }

            std::reverse(alignment.begin(), alignment.end());
            return alignment;
        }

        /**
         * @brief 输出比对结果
         *
         * 将给定的 Kmer 比对结果输出到指定的输出流中。
         *
         * @param alignment Kmer 比对结果
         * @param os 输出流
         */
        template<typename htype>
        void output_alignment(const details::KmerAlignment<htype> & alignment, std::ostream &os) {
            for (const auto & [a, b] : alignment) {
                if (a and b) {
                    VERIFY(a->hash == b->hash);
                    os << *a << '\t' << *b << '\n';
                } else if (a) {
                    os << *a << "\t - \n";
                } else {
                    VERIFY(b);
                    os << "- \t" << *b << "\n";
                }
            }
        }

        /**
         * 该结构体可能用于生物信息学中的序列比对，
         * 其中st1和en1表示第一个序列的对齐范围，st2和en2表示第二个序列的对齐范围，
         * 而is_match表示这两个范围是否匹配。get_size()函数则用于获取对齐的长度
         * 
        */
        struct AlignmentBlock {
            size_t st1 {std::numeric_limits<size_t>::max()}; /* 初始化为size_t的最大值 */
            size_t st2 {std::numeric_limits<size_t>::max()};
            size_t en1 {std::numeric_limits<size_t>::max()};
            size_t en2 {std::numeric_limits<size_t>::max()};
            bool is_match {false};

            size_t get_size() const {
                return std::max(en1 - st1, en2 - st2);
            }
        };

        std::ostream & operator << (std::ostream & os, const AlignmentBlock & block) {
            os << block.st1 << "\t" << block.en1 << "\t" << block.st2 << "\t" << block.en2 << "\n";
            return os;
        }

        using AlignmentBlocks = std::vector<AlignmentBlock>;

        /**
         * @brief 压缩对齐块
         *
         * 将给定的 Kmer 对齐对象压缩成对齐块，并返回结果。
         *
         * @param alignment Kmer 对齐对象
         *
         * @return 对齐块
         */
        template<typename htype>
        AlignmentBlocks compress_alignment(const KmerAlignment<htype> & alignment) {
            AlignmentBlocks blocks;
            auto it = alignment.cbegin();
            while (it != alignment.end()) {
                bool is_match = it->is_match();
                size_t st1 = it->left. has_value() ? it->left. value().pos : std::numeric_limits<size_t>::max();
                size_t st2 = it->right.has_value() ? it->right.value().pos : std::numeric_limits<size_t>::max();
                size_t en1 {st1};
                size_t en2 {st2};
                while (it != alignment.end()) {
                    const AlignedPosHash<htype> pos_hash { *it };
                    if (pos_hash.is_match() != is_match)
                        break;

                    if (is_match) {
                        en1 = pos_hash.left.value().pos;
                        en2 = pos_hash.right.value().pos;
                    } else {
                        if (pos_hash.is_del()) {
                            en2 = pos_hash.right.value().pos;
                            st2 = std::min(st2, en2);
                        } else {
                            VERIFY(pos_hash.is_ins());
                            en1 = pos_hash.left.value().pos;
                            st1 = std::min(en1, st1);
                        }
                    }
                    ++it;
                }
                blocks.push_back({st1, st2, en1, en2, is_match});
            }
            if (not blocks.empty() and blocks.back().is_match) {
                blocks.emplace_back();
            }
            if (not blocks.empty() and not blocks.front().is_match) {
                blocks.erase(blocks.begin());
            }
            return blocks;
        }

        /**
         * @brief 将 Kmer 比对结果转换为块
         *
         * 将给定的 Kmer 比对结果转换为块形式，并返回转换后的块。
         *
         * @param alignment Kmer 比对结果
         * @param tol_gap 允许的最大间隔大小
         * @param k Kmer 的长度
         *
         * @return 转换后的块
         */
        template<typename htype>
        AlignmentBlocks alignment2blocks(const KmerAlignment<htype> & alignment, size_t tol_gap, size_t k) {
            const AlignmentBlocks blocks = compress_alignment(alignment);

            AlignmentBlocks merged_blocks;

            if (blocks.empty()) {
                return merged_blocks;
            }

            size_t st1 = blocks.begin()->st1;
            size_t st2 = blocks.begin()->st2;
            size_t en1 {0};
            size_t en2 {0};
            bool last_add { false };
            for (auto it0 = blocks.cbegin(); it0 != blocks.cend(); std::advance(it0, 2)) {
                VERIFY(it0->is_match);
                auto it = it0;
                ++it;
                if (last_add) {
                    st1 = it0->st1;
                    st2 = it0->st2;
                }
                last_add = false;
                en1 = it0->en1;
                en2 = it0->en2;
                if (it->get_size() > tol_gap) {
                    merged_blocks.push_back({st1, st2, en1, en2, true});
                    last_add = true;
                }
            }
            if (not last_add) {
                merged_blocks.push_back({st1, st2, en1, en2, true});
            }

            for (AlignmentBlock & block : merged_blocks) {
                block.en1 += k;
                block.en2 += k;
            }
            return merged_blocks;
        }
    }

    /**
     * @brief 对两个序列文件进行稀有kmer比对
     *
     * 使用给定的kmer大小、最大稀有频率、容忍间隙和输出目录，对两个序列文件进行稀有kmer比对。
     *
     * @param first_path 第一个序列文件的路径
     * @param second_path 第二个序列文件的路径
     * @param k kmer的大小
     * @param max_rare_freq 最大稀有频率
     * @param tol_gap 容忍间隙
     * @param outdir 输出目录
     * @param logger 日志记录器
     * @param base 哈希函数的基数，默认为239
     */
    template <typename htype>
    void rare_kmer_align(const std::experimental::filesystem::path & first_path,
                         const std::experimental::filesystem::path & second_path,
                         const size_t k, const size_t max_rare_freq, const size_t tol_gap,
                         const std::experimental::filesystem::path & outdir,
                         logging::Logger & logger,
                         const size_t base = 239) {
        io::SeqReader first_reader(first_path); // 读入第一条序列
        std::vector<Contig> first_vec { first_reader.readAllContigs() }; // 序列转换为 Contig 对象数组
        VERIFY(first_vec.size() == 1);
        Contig first { std::move(first_vec[0]) } ;
        logger.info() << "First length " << first.seq.size() << ", name " << first.id
                      << " from " << first_path << std::endl;

        io::SeqReader second_reader(second_path); // 第二条序列
        std::vector<Contig> second_vec { second_reader.readAllContigs() };
        VERIFY(second_vec.size() == 1);
        Contig second { std::move(second_vec[0]) } ;
        logger.info() << "Second length " << second.seq.size() << ", name " << second.id
                      << " from " << second_path << std::endl;

        const RollingHash<htype> hasher(k, base); // rolling hash 对象

        std::set<htype> rare_kmers_1 { details::get_rare_kmers(first, k, max_rare_freq, hasher) };
        logger.info() << "# rare kmers in first string " << rare_kmers_1.size() << "\n";

        std::set<htype> rare_kmers_2 { details::get_rare_kmers(second, k, max_rare_freq, hasher) };
        logger.info() << "# rare kmers in second string " << rare_kmers_2.size() << "\n";

        /* 计算rare_kmers_1和rare_kmers_2的交集,并将结果存储在一个新的unordered_set中 */
        std::vector<htype> rare_kmers_vec;
        std::set_intersection(rare_kmers_1.begin(), rare_kmers_1.end(),
                              rare_kmers_2.begin(), rare_kmers_2.end(),
                              std::back_inserter(rare_kmers_vec));
        std::unordered_set<htype> rare_kmers(rare_kmers_vec.begin(), rare_kmers_vec.end());
        logger.info() << "# rare kmers in both strings " << rare_kmers.size() << "\n";

        details::PosHashVector<htype> pos_hash_vec_1 = details::filter_kmers(first, rare_kmers, k, hasher);
        logger.info() << "|pos hash vec 1| " << pos_hash_vec_1.size() << "\n";

        details::PosHashVector<htype> pos_hash_vec_2 = details::filter_kmers(second, rare_kmers, k, hasher);
        logger.info() << "|pos hash vec 2| " << pos_hash_vec_2.size() << "\n";

        logger.info() << "Alignment computation started\n";
        const details::KmerAlignment<htype> alignment = details::glob_align(pos_hash_vec_1, pos_hash_vec_2, logger);
        logger.info() << "Alignment computation finished\n";

        std::ofstream aligner_output;
        const std::experimental::filesystem::path aligner_fn { outdir / "glob_alignment.tsv" };
        aligner_output.open(aligner_fn);
        details::output_alignment<htype>(alignment, aligner_output);
        aligner_output.close();
        logger.info() << "Alignment exported to " << aligner_fn << "\n";

        const details::AlignmentBlocks blocks = details::alignment2blocks(alignment, tol_gap, k);
        std::ofstream blocks_output;
        const std::experimental::filesystem::path blocks_fn { outdir / "blocks.tsv" };
        blocks_output.open(blocks_fn);
        for (const auto & block : blocks) {
            blocks_output << block;
        }
        blocks_output.close();
        logger.info() << "Matching blocks exported to " << blocks_fn << "\n";
    }

} // End namespace rare_kmer_aligner

