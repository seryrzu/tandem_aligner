//
// Created by Andrey Bzikadze on 05/05/22.
//

#pragma once

#include <queue>
#include <common/logging.hpp>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>
#include "lcp_interval.hpp"
#include "min_interval.hpp"
#include "sparse_aligner.hpp"

namespace tandem_aligner {

/**
 * @brief 最小序列任务,即一个CIGAR串
 * 其中每个CigarFragment是一种修改类型{M,X,I,D}以及对应的长度
 * @param st1 st2 分别表示第一个和第二个序列的起始位置
 * @param len1 len2 分别表示第一个和第二个序列的长度
*/
struct MinSeqTask {
    std::list<CigarFragment>::iterator cigar_it;
    int64_t st1{0}, len1{0};
    int64_t st2{0}, len2{0};
};

class TandemAligner {
    logging::Logger &logger;
    int max_freq{1};
    bool force_highfreq_search{false};
    const std::experimental::filesystem::path output_dir;
    bool no_paths;
    bool bridges;

    /**
     * @brief 读取指定路径的contig
     *
     * 使用给定的路径读取contig，并将其转换为字符串返回。
     *
     * @param path 文件路径
     *
     * @return 返回contig的字符串表示
     * @throws std::exception 如果读取过程中发生错误
     */
    [[nodiscard]] std::string ReadContig(const std::experimental::filesystem::path &path) const {
        io::SeqReader reader(path);
        std::vector<Contig> vec{reader.readAllContigs()};
        VERIFY(vec.size()==1);
        Contig contig{std::move(vec[0])};
        logger.info() << "Length " << contig.seq.size() << ", name "
                      << contig.id
                      << " from " << path << std::endl;
        return contig.str();
    }

    /**
     * @brief 拼接两个字符串，并在中间加上特定符号
     *
     * 将第一个字符串和第二个字符串拼接在一起，并在它们之间加上"$"和"#"符号。
     *
     * @param first 第一个字符串
     * @param second 第二个字符串
     *
     * @return 拼接后的字符串
     */
    [[nodiscard]] std::string ConcatContigs(const std::string &first,
                                            const std::string &second) const {
        std::stringstream concat_stream;
        concat_stream << first << '$' << second << '#';
        return concat_stream.str();
    }

    void RunTask(std::queue<MinSeqTask> &queue,
                 Cigar &main_cigar,
                 const std::string &first,
                 const std::string &second,
                 bool exprt,
                 bool assert_validity = true) const {
        MinSeqTask task = queue.front(); // 获取两个字符串的起始位置和长度
        const std::string first_substr = first.substr(task.st1, task.len1);
        const std::string second_substr = second.substr(task.st2, task.len2);
        const std::string
            concat = ConcatContigs(first_substr, second_substr); // 将两个字符串拼接在一起

        logger.debug() << "Building suffix array...\n";
        const suffix_array::SuffixArray<std::string> suf_arr(concat);
        logger.debug() << "Building LCP array...\n";
        const suffix_array::LCP<std::string> lcp(suf_arr);

        MaxDisjointIntervalFinder
            segment_finder(max_freq, force_highfreq_search, exprt,
                           output_dir/"min_interval_finder");
        logger.debug() << "Computing rare segments...\n";
        const MaxDisjointIntervalCollections
            int_col = segment_finder.Find(lcp, task.len1);

        if (exprt) {
            std::ofstream os(output_dir/"shortest_matches.tsv");
            os << int_col;
        }// 识别出最短锚点

        logger.debug() << "Aligning...\n"; // 根据主要锚点进行对齐
        Cigar cigar = SparseAligner(logger, output_dir, exprt & no_paths, exprt & bridges).Align(int_col,
                                                  first_substr, second_substr); //exprt & for primary alignment 
        auto main_cigar_it = main_cigar.AssignInterval(cigar, task.cigar_it);
        logger.debug() << "Finished alignment\n";
        /*
            代码首先执行一个序列比对，并将结果存储在cigar对象中。
            然后，它遍历cigar中的CIGAR片段，检查每对连续的片段，看它们是否都是删除或插入操作，并且长度大于1。如果是，它会将这些信息添加到队列中。
            同时，它还更新了两个子串的当前位置，基于匹配或错配操作。
        */
        if (cigar.Size() > 2) {
            int64_t i{task.st1}, j{task.st2};
            auto it1 = cigar.cbegin(), it2 = ++cigar.cbegin();
            for (; it2!=cigar.cend(); ++it1, ++it2, ++main_cigar_it) {
                const CigarFragment &fragment1 = *it1;
                const CigarFragment &fragment2 = *it2;
                if ((fragment1.mode==CigarMode::D
                    or fragment1.mode==CigarMode::I) and
                    (fragment2.mode==CigarMode::D
                        or fragment2.mode==CigarMode::I) and
                    std::min(fragment1.length, fragment2.length) > 1) {
                    if (fragment1.mode==CigarMode::D) {
                        queue.push({main_cigar_it,
                                    i, int64_t(fragment1.length),
                                    j, int64_t(fragment2.length)});
                    } else {
                        queue.push({main_cigar_it,
                                    i, int64_t(fragment2.length),
                                    j, int64_t(fragment1.length)});
                    }
                }
                if (fragment1.mode==CigarMode::M
                    or fragment1.mode==CigarMode::X) {
                    i += fragment1.length;
                    j += fragment1.length;
                    continue;
                } else if (fragment1.mode==CigarMode::D) {
                    i += fragment1.length;
                } else {
                    j += fragment1.length;
                }
            }
        }
        if (assert_validity)
            cigar.AssertValidity(first, second);
    }

    /**
     * @brief 分配不匹配的 CIGAR 记录
     *
     * 将给定的两个字符串 first 和 second 的 CIGAR 记录进行比较，
     * 将不匹配的部分重新分配并更新 CIGAR 记录。
     *
     * @param cigar CIGAR 记录对象
     * @param first 第一个字符串
     * @param second 第二个字符串
     */
    void AssignMismatches(Cigar &cigar,
                          const std::string &first,
                          const std::string &second) const {
        if (cigar.Size() < 2)
            return;

        auto it1 = cigar.begin(), it2 = ++cigar.begin();
        int i{0}, j{0};
        std::unordered_map<int, int> counter;
        for (; it2!=cigar.end(); ++it1, ++it2) {
            int64_t length1{it1->length}, length2{it2->length};
            CigarMode mode1{it1->mode}, mode2{it2->mode};
            /* 对初步比对后获得的CIGAR串判断是否存在连续的indel对 */
            if ((mode1==CigarMode::D or mode1==CigarMode::I) and
                (mode2==CigarMode::D or mode2==CigarMode::I) and
                length1==length2) {
                // 存在连续的indel对，且其长度相等
                counter[length1]++; // 记录该长度的square indel对出现的次数
                // 从CIGAR串中删除当前的indel对，并验证删除后的迭代器位置是否相同
                it2 = cigar.Erase(it2);
                it1 = cigar.Erase(it1);
                VERIFY(it1==it2);

                int64_t run_len{1}, k{1};
                bool is_eq{first[i]==second[j]};
                i++, j++;
                while (true) {
                    if (is_eq!=(first[i]==second[j]) or k==length1) {
                        it2 = cigar.Insert(it2,
                                           {run_len, is_eq ? CigarMode::M
                                                           : CigarMode::X});
                        it1 = it2++;
                        if (k==length1)
                            break;
                        is_eq = !is_eq; // 取反以跟踪下一个字符的比较结果
                        run_len = 0; //计算下一个匹配或非匹配的字符数
                    }
                    run_len++, k++, i++, j++;
                }
            } else if (mode1==CigarMode::M or mode1==CigarMode::X) {
                i += length1;
                j += length1;
            } else if (mode1==CigarMode::D) {
                i += length1;
            } else {
                j += length1;
            }
        }
        /* 识别出来长度相等的indel对，说明这对序列中可能存在匹配和错误匹配 */
        for (auto [length, cnt] : counter) {
            logger.trace() << "Square Indel-block of length " << length
                << " appears " << cnt << " times\n";
        }
        cigar.AssertValidity(first, second);
    }

 public:
    TandemAligner(logging::Logger &logger,
                  std::experimental::filesystem::path output_dir,
                  const int max_freq,
                  const bool force_highfreq_search,
                  bool no_paths,
                  bool bridges) :
        logger{logger}, output_dir{std::move(output_dir)}, max_freq{max_freq},
        force_highfreq_search{force_highfreq_search}, no_paths{no_paths}, bridges{bridges} {}

    void Find(const std::experimental::filesystem::path &first_path,
              const std::experimental::filesystem::path &second_path) const {
        io::SeqReader first_reader(first_path);
        std::vector<Contig> first_vec{first_reader.readAllContigs()};
        VERIFY(first_vec.size()==1);
        const std::string first = ReadContig(first_path);
        const std::string second = ReadContig(second_path);

        Cigar cigar;
        std::queue<MinSeqTask> queue; // 队列中存储了待处理的序列比对结果CIGAR
        bool matches_exported{false};
        queue.push({cigar.begin(), 0, (int64_t) first.size(), 0,
                    (int64_t) second.size()});
        logger.info() << "Running primary alignment...\n";
        if (no_paths) std::ofstream no_paths(output_dir/"no_paths.csv"); //just empty the file 
        if (bridges) std::ofstream no_paths(output_dir/"bridges.txt"); //just empty the file 
        RunTask(queue, cigar, first, second,
                /*export_matches*/ true);
        logger.info() << "Number of indel-blocks " << queue.size() << "\n";
        queue.pop();
        logger.info() << "Finished running primary alignment\n";
        cigar.Summary(logger);
        
        /* 对square indel 对进行迭代处理 */
        logger.info() << "Running recursive alignments...\n";

        {
            std::string cigar_outfile = output_dir/"cigar_primary.txt";
            std::ofstream cigar_os(cigar_outfile);
            cigar_os << cigar;
            logger.info() << "Primary cigar exported to " << cigar_outfile
                          << "\n\n";
        }
            /* 迭代处理的部分，对indel块重复使用RunTask*/
        for (; not queue.empty(); queue.pop()) {
            RunTask(queue, cigar, first, second,
                /*export_matches*/ false,
                /*assert_validity*/ false);
        }
        cigar.AssertValidity(first, second);
        logger.info() << "Finished running recursive alignment\n";
        cigar.Summary(logger);

        {
            std::string cigar_outfile = output_dir/"cigar_recursive.txt";
            std::ofstream cigar_os(cigar_outfile);
            cigar_os << cigar;
            logger.info() << "Cigar after recursion exported to "
                          << cigar_outfile << "\n\n";
        }

        logger.info() << "Assigning mismatches...\n";
        AssignMismatches(cigar, first, second);
        logger.info() << "Finished assigning mismatches\n";
        cigar.Summary(logger);

        {
            std::string cigar_outfile = output_dir/"cigar.txt";
            std::ofstream cigar_os(cigar_outfile);
            cigar_os << cigar;
            logger.info() << "Cigar w/ mismatches exported to " << cigar_outfile
                          << "\n";
        }
    }
};

} // namespace tandem_aligner
