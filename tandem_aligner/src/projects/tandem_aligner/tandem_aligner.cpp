//
// Created by Andrey Bzikadze on 05/05/22.
//
/*
#当前文件 #当前文件 #tandem_aligner #py #common #utils 
#rare_kmer_aligner #sequences #suffix_array #common #rolling_hash #min_queue #utils
*/
#include <iomanip>
#include <iostream>

#include <common/cl_parser.hpp>
#include <common/logging.hpp>

#include "tandem_aligner.hpp"

int main(int argc, char **argv) {
    CLParser parser{{"output-dir=", "first=", "second=", "max_count=50",
                     "debug", "force_highfreq_search", "no_paths", "bridges"}, {},
                    {"o=output-dir", "f=first", "s=second", "f=max_count"}};
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cerr << "Incorrect parameters" << std::endl;
        std::cerr << parser.check() << std::endl;
        return 1;
    }

    const std::experimental::filesystem::path
        output_dir{parser.getValue("output-dir")};
    ensure_dir_existance(output_dir);

    bool debug = parser.getCheck("debug");
    logging::LoggerStorage ls{output_dir, "tandem_aligner"};
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(),
                      debug ? logging::debug : logging::trace);

    auto time_point{std::chrono::system_clock::now()};
    std::time_t now = std::chrono::system_clock::to_time_t(time_point);
    logger << "Time: " << std::put_time(std::localtime(&now), "%c %Z")
           << std::endl;

    for (size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    logging::logGit(logger, output_dir/"version.txt");

    const std::experimental::filesystem::path first_path =
        std::experimental::filesystem::canonical(parser.getValue("first"));
    const std::experimental::filesystem::path second_path =
        std::experimental::filesystem::canonical(parser.getValue("second"));

    int max_freq = std::stoi(parser.getValue("max_count"));

    bool force_highfreq_search = parser.getCheck("force_highfreq_search");
    bool no_paths = parser.getCheck("no_paths");
    bool bridges = parser.getCheck("bridges");
    tandem_aligner::TandemAligner(logger,
                                  output_dir,
                                  max_freq,
                                  force_highfreq_search,
                                  no_paths,
                                  bridges)
        .Find(first_path, second_path);

    logger.info() << "Thank you for using TandemAligner!\n";
    return 0;
}
