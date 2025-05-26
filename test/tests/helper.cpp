#include <string>
#include <sstream>
#include <iostream>
#include <functional>
#include <fstream>
#include <filesystem>

#include <HFold.hpp>
#include "helper.hpp"

std::string capture_console_output(const std::function<void()>& func) {
    std::stringstream buffer;
    std::streambuf* oldCout = std::cout.rdbuf();  // Save cout buffer
    std::cout.rdbuf(buffer.rdbuf());              // Redirect to our buffer

    func();                                       // Run the code that prints to console

    std::cout.rdbuf(oldCout);                     // Restore cout
    return buffer.str();                          // Return captured output
}

std::string get_output(const std::string& file_path, bool delete_after_read) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open output file: " + file_path);
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();

    // Delete the file after reading
    if (delete_after_read) {
        std::error_code ec;
        std::filesystem::remove(file_path, ec);
        if (ec) {
            throw std::runtime_error("Failed to delete file: " + file_path + " - " + ec.message());
        }
    }

    return buffer.str();
}


std::vector<std::vector<Result>> hfold_test(HFoldParams& p) {
    if (!p.fileO.empty()){
        // delete the output file if it exists. (our code appends instead of overwriting, so having an existing file is problematic)
        std::error_code ec;
        if (std::filesystem::exists(p.fileO)) {
            std::filesystem::remove(p.fileO, ec);
            if (std::filesystem::exists(p.fileO)) {
                throw std::runtime_error("Failed to delete existing output file: " + p.fileO);
            }
        }
    }


    if (!p.paramFile.empty()) p.paramFile = "../../params/" +  p.paramFile;
    bool input_structure_given = !p.restricted.empty();


    std::vector<RNAEntry> inputs = get_all_inputs(p.fileI, p.sequence, p.restricted);
    std::vector<std::vector<Result>> all_results;
    for (RNAEntry& current : inputs){
        preprocess_sequence(current.sequence, current.structure, p.noConv_given);
        load_energy_parameters(p.paramFile, current.sequence);

        std::vector<Hotspot> hotspots = build_hotspots(
            current.sequence,
            current.structure,
            p.suboptCount
        );

        std::vector<Result>  results  = fold_hotspots(
            current.sequence,
            hotspots, p.pk_free,
            p.pk_only, p.dangles,
            input_structure_given
        );

        output_results(
            current.sequence,
            results, p.fileO,
            p.suboptCount,
            current.name,
            inputs.size()
        );

        all_results.push_back(results);
    }
    return all_results;
}


const std::string traceback_cases_output =
    ">VP_case1\n"
    "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n"
    "(((((.......[[[[[[[[[[........)))))......]]]]]]]]]] (-20.16)\n"
    ">VP_case2\n"
    "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n"
    "[[[[[.......((((((((((........]]]]]......)))))))))) (-20.16)\n"
    ">VP_case3\n"
    "CCCCCCCCCCAAAAAAGGGGGAAAGGGGGGGGGGAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n"
    "((((((((((......[[[[[...))))))))))....((((((((((........]]]]]......)))))))))) (-34.62)\n"
    ">VP_case1_and_case2\n"
    "CCCCCCCCCCAAAAAAGGGGGAAAGGGGGGGGGGAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n"
    "[[[[[[[[[[......(((((...]]]]]]]]]]....[[[[[[[[[[........)))))......]]]]]]]]]] (-34.62)\n"
    ">Multiloop_with_pseudoknot\n"
    "AAAAAAGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCCUUUUUU\n"
    "(((((((((((.......[[[[[[[[[[........)))))......]]]]]]]]]])))))) (-8.71)\n"
    ">Multiloop_that_spans_a_band\n"
    "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAAAACCCCCC\n"
    "((((((.....[[[[[[....[[[[[[))))))......]]]]]]....((((((......))))))......]]]]]] (-25.71)\n"
    ">Multiloop_that_spans_a_band_with_extra\n"
    "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC\n"
    "((((((.....[[[[[[....[[[[[[))))))......]]]]]]....((((((......))))))...(((...)))]]]]]] (-26.14)\n"
    ">Multiloop_that_spans_a_band_with_extra2\n"
    "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC\n"
    "...........((((((....((((((............))))))....((((((......))))))...((....)).)))))) (-28.51)\n";


const std::string input_file_and_sequence_output =     
    ">Console Sequence\n"
    "AUGCUAGC\n"
    "........ (0)\n"
    ">VP_case1\n"
    "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n"
    "(((((.......[[[[[[[[[[........)))))......]]]]]]]]]] (-20.16)\n"
    ">VP_case2\n"
    "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n"
    "[[[[[.......((((((((((........]]]]]......)))))))))) (-20.16)\n"
    ">VP_case3\n"
    "CCCCCCCCCCAAAAAAGGGGGAAAGGGGGGGGGGAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n"
    "((((((((((......[[[[[...))))))))))....((((((((((........]]]]]......)))))))))) (-34.62)\n"
    ">VP_case1_and_case2\n"
    "CCCCCCCCCCAAAAAAGGGGGAAAGGGGGGGGGGAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n"
    "[[[[[[[[[[......(((((...]]]]]]]]]]....[[[[[[[[[[........)))))......]]]]]]]]]] (-34.62)\n"
    ">Multiloop_with_pseudoknot\n"
    "AAAAAAGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCCUUUUUU\n"
    "(((((((((((.......[[[[[[[[[[........)))))......]]]]]]]]]])))))) (-8.71)\n"
    ">Multiloop_that_spans_a_band\n"
    "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAAAACCCCCC\n"
    "((((((.....[[[[[[....[[[[[[))))))......]]]]]]....((((((......))))))......]]]]]] (-25.71)\n"
    ">Multiloop_that_spans_a_band_with_extra\n"
    "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC\n"
    "((((((.....[[[[[[....[[[[[[))))))......]]]]]]....((((((......))))))...(((...)))]]]]]] (-26.14)\n"
    ">Multiloop_that_spans_a_band_with_extra2\n"
    "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC\n"
    "...........((((((....((((((............))))))....((((((......))))))...((....)).)))))) (-28.51)\n";