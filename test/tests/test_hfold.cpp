#include <gtest/gtest.h>
#include <string>
#include <tuple>
#include <filesystem>
#include <sstream>
#include <iostream>
#include <functional>

#include <HFold.hpp>
#include <W_final.hpp>
#include "helper.hpp"

TEST(SimpleHFold, ZeroSequence){
    HFoldParams params;
    params.sequence = "AUGCUAGC";
    Result result = hfold_test(params)[0][0];
    EXPECT_EQ(result.get_final_structure(), "........");
    EXPECT_EQ(result.get_final_energy(), 0);
}

TEST(SimpleHFold, Positive){ // Normally positive, but is converted to 0
    HFoldParams params;
    params.sequence = "AAAAGGAAUUUC";
    Result result = hfold_test(params)[0][0];
    EXPECT_EQ(result.get_final_structure(), "............");
    EXPECT_EQ(result.get_final_energy(), 0);
}

TEST(SimpleHFold, Knotted){
    HFoldParams params;
    params.sequence = "AGGGCUAUCCUU";
    Result result = hfold_test(params)[0][0];
    EXPECT_EQ(result.get_final_structure(), "((((..[))))]");
    EXPECT_EQ(result.get_final_energy(), -0.49);
}

TEST(SimpleHFold, KnottedNotDensity2){ // This is knotted in knotty but not in HFold
    HFoldParams params;
    params.sequence = "CCCCCCCCCCAAAAAAGGGGGAAAGGGGGGGGGGAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    Result result = hfold_test(params)[0][0];
    EXPECT_EQ(result.get_final_structure(), "....(((((.......)))))...((((((((((......(((((...........)))))......))))))))))");
    EXPECT_EQ(result.get_final_energy(), -32.4);
}

TEST(SimpleHFold, Unknotted){
    HFoldParams params;
    params.sequence = "CCAAGGAAGGCCAA";
    Result result = hfold_test(params)[0][0];
    EXPECT_EQ(result.get_final_structure(), "....((....))..");
    EXPECT_EQ(result.get_final_energy(), -1.48);
}

TEST(SimpleHFold, Suboptimal){
    HFoldParams params;
    params.sequence = "AGGGCUAUCCUU";
    params.suboptCount = 2;
    std::vector<Result> result = hfold_test(params)[0];

    EXPECT_EQ(result[0].get_final_structure(), "((((..[))))]");
    EXPECT_EQ(result[0].get_final_energy(), -0.49);
    EXPECT_EQ(result[1].get_final_structure(), "((((....))))");
    EXPECT_EQ(result[1].get_final_energy(), -0.33);
}

TEST(SimpleHFold, DanglesAndRestriction) {
    HFoldParams params;
    params.sequence = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    params.restricted = "(((((.........................)))))................";
    params.dangles = 1;
    Result result = hfold_test(params)[0][0];
    EXPECT_EQ(result.get_final_structure(), "(((((.......[[[[[[[[[[........)))))......]]]]]]]]]]");
    EXPECT_EQ(result.get_final_energy(), -20.16);
}

TEST(FileHFold, InputFileOnly){
    HFoldParams params;

    params.fileI = "input/traceback_cases.txt";
    std::string output = capture_console_output([&](){
        hfold_test(params);
    });

    // This has weird formating to test the parser
    params.fileI = "input/traceback_cases_bad_formating.txt"; 
    std::string output2 = capture_console_output([&](){
        hfold_test(params);
    });

    EXPECT_EQ(output , traceback_cases_output);
    EXPECT_EQ(output2, traceback_cases_output);
}

TEST(FileHFold, InputFileAndSequence) {
    HFoldParams params;
    params.sequence = "AUGCUAGC";
    params.fileI = "input/traceback_cases.txt";
    std::string output = capture_console_output([&](){
        hfold_test(params);
    });

    EXPECT_EQ(output , input_file_and_sequence_output);
}

TEST(FileHFold, InputAndOutputFile) {
    HFoldParams params;
    params.fileI = "input/traceback_cases.txt";
    params.fileO = "traceback_cases_out.txt";
    hfold_test(params);
    std::string output   = get_output(params.fileO);
    std::string expected = get_output("expected_output/traceback_cases_expected.txt" );
    EXPECT_EQ(output , expected);
}


/* FAILS ON MAC. WORKING ON FIX */
// TEST(FileHFold, InputAndOutputFileWithSuboptimal) {
//     HFoldParams params;
//     params.fileI = "../../test/tests/input/traceback_cases.txt";
//     params.fileO = "traceback_cases_subopt_out.txt";
//     params.suboptCount = 3;
//     hfold_test(params);
//     std::string output   = get_output(params.fileO);
//     std::string expected = get_output("../../test/tests/expected_output/traceback_cases_subopt_expected.txt" );
//     EXPECT_EQ(output , expected);
// }