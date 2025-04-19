#include <gtest/gtest.h>
#include <string.h>
#include <tuple>
#include <HFold.hpp>
#include <sparse_tree.hpp>
#include <W_final.hpp>
#include <filesystem>

struct FoldResult {
    std::string name;
    std::string sequence;
    std::string structure;
    double energy;
    
};

std::vector<std::vector<Result>> hfold_test(
    std::string sequence = "",
    std::string restricted = "",
    double energy = 0.0,
    sparse_tree* tree = nullptr,
    bool pk_free = false,
    bool pk_only = false,
    int dangles = 2,
    int suboptCount = 1,
    std::string fileI = "",
    std::string fileO= "",
    std::string paramFile = "../../params/rna_DirksPierce09.par",
    bool noConv_given = false
){  
    if (!fileI.empty()) fileI = "./input/" + fileI;
    if (!fileO.empty()) fileO = "./output/" + fileI;
    bool input_structure_given = !restricted.empty();


    std::vector<RNAEntry> inputs = get_all_inputs(fileI, sequence, restricted);
    std::vector<std::vector<Result>> all_results;
    for (RNAEntry current : inputs){
        preprocess_sequence(current.sequence, current.structure, noConv_given);
        load_energy_parameters(paramFile, current.sequence);
        std::vector<Hotspot> hotspots = build_hotspots(current.sequence, current.structure, suboptCount);
        std::vector<Result>  results  = fold_hotspots(current.sequence, hotspots, pk_free, pk_only, dangles, input_structure_given);
        output_results(current.sequence, results, fileO, suboptCount, current.name, inputs.size());
        all_results.push_back(results);
    }
    return all_results;
}

TEST(SimpleHFold, ZeroSequence){
    std::string sequence = "AUGCUAGC";
    Result result = hfold_test(sequence)[0][0];
    EXPECT_EQ(result.get_final_structure(), "........");
    EXPECT_EQ(result.get_final_energy(), 0);
}

TEST(SimpleHFold, Positive){ // Normally positive, but is converted to 0
    std::string sequence = "AAAAGGAAUUUC";
    Result result = hfold_test(sequence)[0][0];
    EXPECT_EQ(result.get_final_structure(), "............");
    EXPECT_EQ(result.get_final_energy(), 0);
}

TEST(SimpleHFold, Knotted){
    std::string seq = "AGGGCUAUCCUU";
    Result result = hfold_test(seq)[0][0];
    EXPECT_EQ(result.get_final_structure(), "((((..[))))]");
    EXPECT_EQ(result.get_final_energy(), -0.49);
}

TEST(SimpleHFold, KnottedNotDensity2){
    std::string sequence = "CCCCCCCCCCAAAAAAGGGGGAAAGGGGGGGGGGAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    Result result = hfold_test(sequence)[0][0];
    EXPECT_EQ(result.get_final_structure(), "....(((((.......)))))...((((((((((......(((((...........)))))......))))))))))");
    EXPECT_EQ(result.get_final_energy(), -32.4);
}

TEST(SimpleHFold, Unknotted){
    std::string sequence = "CCAAGGAAGGCCAA";
    Result result = hfold_test(sequence)[0][0];
    EXPECT_EQ(result.get_final_structure(), "....((....))..");
    EXPECT_EQ(result.get_final_energy(), -1.48);
}



