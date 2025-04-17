#include <gtest/gtest.h>
#include <string.h>
#include <tuple>
#include <HFold.hpp>
#include <sparse_tree.hpp>
#include <W_final.hpp>

void load_parameters(std::string param_file){
    if (exists(param_file)) {
        vrna_params_load(param_file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
    } else {
        FAIL() << "Parameter file not found: " << param_file;
    }

    // Scale parameters (matches main behavior)
    vrna_param_s* params = scale_parameters();
}

struct FoldResult {
    std::string final_structure;
    double energy;
};

FoldResult simple_hfold(std::string seq, std::string structure, double energy = 0.0, sparse_tree* tree = nullptr, bool pk_free = false, bool pk_only = false, int dangles = 2){
    load_parameters("../../params/rna_DirksPierce09.par");
    sparse_tree local_tree(structure, seq.length());

    // if the tree is not provided, create a new one
    if (tree == nullptr) {
        tree = &local_tree;
    }

    std::string final_structure = hfold(seq, structure, energy, *tree, pk_free, pk_only, dangles);
    return FoldResult{final_structure, energy};
}

TEST(SimpleHFold, ZeroSequence){
    std::string seq = "AUGCUAGC";
    std::string structure = "..........";
    FoldResult result = simple_hfold(seq, structure);
    EXPECT_EQ(result.final_structure, "........");
    EXPECT_EQ(result.energy, 0);
}

TEST(SimpleHFold, Positive){
    std::string seq = "AAAAGGAAUUUC";
    std::string structure = ".(((....))).";
    FoldResult result = simple_hfold(seq, structure);
    EXPECT_EQ(result.final_structure, ".(((....))).");
    EXPECT_EQ(result.energy, 0.87);
}

TEST(SimpleHFold, Knotted){
    std::string seq = "AGGGCUAUCCUU";
    std::string structure = "((((...)))).";
    FoldResult result = simple_hfold(seq, structure);
    EXPECT_EQ(result.final_structure, "((((..[))))]");
    EXPECT_EQ(result.energy, -0.49);
}

TEST(SimpleHFold, KnottedNotDensity2){
    std::string seq = "CCCCCCCCCCAAAAAAGGGGGAAAGGGGGGGGGGAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "........................((((((((((.................................))))))))))";
    FoldResult result = simple_hfold(seq, structure);
    EXPECT_EQ(result.final_structure, "....(((((.......)))))...((((((((((......(((((...........)))))......))))))))))");
    EXPECT_EQ(result.energy, -32.4);
}

TEST(SimpleHFold, Unknotted){
    std::string seq = "CCAAGGAAGGCCAA";
    std::string structure = "................";
    FoldResult result = simple_hfold(seq, structure);
    EXPECT_EQ(result.final_structure, "....((....))..");
    EXPECT_EQ(result.energy, -1.48);
}



