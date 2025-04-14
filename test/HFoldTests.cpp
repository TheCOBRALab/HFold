#include <gtest/gtest.h>
#include <string.h>
#include "../src/HFold.hh"
#include "../src/sparse_tree.hh"
#include "../src/W_final.hh"


void load_parameters(std::string param_file){
    if (exists(param_file)) {
        vrna_params_load(param_file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
    } else {
        FAIL() << "Parameter file not found: " << param_file;
    }

    // Scale parameters (matches main behavior)
    vrna_param_s* params = scale_parameters();
}


TEST(SimpleHFold, ZeroSequence){
    load_parameters("../../params/rna_DirksPierce09.par");
    std::string seq = "AUGCUAGC";
    std::string structure = "..........";
    double energy;
    sparse_tree tree(structure, seq.length());
    std::string final_structure = hfold(seq, structure, energy, tree, false, false, 2);
    
    EXPECT_EQ(final_structure, "........");
    EXPECT_EQ(energy, 0);
}

TEST(SimpleHFold, Knotted){
    load_parameters("../../params/rna_DirksPierce09.par");
    std::string seq = "AGGGCUAUCCUU";
    std::string structure = "((((...)))).";
    double energy;
    sparse_tree tree(structure, seq.length());
    std::string final_structure = hfold(seq, structure, energy, tree, false, false, 2);
    
    EXPECT_EQ(final_structure, "((((..[))))]");
    EXPECT_EQ(energy, -0.49);
}

