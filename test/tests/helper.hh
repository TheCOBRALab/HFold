#ifndef HELPER_H
#define HELPER_H
#include <string>
#include <SHAPE.hh>

struct HFoldParams {
    std::string sequence = "";
    std::string restricted = "";
    double energy = 0.0;
    bool pk_free = false;
    bool pk_only = false;
    int dangles = 2;
    int suboptCount = 1;
    SHAPEData shape_data{"", (cand_pos_t) sequence.size()};
    std::string fileI = "";
    std::string fileO = "";
    std::string paramFile = "rna_DirksPierce09.par";
    bool noConv_given = false;
};

std::string capture_console_output(const std::function<void()>& func);
std::string get_output(const std::string& file_path, bool delete_after_read = false);
std::vector<std::vector<Result>> hfold_test(HFoldParams& p);

// terminal output strings
extern const std::string traceback_cases_output;
extern const std::string input_file_and_sequence_output;
extern const std::string input_and_output_file_output;

#endif