#include "Result.hpp"
#include "cmdline.hpp"
#include "W_final.hpp"
#include "HFold.hpp"
#include "Hotspot.hpp"

#include <iostream>
#include <fstream>
#include <filesystem>
#include <algorithm>


int main(int argc, char* argv[]) {
    args_info a;
    if (cmdline_parser(argc, argv, &a) != 0) return 1;

	/* ────── Command‑line → variables ────── */
    std::string sequence   = getSequence(a);
    std::string restricted = a.input_structure_given ? a.input_structure_arg   : ""; // Initial structure
    std::string fileI      = a.input_file_given      ? a.input_file_arg        : "";
    std::string fileO      = a.output_file_given     ? a.output_file_arg       : "";
    int suboptCount        = a.opt_given             ? a.opt_arg               : 1;
    bool pk_free           = a.pk_free_given;
    bool pk_only           = a.pk_only_given;
    int dangles            = a.dangles_given         ? a.dangles_arg           : 1;
    std::string paramFile  = a.paramFile_given       ? a.paramFile_arg         : std::string(PARAMS_DIR) + "/rna_DirksPierce09.par";

    std::vector<RNAEntry> inputs = get_all_inputs(fileI, sequence, restricted);

    for (RNAEntry& current : inputs){
        preprocess_sequence(current.sequence, current.structure, a.noConv_given);
        load_energy_parameters(paramFile, current.sequence,a.paramFile_given);

        std::vector<Hotspot> hotspots = build_hotspots(
            current.sequence,
            current.structure,
            suboptCount
        );

        std::vector<Result>  results  = fold_hotspots(
            current.sequence,
            hotspots, pk_free,
            pk_only, dangles,
            a.input_structure_given
        );

        output_results(current.sequence, 
            results,
            fileO,
            suboptCount,
            current.name,
            inputs.size()
        );
    }

    cmdline_parser_free(&a);
    return 0;
}
