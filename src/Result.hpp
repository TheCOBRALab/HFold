#ifndef RESULT_HEADER
#define RESULT_HEADER

#include <string>

class Result{
    public:
        //constructor
        Result(std::string sequence, std::string restricted, double restricted_energy, std::string final_structure, double final_energy);
        //destructor
        ~Result();

        //getter
        std::string get_sequence() const;
        std::string get_restricted() const;
        std::string get_final_structure() const;
        double get_restricted_energy() const;
        double get_final_energy() const;

        struct Result_comp{
		bool operator ()(Result &x, Result &y) const {
			if (x.get_final_energy() != y.get_final_energy())
                return x.get_final_energy() < y.get_final_energy();
            return x.get_restricted_energy() < y.get_restricted_energy();
		}
		} result_comp;
        
    private:
        std::string sequence;
        std::string restricted;
        double restricted_energy;
        std::string final_structure;
        double final_energy;
};



#endif
