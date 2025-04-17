#include "Result.hpp"

//constructor
Result::Result(std::string sequence,std::string restricted,double restricted_energy, std::string final_structure, double final_energy){
    this->sequence = sequence;
    this->restricted = restricted;
    this->restricted_energy = restricted_energy;
    this->final_structure = final_structure;
    this->final_energy = final_energy;

}

//destructor
Result::~Result(){
   
}

std::string Result::get_sequence() const {
    return this->sequence;
}
std::string Result::get_restricted() const {
    return this->restricted;
}
std::string Result::get_final_structure() const {
    return this->final_structure;
}
double Result::get_final_energy() const {
    return this->final_energy;
}

double Result::get_restricted_energy() const {
    return this->restricted_energy;
}
