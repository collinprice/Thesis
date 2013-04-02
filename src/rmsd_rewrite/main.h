#include "molecule.h"

// Structs
struct Optim {
	std::string main_program;
	std::string calc_data_file;
	std::string exp_data_file;
	int num_of_calc_data;
	int num_of_exp_data;
	double calc_data_rescale;
	double exp_data_rescale;
	double calc_y_shift;
	double exp_y_shift;
	double x_min;
	double x_max;
	int cycle_number;
	double dr_max; //Drmax0
	double tolerance;
	double temperature;
	bool annealing;
	std::string temperature_program; // expo / lin
	double k_temp;
	int total_annealings;
	int number_of_t_points;
	bool give_atom_list;
	bool give_bond_list;
	bool use_harmonic_restriction;
	double universal_harmonic_constant;
	int number_of_atoms;
	std::string feff;
	std::string ifeffit;
};

/*
	Evaluate the total energy of a molecule.
*/
double evaluate(Optim settings);

double calc_rmsd(Optim settings);

/*
	Move the atoms within a molecule.
*/
Molecule manipulate(Optim settings, Molecule molecule);

/*
	Read in configuration file.
*/
Optim readConfig(std::string filename);

Molecule compare(Molecule m1, Molecule m2);
double unifRand();

void createFEFF(std::vector<Atom> atoms, std::string filename);
bool contains_atom(std::vector<Atom> &atoms, Atom &atom);
