#include <iostream>
#include <fstream>
#include <map>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <utility>
#include "main.h"
#include "ifeffit.h"

double evaluate(Optim settings) {

	int n_new;
	double cost;

	std::vector<double> xd;
	std::vector<double> yd;
	std::vector<double> xe;
	std::vector<double> ye;

	std::vector<double> yd_new;
	std::vector<double> ye_new;
	std::vector<double> x_new;

	int n_data = settings.num_of_calc_data;
	int n_exp = settings.num_of_exp_data;

	// std::cout << "Reading data files: " << settings.calc_data_file << ", " << settings.exp_data_file << std::endl;

	std::ifstream data_file, exp_file;
	data_file.open(settings.calc_data_file.c_str());
	exp_file.open(settings.exp_data_file.c_str());

	// Read in calculated data
	for (int i = 0; i < n_data; ++i) {
		std::string temp_x, temp_y;

		data_file >> temp_x >> temp_y;
		xd.push_back(atof(temp_x.c_str()));
		yd.push_back(atof(temp_y.c_str()));
	}

	// Read in experimental data
	for (int i = 0; i < n_exp; ++i) {
		std::string temp_x, temp_y;

		exp_file >> temp_x >> temp_y;
		xe.push_back(atof(temp_x.c_str()));
		ye.push_back(atof(temp_y.c_str()));
	}

	// Clean up files.
	data_file.close();
	exp_file.close();

	/*
	 If the highest x-value in the calculated data is larger then 
	 the highest x-value in the experimental data.
	 */
	if ( std::max( xe[n_exp-1], xd[n_data-1] ) == xe[n_exp-1] ) {
		
		// Init ye_new
		for (int i = 0; i < n_data; ++i) {
			ye_new.push_back(0);
		}

		n_new = n_data;
		x_new = xd;
		yd_new = yd;

		for (int i = 0; i < n_data; ++i) {
			for (int j = 0; j < n_exp-1; ++j) {
				if (x_new[i] <= xe[j+1] && x_new[i] >= xe[j]) {
					ye_new[i] = ye[j]+(x_new[i]-xe[j])*(ye[j+1]-ye[j])/(xe[j+1]-xe[j]);
				}
			}
		}
	} else {
		std::cout << "top bottom" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Scale new values
	for (int i = 0; i < (int)ye_new.size(); ++i) {
		ye_new[i] = ye_new[i] * settings.exp_data_rescale + settings.exp_y_shift;
	}

	for (int i = 0; i < (int)yd_new.size(); ++i) {
		yd_new[i] = yd_new[i] * settings.calc_data_rescale + settings.calc_y_shift;
	}

	// Calculate the cost function
	cost = 0;
	for (int i = 0; i < n_new; ++i) {
		if(x_new[i] > settings.x_min && x_new[i] < settings.x_max) {
			cost = cost + pow(ye_new[i]-yd_new[i],2);
		}
	}

	// std::cout << std::setprecision(15) << "cost " << cost << std::endl;

	return cost;
}

Molecule manipulate(Optim settings, Molecule molecule) {

	const int random_index = unifRand() * molecule.size();
	
	Atom atom = molecule.atoms[random_index];

	atom.x += (2 * unifRand() - 1) * settings.dr_max;
	atom.y += (2 * unifRand() - 1) * settings.dr_max;
	atom.z += (2 * unifRand() - 1) * settings.dr_max;

	molecule.atoms[random_index] = atom;

	return molecule;
}

Molecule compare(Molecule m1, Molecule m2) {
	// std::cout << m1.cost << " - " << m2.cost << std::endl;
	return m1.cost > m2.cost ? m2 : m1;
}

Optim readConfig(std::string filename) {

	Optim settings;
	std::string temp;

	std::ifstream input(filename.c_str());

	if (!input.is_open()) {
		std::cout << "Config file not found." << std::endl;
		exit(EXIT_FAILURE);
	}

	input >> settings.main_program;
	input.ignore(256, '\n');
	input >> settings.calc_data_file;
	input.ignore(256, '\n');
	input >> settings.exp_data_file;
	input.ignore(256, '\n');

	input >> temp;
	settings.num_of_calc_data = atoi(temp.c_str());
	input.ignore(256, '\n');
	input >> temp;
	settings.num_of_exp_data = atoi(temp.c_str());
	input.ignore(256, '\n');

	input >> temp;
	settings.calc_data_rescale = atof(temp.c_str());
	input.ignore(256, '\n');
	input >> temp;
	settings.exp_data_rescale = atof(temp.c_str());
	input.ignore(256, '\n');

	input >> temp;
	settings.calc_y_shift = atof(temp.c_str());
	input.ignore(256, '\n');
	input >> temp;
	settings.exp_y_shift = atof(temp.c_str());
	input.ignore(256, '\n');

	input >> temp;
	settings.x_min = atof(temp.c_str());
	input.ignore(256, '\n');
	input >> temp;
	settings.x_max = atof(temp.c_str());
	input.ignore(256, '\n');

	input >> temp;
	settings.cycle_number = atoi(temp.c_str());
	input.ignore(256, '\n');

	input >> temp;
	settings.dr_max = atof(temp.c_str());
	input.ignore(256, '\n');
	input >> temp;
	settings.tolerance = atof(temp.c_str());
	input.ignore(256, '\n');
	input >> temp;
	settings.temperature = atof(temp.c_str());
	input.ignore(256, '\n');

	input >> temp;
	settings.annealing = temp.compare("Y") == 0 ? true : false;
	input.ignore(256, '\n');

	input >> settings.temperature_program;
	input.ignore(256, '\n');

	input >> temp;
	settings.k_temp = atof(temp.c_str());
	input.ignore(256, '\n');

	input >> temp;
	settings.total_annealings = atoi(temp.c_str());
	input.ignore(256, '\n');
	input >> temp;
	settings.number_of_t_points = atoi(temp.c_str());
	input.ignore(256, '\n');

	input >> temp;
	settings.give_atom_list = temp.compare("Y") == 0 ? true : false;
	input.ignore(256, '\n');
	input >> temp;
	settings.give_bond_list = temp.compare("Y") == 0 ? true : false;
	input.ignore(256, '\n');
	input >> temp;
	settings.use_harmonic_restriction = temp.compare("Y") == 0 ? true : false;
	input.ignore(256, '\n');

	input >> temp;
	settings.universal_harmonic_constant = atof(temp.c_str());
	input.ignore(256, '\n');

	input >> settings.feff;
	input.ignore(256, '\n');

	input >> settings.ifeffit;
	input.ignore(256, '\n');

	input.close();

	return settings;
}

double unifRand() {
	return rand() / double(RAND_MAX);
}


int main() {

	Optim settings = readConfig("optim.inp");
	Molecule molecule = Molecule::read_from_disk("coords.xyz");

	std::string exafs_program = "bash " + settings.main_program;

	std::ofstream results("results.dat");
	results << "Current, Best" << std::endl;

	generate_feff(molecule.atoms,settings.feff,settings.ifeffit);
	system(exafs_program.c_str());
	molecule.cost = evaluate(settings);


	// Simulated Annealing - From Sheridan

	int generation = 0;
	int max_generations = 1000;
	double initial_temperature = 1000;
	double current_temperature = initial_temperature;
	Molecule best_molecule = molecule;
	Molecule current_molecule = molecule;

	while (generation < max_generations) {

		std::cout << "Generation: " << generation << std::endl;

		// Generate and evaluate new molecule
		Molecule new_molecule = manipulate(settings,current_molecule);
		generate_feff(new_molecule.atoms,settings.feff,settings.ifeffit);
		system(exafs_program.c_str());
		new_molecule.cost = evaluate(settings);

		if (new_molecule.cost < current_molecule.cost) {
			current_molecule = new_molecule;
			if (current_molecule.cost < best_molecule.cost) {
				best_molecule = current_molecule;
			}
		} else {
			if (unifRand() < exp((new_molecule.cost-current_molecule.cost)/current_temperature)) {
				current_molecule = new_molecule;
			}
		}

		results << std::setprecision(15) << best_molecule.cost << "," << current_molecule.cost << std::endl;
		std::cout << std::setprecision(15) << "Best:" << best_molecule.cost << ", Current: " << current_molecule.cost << std::endl;

		generation++;
		current_temperature = initial_temperature*pow(1/initial_temperature, generation/max_generations);
	}

	results.close();
	return 0;
}
