#include <iostream>
#include <fstream>
#include <map>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <utility>
#include <limits>
#include "main.h"
#include "ifeffit.h"

double calc_rmsd(Optim settings) {

	std::vector< std::pair<double, double> > calc_data;
	std::vector< std::pair<double, double> > exp_data;

	std::ifstream calc_file(settings.calc_data_file.c_str());
	std::ifstream exp_file(settings.exp_data_file.c_str());

	double rmsd = 0;
	if (calc_file.is_open() && exp_file.is_open()) {

		// Read in data from files.
		std::string x,y;
		while(calc_file.good()) {
			calc_file >> x >> y;
			calc_data.push_back(std::make_pair(atof(x.c_str()), atof(y.c_str())));
		}
		while(exp_file.good()) {
			exp_file >> x >> y;
			exp_data.push_back(std::make_pair(atof(x.c_str()), atof(y.c_str())));
		}
/*
		// Calculate steps in data.
		double diff_calc;
		double diff_exp;
		diff_calc = calc_data[1].first - calc_data[0].first;
		diff_exp = exp_data[1].first - exp_data[0].first;

		// Find initial start positions for each data.
*/
		int calc_index = 0;
		int exp_index = 0;
		while (calc_index < (int)calc_data.size() && exp_index < (int)exp_data.size()) {

			if (fabs(calc_data[calc_index].first - exp_data[exp_index].first) < std::numeric_limits<double>::epsilon()) {
				// std::cout << calc_data[calc_index].first << ", " << exp_data[exp_index].first << std::endl;

				if (calc_data[calc_index].first > settings.x_min && calc_data[calc_index].first < settings.x_max) {
					// std::cout << calc_data[calc_index].first << ", " << exp_data[exp_index].first << std::endl;
					double scale_calc = calc_data[calc_index].second * settings.calc_data_rescale + settings.calc_y_shift;
					double scale_exp = exp_data[exp_index].second * settings.exp_data_rescale + settings.exp_y_shift;
					rmsd += pow(scale_calc - scale_exp,2);
				}
				++calc_index;
				++exp_index;
			} else if (calc_data[calc_index].first < exp_data[exp_index].first) {
				++calc_index;
			} else {
				++exp_index;
			}
		}

	} else {
		rmsd = -1; // Error
	}

	calc_file.close();
	exp_file.close();
	return rmsd;
}

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

#include <sys/time.h>
#include <ctime>

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

int main() {

	// Molecule mole = Molecule::read_from_disk("coords.xyz");
	// int num = generate_feff(mole.atoms, "tester");


	Optim settings = readConfig("optim.inp");
	Molecule molecule = Molecule::read_from_disk("coords.xyz");

	std::string exafs_program = "bash " + settings.main_program;

	// std::ofstream results("results.dat");

	generate_feff(molecule.atoms,settings.feff,settings.ifeffit);
	system(exafs_program.c_str());

	timestamp_t t0 = get_timestamp();
    double eold = evaluate(settings);
    timestamp_t t1 = get_timestamp();
    double olds = (t1 - t0) / 1000000.0L;
	
	t0 = get_timestamp();
    double enew = calc_rmsd(settings);
    t1 = get_timestamp();
    double news = (t1 - t0) / 1000000.0L;

	std::cout << std::setprecision(15) << "original = " << olds << std::endl;
	std::cout << std::setprecision(15) << "new = " << news << std::endl;
	std::cout << std::setprecision(15) << "original = " << eold << std::endl;
	std::cout << std::setprecision(15) << "new = " << enew << std::endl;
	// std::cout << "Generation: 0" << std::endl;
	
	/*
	molecule.cost = evaluate(settings);

	// Generations
	for (int i = 0; i < 10; ++i) {
		
		std::cout << "Generation: " << (i+1) << std::endl;

		// Manipulate the molecule to create a new one.
		Molecule moved_molecule = manipulate(settings,molecule);

		// Write new molecule to disk for exafs program.
		// moved_molecule.write_to_disk("coords.xyz");

		generate_feff(moved_molecule.atoms,settings.feff,settings.ifeffit);
		system(exafs_program.c_str());

		// Evaluate the cost of this new molecule.
		moved_molecule.cost = evaluate(settings);

		// If the new molecule is better than the old we will save that.
		molecule = compare(molecule,moved_molecule);

		results << std::setprecision(15) <<	 molecule.cost << std::endl;
	}

	// Write out final positions.
	molecule.write_to_disk("coords.xyz");

	results.close();
	*/
	return 0;
}
