#include "molecule.h"
#include <fstream>
#include <stdlib.h>
#include <iostream>

Molecule::Molecule(std::vector<Atom> atoms) {
	this->atoms = atoms;
	this->isEvaluated = false;
	this->cost = 0;
}

int Molecule::size() {
	return (int)this->atoms.size();
}

// Utility
bool Molecule::write_to_disk(std::string filename) {
	std::ofstream output(filename.c_str());
	if (output.is_open()) {
		output << this->size() << std::endl;
		output << std::endl;

		for (std::vector<Atom>::iterator i = this->atoms.begin(); i != this->atoms.end(); ++i) {
			output << (*i).file_string() << std::endl;
		}

		output.close();
	}

	return false;
}

Molecule Molecule::read_from_disk(std::string filename) {

	std::ifstream input(filename.c_str());

	if (!input.is_open()) {
		std::cout << "Molecule file not found." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Read number of atoms.
	std::string temp;
	int count;
	input >> temp;
	count = atoi(temp.c_str());

	// Read individual atoms.
	std::vector<Atom> atoms;
	for (int i = 0; i < count; ++i) {
		std::string name;
		double x,y,z;

		input >> name;

		input >> temp;
		x = atof(temp.c_str());
		input >> temp;
		y = atof(temp.c_str());
		input >> temp;
		z = atof(temp.c_str());

		atoms.push_back(Atom(name,x,y,z));
	}

	input.close();

	return Molecule(atoms);
}