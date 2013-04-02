#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <map>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <unistd.h>
#include "atom.h"

int generate_feff(std::vector< Atom > atoms, std::string feff, std::string ifeffit) {

	// Stored as index of metal
	std::vector< std::pair< int, bool > > metals;

	// Stored as atomic number
	std::map< std::string, std::pair< int, int > > unique_atoms;

	// Record all metals.
	for (int i = 0; i < (int)atoms.size(); ++i) {
		if (atoms[i].is_metal()) {
			metals.push_back(std::make_pair(i, false));
		}
	}

	// Find unique atoms
	// int counter = 1;
	for (int i = 0; i < (int)atoms.size(); ++i) {
		if (unique_atoms.find(atoms[i].name) == unique_atoms.end()) {
			unique_atoms[atoms[i].name] = std::make_pair(0,atoms[i].atomic_number);
		}
	}

	// Search metal list to see if any atoms are individuals.
	for (int i = 0; i < (int)metals.size(); ++i) {
		int atomic_number = atoms[metals[i].first].atomic_number;
		int counter = 0;
		for (int j = 0; j < (int)atoms.size(); ++j) {
			if (atoms[j].atomic_number == atomic_number) {
				counter++;
			}
		}
		if (counter == 1) {
			metals[i].second = true;
		}
	}

	// Loop foreach metal
	for (int i = 0; i < (int)metals.size(); ++i) {
		
		std::ostringstream oss;
		oss << i;

		std::string name = std::string("feff-") + oss.str() + std::string(".inp");
		std::ofstream ofile(name.c_str());

		if (ofile.is_open() && ofile.good()) {
			ofile << "TITLE My Molecule" << std::endl;
			ofile << std::endl;
			ofile << "POTENTIALS" << std::endl;
			ofile << "*	ipot	z	tag" << std::endl;

			// Write absorbing atom
			Atom absorbing_atom = atoms[metals[i].first];
			ofile << 0 << "\t" << absorbing_atom.atomic_number << "\t" << absorbing_atom.name << std::endl;

			// Write all unique atoms
			int counter = 1;
			for (std::map< std::string, std::pair< int, int > >::iterator iter = unique_atoms.begin(); iter != unique_atoms.end(); ++iter) {
				
				// Skip atom if this atom is a selected metal and unique to the molecule.
				if (metals[i].second && absorbing_atom.atomic_number == iter->second.second) {
					continue;
				}
				iter->second.first = counter++;
				ofile << iter->second.first << "\t" << iter->second.second << "\t" << iter->first << std::endl;
			}

			ofile << std::endl;
			ofile << "CONTROL 1 1 1 1 1 1" << std::endl;
			ofile << "NLEG 8" << std::endl;
			ofile << "PRINT 0 0 0 3" << std::endl;
			ofile << "EXCHANGE 0 0.0 0.0" << std::endl;
			ofile << "CRITERIA 4.0 2.5" << std::endl;
			ofile << std::endl;
			ofile << std::endl;
			ofile << "ATOMS" << std::endl;
			ofile << "*	x	y	z	ipot	tag" << std::endl;

			for (int j = 0; j < (int)atoms.size(); ++j) {

				// If current index equals the index of the selected metal.
				if (j == metals[i].first) {
					ofile << std::setprecision(15) << "\t" << atoms[j].x << "\t" << atoms[j].y << "\t" << atoms[j].z << "\t" << 0 << "\t" << atoms[j].name << std::endl;
				} else {
					ofile << std::setprecision(15) << "\t" << atoms[j].x << "\t" << atoms[j].y << "\t" << atoms[j].z << "\t" << unique_atoms[atoms[j].name].first << "\t" << atoms[j].name << std::endl;
				}
			}

			ofile << "END" << std::endl;
			ofile.close();
		} else {
			std::cout << "could not open" << std::endl;
		}
	}

	// Create script file to process feff files.
	std::ofstream fscript("fscript.sh");

	fscript << std::endl << "rm -rf ";
	for (int i = 0; i < (int)metals.size(); ++i) {
		fscript << i << " ";
	}
	fscript << std::endl;

	for (int i = 0; i < (int)metals.size(); ++i) {

		fscript << "mkdir " << i << std::endl;
		fscript << "mv " << std::string("feff-") << i << std::string(".inp") << " " << "./" << i << "/feff.inp" << std::endl;
		fscript << "cp chi.chi3 ./" << i << std::endl;
		fscript << "cd " << i << std::endl;
		fscript << feff << " > /dev/null" << std::endl;
		fscript << "bash ../process_ifeffit.sh" << std::endl;
		fscript << ifeffit << " -q process.iff" << std::endl;
		fscript << "cd .." << std::endl;
		fscript << std::endl;
	}

	std::string first_part = "";
	std::string second_part = "";
	for (int i = 0; i < (int)metals.size(); ++i) {

		fscript << "folder" << i << "=($(grep -v \"#\" ./" << i << "/my_chi.chi3))" << std::endl;

		std::ostringstream oss;
		oss << i;
		first_part += "${folder" + oss.str() + "[$i+1]}";

		std::ostringstream oss2;
		oss2 << (i+1);
		second_part += "$" + oss2.str();
		if (i < (int)metals.size()-1) {
			first_part += " + ";
			second_part += " + ";
		}
	}

	fscript << std::endl << "rm -rf dod.dat" << std::endl;

	fscript << std::endl << "tLen=${#folder0[@]}" << std::endl << std::endl;
	fscript << "for (( i=0; i<${tLen}; i=i+2 ));" << std::endl;
	fscript << "do" << std::endl;
	fscript << "\ttotal=$( echo \"" + first_part + "\" | awk -F \"+\" '{print " + second_part + " }')" << std::endl;
	fscript << "\tavg=$( echo $total / " << (int)metals.size() << " | awk -F \"/\" '{print $1 / $2}')" << std::endl;
	fscript << "\t$(echo \"${folder0[$i]} $avg\" >> dod.dat)" << std::endl;
	fscript << "done" << std::endl;

	fscript << std::endl << "rm -rf ";
	for (int i = 0; i < (int)metals.size(); ++i) {
		fscript << i << " ";
	}
	fscript << std::endl;

	fscript.close();

	return (int)metals.size();
}