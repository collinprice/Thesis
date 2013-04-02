#include "atom.h"
#include <sstream>
#include <iomanip>
#include <string.h>
#include <algorithm>

Atom::Atom(std::string name, double x, double y, double z) {

	this->name = name;
	this->x = x;
	this->y = y;
	this->z = z;

	for (int i = 0; i < this->periodic_table_size; ++i) {
		if (strcmp(this->periodic_table[i], this->name.c_str()) == 0) {
			this->atomic_number = (i+1);
		}
	}
}

const char* Atom::periodic_table[] = {
	"H", "He","Li","Be","B", "C", "N", "O", "F", "Ne",
	"Na","Mg","Al","Si","P", "S", "Cl","Ar","K", "Ca",
	"Sc","Ti","V", "Cr","Mn","Fe","Co","Ni","Cu","Zn",
	"Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y", "Zr",
	"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
	"Sb","Te","I", "Xe","Cs","Ba","La","Ce","Pr","Nd",
	"Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
	"Lu","Hf","Ta","W", "Re","Os","Ir","Pt","Au","Hg",
	"Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
	"Pa","U", "Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
	"Md","No","Lr"
};

const int Atom::periodic_table_metals[] = {
	3, 4,
	11, 12, 13,
	19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
	37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
	55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
	87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103
};

std::string Atom::to_string() {

	std::ostringstream strs;
	strs << this->name << " (" << this->x << ", " << this->y << ", " << this->z << ")";

	return strs.str();
}

std::string Atom::file_string() {

	std::ostringstream strs;
	strs << std::setprecision(14) << this->name << "\t" << this->x << "\t" << this->y << "\t" << this->z;

	return strs.str();
}

bool Atom::is_metal() {
	int size = sizeof(this->periodic_table_metals)/sizeof(this->periodic_table_metals[0]);
	for (int i = 0; i < size; ++i) {
		if (this->atomic_number == this->periodic_table_metals[i]) return true;
	}

	return false;
}

std::string Atom::get_name(int atomic_number) {
	if (atomic_number <= Atom::periodic_table_size && atomic_number >= 0) {
		return Atom::periodic_table[atomic_number];
	}

	return "";
}