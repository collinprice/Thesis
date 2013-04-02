#include <string>

class Atom {

	public:
		std::string name;
		double x,y,z;
		int atomic_number;

		Atom(std::string name, double x, double y, double z);
		std::string to_string();
		std::string file_string();
		bool is_metal();
		static std::string get_name(int atomic_number);

	private:
		const static int periodic_table_size = 103;
		const static char* periodic_table[];

		const static int periodic_table_metals[];
};
