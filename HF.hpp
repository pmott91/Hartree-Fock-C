// Build a struct to store atomic number and coordinates.
struct Atom {
       int  atomic_number;
       double x, y, z;
       };
//Declare the function read_xyz which will collect coordinates 
//and number of atoms from an XYZ file.
std::vector<Atom> read_xyz(const std::string& xyz_file,int& atomn);

