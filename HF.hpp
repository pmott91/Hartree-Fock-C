// Build a struct to store atomic number and coordinates.
struct Atom {
       int  atomic_number;
       double x, y, z;
       };
//Declare the function read_xyz which will collect coordinates 
//and number of atoms from an XYZ file.
std::vector<Atom> read_xyz(const std::string& xyz_file,int& atomn);

//Build a struct to store variables associated with a Gaussian.

struct gauss{
       double alph;
       double x, y, z;
       };

//Declare the function gauss_product which will calculate
//the product of two Gaussians.

std::vector<double> gauss_product(gauss gauss_a,gauss gauss_b);

double overlap(gauss GA, gauss GB);

double kinetic(gauss KA, gauss KB);

double fo(double t);

double potential(gauss NA, gauss NB, std::vector<double> CAtom);

double multi(gauss MA, gauss MB, gauss MC, gauss MD);
