#include "Eigen/Eigenvalues"
#include "Eigen/Dense"

// Build a struct to store atomic number and coordinates.
struct Atom {
       int  atomic_number;
       double x, y, z;
       };
//Declare the function read_xyz which will collect coordinates 
//and number of atoms from an XYZ file.
std::vector<Atom> read_xyz(const std::string& xyz_file,int& atomn);

//Declare the function max_quantum_number which will determine
//the maximum quantum number for each atom.

std::vector<int> max_quantum_number(std::vector<int>& atomz,int& atomn);

//Build a struct to store variables associated with a Gaussian.

struct gauss{
       double alph;
       double x, y, z;
       };

//Declare the function gauss_product which will calculate
//the product of two Gaussians.

std::vector<double> gauss_product(gauss gauss_a,gauss gauss_b);

//Declare the function overlap which will
//calculate the overlap integral.
double overlap(gauss GA, gauss GB);

//Declare the function kinetic which will
//calculate the kinetic integral.
double kinetic(gauss KA, gauss KB);

//Declare the Fo function.
double fo(double t);

//Declare the potential function.
double potential(gauss NA, gauss NB, Atom Catom);

//Declare the multi-electron tensor
double multi(gauss MA, gauss MB, gauss MC, gauss MD);

//Declare the succ_dens iterator.
double succ_dens(Eigen::MatrixXd ptilde, Eigen::MatrixXd p, int B);
