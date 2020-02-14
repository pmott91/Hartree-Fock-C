#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "HF.hpp"


using std::fixed;
using std::setprecision;

//Define main program
int main(int argc, char *argv[]) {


//Determine the source of the input data.
//The data file can be given via the command line.
//By default, the program reads in h2.xyz.
const auto xyz_file = (argc > 1) ? argv[1] : "h2.xyz";
int atomn=0;

//Define the variable 'atoms' by calling read_xyz.
std::vector<Atom> atoms = read_xyz(xyz_file,atomn);

//Define floating point format for output.
std::cout << std::fixed;

//Print the input data from the XYZ file.
int i;
for (i = 0; i < atomn; i++){
std::cout << atoms[i].atomic_number << "  ";
std::cout << setprecision(9) << atoms[i].x << "  ";
std::cout << atoms[i].y << "  ";
std::cout << atoms[i].z << "\n";
}

//Basis Set Setup

int STOnG = 3;

// Max Quantum Number

int H =1;

//Slater Zeta value

double Zeta = 1.24;

//Gaussian contraction coefficients

double D[2][3];

D[0][0]=0.444635;
D[0][1]=0.535328;
D[0][2]=0.154329;
D[1][0]=0.700115;
D[1][1]=0.399513;
D[1][2]=-0.0999672;

//Gaussian orbital exponents

double alpha[2][3];

alpha[0][0]=0.109818;
alpha[0][1]=0.405771;
alpha[0][2]=2.22766;
alpha[1][0]=0.0751386;
alpha[1][1]=0.231031;
alpha[1][2]=0.994203;

//Basis set size

int B=2;

//Number of electrons

int N =2;

//Compute the nuclear-repulsion energy.

double enuc = 0.0;

for (int k = 0; k < atomn; k++) {

       for (int l = k +1; l < atomn; l++) {
       double xkl = atoms[k].x - atoms[l].x;
       double ykl = atoms[k].y - atoms[l].y;
       double zkl = atoms[k].z - atoms[l].z;
       double r2 = xkl*xkl + ykl*ykl + zkl*zkl;
       double r = sqrt(r2);       
       
       enuc += atoms[k].atomic_number * atoms[l].atomic_number/r;

       }
}

std::cout << "The nuclear repulsion energy is: " << enuc << std::endl;

//Compute the overlap integrals




//Initialize matrices
/*
double S[B][B];
double T[B][B];
double V[B][B];
double multi_electron_tensor[B][B][B][B];

int j;
double Za[atomn];
std::vector<double> Ra;
std::vector<double> d_vec;
for (j=0; j < atomn; j++){
Za[j] = atoms[j].atomic_number;
Ra = {atoms[j].x,atoms[j].y,atoms[j].z};

d_vec = {D[0][0],D[0][1],D[0][2]};


*/



}


