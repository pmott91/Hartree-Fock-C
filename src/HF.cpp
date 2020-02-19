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

// For each atom determine the max quantum number.
int i2;
std::vector<int> atomz;
for (i2 = 0; i2 < atomn; i2++){
atomz.push_back(atoms[i2].atomic_number);
}

std::vector<int> nlist = max_quantum_number(atomz,atomn);

int i3;
for (i3 = 0; i3 < atomn; i3++){
std::cout << "Max quantum number for atom " << i3+1 << " is " << nlist[i3] << std::endl;
}

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

int B=0;
for (int i4=0; i4 < atomn; i4++) {
B+=atomz[i4];
}

std::cout << "Basis set size is " << B << std::endl;

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

//Initialize matrices

double S[B][B];
double T[B][B];
double V[B][B];
double multi_electron_tensor[B][B][B][B];

//Iterate through atoms.

int j;
double Za[atomn];
double Ra[atomn][3];
std::vector<double> d_vec_m;
std::vector<double> alpha_vec_m;
double Zb[atomn];
double Rb[atomn][3];
std::vector<double> d_vec_n;
std::vector<double> alpha_vec_n;
int a;
int b;
for (j=0; j < atomn; j++){
       
       //For each atom, get the charge and center.
       Za[j] = atoms[j].atomic_number;
       Ra[j][0]=atoms[j].x;
       Ra[j][1]=atoms[j].y;
       Ra[j][2]=atoms[j].z;

       //For each quantum number, get the contraction
       //coefficients, then get zeta,
       //then scale the exponents accordingly. 

       for(int j2=0; j2 < nlist[j]; j2++){

              d_vec_m = {D[j2][0],D[j2][1],D[j2][2]};
              alpha_vec_m = {alpha[j2][0]*pow(Zeta,2), alpha[j2][1]*pow(Zeta,2), alpha[j2][2]*pow(Zeta,2)};

              for (int k=0; k < STOnG; k++) {
       
                     for (int l=0; l < atomn; l++) {
                            Zb[l] = atoms[l].atomic_number;
                            Rb[l][0] = atoms[l].x;
                            Rb[l][1] = atoms[l].y;
                            Rb[l][2] = atoms[l].z;

                            for (int l2=0; l2 < nlist[l]; l2++){

                                   d_vec_n = {D[l2][0],D[l2][1],D[l2][2]};
                                   alpha_vec_n = {alpha[l2][0]*pow(Zeta,2), alpha[l2][1]*pow(Zeta,2), alpha[l2][2]*pow(Zeta,2)};
       
                                   for (int q=0; q < STOnG; q++) {
                                          a = (j+1)*(j2+1)-1;
                                          b = (l+1)*(l2+1)-1;

                                          //std::cout << d_vec_m[k] << "\t" << d_vec_n[q] << std::endl;

                                          //Define gaussians for the overlap matrix.
                                          gauss GA;
                                          GA.alph=alpha_vec_m[k];
                                          GA.x=Ra[j][0];
                                          GA.y=Ra[j][1];
                                          GA.z=Ra[j][2];
                                          gauss GB;
                                          GB.alph=alpha_vec_n[q];
                                          GB.x=Rb[l][0];
                                          GB.y=Rb[l][1];
                                          GB.z=Rb[l][2];
                                          //std::cout << GA.alph << "\t" << GA.x << "\t" << GA.y << "\t" << GA.z << std::endl; 
                                          //Fill the overlap matrix.
                                          S[a][b] += d_vec_m[k]*d_vec_n[q]*overlap(GA,GB);
                                          //std::cout << "Overlap element S[" << a << "][" << b << "] is " << S[a][b] << std::endl;

                                          //Define the gaussians for the kinetic matrix.
                                          gauss KA;
                                          KA.alph=alpha_vec_m[k];
                                          KA.x=Ra[j][0];
                                          KA.y=Ra[j][1];
                                          KA.z=Ra[j][2];
                                          gauss KB;
                                          KB.alph=alpha_vec_n[q];
                                          KB.x=Rb[l][0];
                                          KB.y=Rb[l][1];
                                          KB.z=Rb[l][2];
                                          
                                          //Fill the kinetic matrix.
                                          T[a][b] = d_vec_m[k]*d_vec_n[q]*kinetic(KA,KB);
                                          //std::cout << "Kinetic element T[" << a << "][" << b << "] is " << T[a][b] << std::endl;
                                          //Define gaussians for the potential matrix.                                          
                                          for(int i5 = 0; i5 < atomn; i5++){
                                          gauss NA; 
                                          NA.alph = alpha_vec_m[k];
                                          NA.x = Ra[j][0];
                                          NA.y = Ra[j][1];
                                          NA.z = Ra[j][2];
                                          gauss NB; 
                                          NB.alph = alpha_vec_n[q];
                                          NB.x = Rb[l][0];
                                          NB.y = Rb[l][1];
                                          NB.z = Rb[l][2];

                                          //Populate the potential matrix.
                                          V[a][b] += d_vec_m[k]*d_vec_n[q]*potential(NA,NB,atoms[i5]);
                                          //std::cout << "Potential element V[" << a << "][" << b << "] is " <<  V[a][b] << std::endl;
                                          }

                                   }
                            }
                     }
              }
       }
}
for (int i6 = 0; i6 < 2; i6++){
       for(int j6 = 0; j6 < 2; j6++){
       std::cout << "Overlap element S[" << i6 << "][" << j6 << "] is " << S[i6][j6] << std::endl;
       std::cout << "Kinetic element T[" << i6 << "][" << j6 << "] is " << T[i6][j6] << std::endl;
       std::cout << "Potential element P[" << i6 << "][" << j6 << "] is " << V[i6][j6] << std::endl;

       }
}






}
