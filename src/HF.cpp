#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "HF.hpp"

#include "Eigen/Eigenvalues"

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

for(int s1 = 0; s1 < B; s1++){
       for(int s2 = 0; s2 < B; s2++){
              S[s1][s2]=0;
              T[s1][s2]=0;
              V[s1][s2]=0;
       }
}

double multi_electron_tensor[B][B][B][B];
for(int m1 = 0; m1 < B; m1++) {
       for(int m2 = 0; m2 < B; m2++){
              for(int m3=0; m3 < B; m3++){
                     for(int m4=0; m4 < B; m4++){
                            multi_electron_tensor[m1][m2][m3][m4] = 0;
                     }
              }
       }
}
//Iterate through atoms.

double Za[atomn];
double Ra[atomn][3];
std::vector<double> d_vec_m;
std::vector<double> alpha_vec_m;
double Zb[atomn];
double Rb[atomn][3];
std::vector<double> d_vec_n;
std::vector<double> alpha_vec_n;
double Zc[atomn];
double Rc[atomn][3];
std::vector<double> d_vec_k;
std::vector<double> alpha_vec_k;
double Zd[atomn];
double Rd[atomn][3];
std::vector<double> d_vec_l;
std::vector<double> alpha_vec_l;
int a=0;
int b=0;
int c=0;
int d=0;
for (int j=0; j < atomn; j++){
       
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
                                          T[a][b] += d_vec_m[k]*d_vec_n[q]*kinetic(KA,KB);
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
                                          for(int j3 = 0; j3 < atomn; j3++){
                                                 Zc[j3] = atoms[j3].atomic_number;
                                                 Rc[j3][0] = atoms[j3].x;
                                                 Rc[j3][0] = atoms[j3].y;
                                                 Rc[j3][0] = atoms[j3].z;
                                                 for (int j4 = 0; j4 < nlist[j3]; j4++){
                                                        d_vec_k = {D[j4][0],D[j4][1],D[j4][2]};
                                                        alpha_vec_k = {alpha[j4][0]*pow(Zeta,2), alpha[j4][1]*pow(Zeta,2), alpha[j4][2]*pow(Zeta,2)};
                                                        for(int k2=0; k2 < STOnG; k2++){
                                                               for (int i9 = 0; i9 < atomn; i9++){
                                                               Zd[k2]=atoms[k2].atomic_number;
                                                               Rd[k2][0]=atoms[i9].x;
                                                               Rd[k2][1]=atoms[i9].y;
                                                               Rd[k2][2]=atoms[i9].z;
                                                                      for(int l3=0; l3 < nlist[k2]; l3++){
                                                                             d_vec_l = {D[l3][0],D[l3][1],D[l3][2]};
                                                                             alpha_vec_l = {alpha[l3][0]*pow(Zeta,2), alpha[l3][1]*pow(Zeta,2), alpha[l3][2]*pow(Zeta,2)};
                                                                             for(int l4 = 0; l4 < STOnG; l4++){
                                                                                    c = (j3+1)*(j4+1)-1;
                                                                                    d = (k2+1)*(l3+1)-1;
                                                                                    //std::cout << a << "\t" << b << "\t" << c << "\t" << d << std::endl;
                                                                                    gauss MA;
                                                                                    MA.alph = alpha_vec_m[k];
                                                                                    MA.x = Ra[j][0];
                                                                                    MA.y = Ra[j][1];
                                                                                    MA.z = Ra[j][2];
                                                                                    gauss MB;
                                                                                    MB.alph = alpha_vec_n[q];
                                                                                    MB.x = Rb[l][0];
                                                                                    MB.y = Rb[l][1];
                                                                                    MB.z = Rb[l][2];
                                                                                    gauss MC;
                                                                                    MC.alph = alpha_vec_k[k2];
                                                                                    MC.x = Rc[i9][0];
                                                                                    MC.y = Rc[i9][1];
                                                                                    MC.z = Rc[i9][2];
                                                                                    gauss MD;
                                                                                    MD.alph = alpha_vec_l[l4];
                                                                                    MD.x = Rd[j3][0];
                                                                                    MD.y = Rd[j3][1];
                                                                                    MD.z = Rd[j3][2];
                                                                                    multi_electron_tensor[a][b][c][d]+=d_vec_m[k]*d_vec_n[q]*d_vec_k[k2]*d_vec_l[l4]*multi(MA,MB,MC,MD);
//                                                                                   std::cout << "Multi-Electron Tensor element Multi[" << a << "][" << b << "][" << c << "][" << d << "] is " << multi_electron_tensor[a][b][c][d] << std::endl;
             
                                                                             } 
                                                                      }
                                                               }          
                                                        }
                                                 }
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
for (int i7 = 0; i7 < 2; i7++){
       for (int j7 =0; j7 < 2; j7++){
              for (int i8 =0; i8 < 2; i8 ++){
                     for (int j8 =0; j8 < 2; j8++){
       std::cout << "Multi-Electron Tensor element Multi[" << i7 << "][" << j7 << "][" << i8 << "][" << j8 << "] is " << multi_electron_tensor[i7][j7][i8][j8] << std::endl;
                     }       
              } 
       }      
}
//Form the Core Hamiltonian.
double Hcore[B][B];
for(int hc1 =0; hc1 < B; hc1++){
       for(int hc2=0; hc2 < B; hc2++){
              Hcore[hc1][hc2]=0;
       }
}

for (int h1 = 0; h1 < B; h1++){
       for (int h2 = 0; h2 < B; h2++){
              Hcore[h1][h2] = T[h1][h2]+V[h1][h2];
              std::cout << "Core Hamiltonian element HCore[" << h1 << "][" << h2 << "] is " << Hcore[h1][h2] << std::endl;     
       }
}

//Perform the symmetric orthoganilization of basis.


}
