#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "HF.hpp"

#include "Eigen/Eigenvalues"
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"

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
const std::vector<Atom> atoms = read_xyz(xyz_file,atomn);

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
B+=nlist[i4];
}

std::cout << std::endl << "Basis set size is " << B << std::endl << std::endl;

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
       
       enuc += (atoms[k].atomic_number * atoms[l].atomic_number)/r;

       }
}

std::cout << "The nuclear repulsion energy is: " << enuc << std::endl << std::endl;

//Initialize matrices

Eigen::MatrixXd S(B,B);
Eigen::MatrixXd T(B,B);
Eigen::MatrixXd V(B,B);

for(int s1 = 0; s1 < B; s1++){
       for(int s2 = 0; s2 < B; s2++){
              S(s1,s2)=0;
              T(s1,s2)=0;
              V(s1,s2)=0;
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
                                          //Fill the overlap matrix.
                                          S(a,b) += d_vec_m[k]*d_vec_n[q]*overlap(GA,GB);
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
                                          T(a,b) += d_vec_m[k]*d_vec_n[q]*kinetic(KA,KB);
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
                                                 V(a,b) += d_vec_m[k]*d_vec_n[q]*potential(NA,NB, atoms[i5]);
                                          }
                                          for(int j3 = 0; j3 < atomn; j3++){
                                                 Zc[j3] = atoms[j3].atomic_number;
                                                 Rc[j3][0] = atoms[j3].x;
                                                 Rc[j3][1] = atoms[j3].y;
                                                 Rc[j3][2] = atoms[j3].z;
                                                  
                                                 for (int j4 = 0; j4 < nlist[j3]; j4++){
                                                        d_vec_k = {D[j4][0],D[j4][1],D[j4][2]};
                                                        alpha_vec_k = {alpha[j4][0]*pow(Zeta,2), alpha[j4][1]*pow(Zeta,2), alpha[j4][2]*pow(Zeta,2)};
                                                        for(int k2=0; k2 < STOnG; k2++){
                                                               for (int i9 = 0; i9 < atomn; i9++){
                                                               Zd[i9]=atoms[i9].atomic_number;
                                                               Rd[i9][0]=atoms[i9].x;
                                                               Rd[i9][1]=atoms[i9].y;
                                                               Rd[i9][2]=atoms[i9].z;
                                                                      for(int l3=0; l3 < nlist[i9]; l3++){
                                                                             d_vec_l = {D[l3][0],D[l3][1],D[l3][2]};
                                                                             alpha_vec_l = {alpha[l3][0]*pow(Zeta,2), alpha[l3][1]*pow(Zeta,2), alpha[l3][2]*pow(Zeta,2)};
                                                                             for(int l4 = 0; l4 < STOnG; l4++){
                                                                                    c = (j3+1)*(j4+1)-1;
                                                                                    d = (i9+1)*(l3+1)-1;
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
                                                                                    MC.x = Rc[j3][0];
                                                                                    MC.y = Rc[j3][1];
                                                                                    MC.z = Rc[j3][2];
                                                                                    gauss MD;
                                                                                    MD.alph = alpha_vec_l[l4];
                                                                                    MD.x = Rd[i9][0];
                                                                                    MD.y = Rd[i9][1];
                                                                                    MD.z = Rd[i9][2];
                                                                                    multi_electron_tensor[a][b][c][d]+=d_vec_m[k]*d_vec_n[q]*d_vec_k[k2]*d_vec_l[l4]*multi(MA,MB,MC,MD);
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
std::cout << "Here is the Overlap Matrix:" << std::endl << S << std::endl << std::endl;
std::cout << "Here is the Kinetic Matrix:" << std::endl << T << std::endl << std::endl;
std::cout << "Here is the Potential Matrix:" << std::endl << V << std::endl << std::endl;
for (int i7 = 0; i7 < 2; i7++){
       for (int j7 =0; j7 < 2; j7++){
              for (int i8 =0; i8 < 2; i8 ++){
                     for (int j8 =0; j8 < 2; j8++){
       std::cout << "Multi-Electron Tensor element Multi[" << i7 << "][" << j7 << "][" << i8 << "][" << j8 << "] is " << multi_electron_tensor[i7][j7][i8][j8] << std::endl;
                     }       
              } 
       }      
}
std::cout << std::endl;
//Form the Core Hamiltonian.
Eigen::MatrixXd Hcore(B,B);
for(int hc1 =0; hc1 < B; hc1++){
       for(int hc2=0; hc2 < B; hc2++){
              Hcore(hc1,hc2)=0;
       }
}

Hcore=T+V;

std::cout << "Core Hamiltonian is:" << std::endl << Hcore << std::endl << std::endl;

std::cout << std::endl;
//Perform the symmetric orthoganilization of basis.

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
//std::cout << "The eigenvalues of S are:" << std::endl << es.eigenvalues() << std::endl << std::endl;
//std::cout << "The eigenvectors of S are:" << std::endl << es.eigenvectors() << std::endl << std::endl;
auto U = es.eigenvectors();
//std::cout << "U matrix is:" << std::endl << U << std::endl << std::endl;
auto Utran = U.adjoint();
//std::cout << "U transpose is:" << std::endl << Utran << std::endl << std::endl;

auto USU = Utran*S*U;
//std::cout << "The diagonalized overlap matrix, s, is:" << std::endl << USU << std::endl << std::endl;

auto sinv = USU.inverse();
//std::cout << "The inverse of this matrix is:" << std::endl << sinv << std::endl << std::endl;


auto sinvsq = sinv.sqrt();
auto Xmat = U*sinvsq*Utran;
//std::cout << "The X matrix is:" << std::endl << Xmat << std::endl << std::endl;

//Initial guess at P.

Eigen::MatrixXd P(B,B);
for(int piter =0; piter < B; piter++){
       for(int piter2=0; piter2 < B; piter2++){
              P(piter,piter2)=0;
       }
}
Eigen::MatrixXd Pprevious(B,B);
for(int priter =0; priter < B; priter++){
       for(int priter2=0; priter2 < B; priter2++){
              P(priter,priter2)=0;
       }
}

double thresh = 100;
int ccount =0;
Eigen::MatrixXd evalfprime;
Eigen::MatrixXd Cmat;

while(thresh > .0001){

//Calculate the Fock matrix with a guess.

Eigen::MatrixXd G(B,B);

       for(int gitr1=0; gitr1 < B; gitr1++){
              for(int gitr2=0; gitr2 < B; gitr2++){
                     G(gitr1,gitr2)=0;
              }
       }  

       for(int gitra = 0; gitra < B; gitra++){
              for(int gitrb=0; gitrb < B; gitrb++){
                     for(int gitrc=0; gitrc < B; gitrc++){
                            for(int gitrd =0; gitrd < B; gitrd++){
                                   G(gitra,gitrb) += P(gitrc,gitrd)*(multi_electron_tensor[gitra][gitrb][gitrc][gitrd] - 0.5*multi_electron_tensor[gitra][gitrc][gitrd][gitrb]);
                            }
                     }
              }
       }

auto Fock = Hcore+G;

std::cout << "Fock Matrix is:" << std::endl << Fock << std::endl << std::endl;

Eigen::MatrixXd Fockprime = Xmat.transpose()*Fock*Xmat;

std::cout << "The Fock' matrix is:" << std::endl << Fockprime << std::endl << std::endl;

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ef(Fockprime);

//std::cout << "The eigenvectors of the Fock' matrix are" << std::endl << std::endl << ef.eigenvectors() << std::endl;
//std::cout << "The eigenvalues of the Fock' matrix are" << std::endl << std::endl << ef.eigenvalues() << std::endl;

evalfprime = ef.eigenvalues();
Eigen::MatrixXd cprime = ef.eigenvectors();

Cmat = Xmat * cprime;

//std::cout << "The C matrix is: " << std::endl << Cmat << std::endl << std::endl; 

for(int npitr1 = 0; npitr1 < B; npitr1++){
       for(int npitr2 = 0; npitr2 < B; npitr2++){
              for(int npitr3 = 0; npitr3 < (N/2); npitr3++){
                     P(npitr1,npitr2)=2*Cmat(npitr1,npitr3)*Cmat(npitr2,npitr3);
              }
       }
}

//std::cout << "The P matrix is:" << std::endl << P << std::endl << std::endl;

thresh = succ_dens(P,Pprevious,B);

Pprevious = P;

ccount+=1;

}

std::cout << "The HF algorithm took " << ccount << " iterations to converge" << std::endl << std::endl;

std::cout << "The orbital energies are " << std::endl << evalfprime << std::endl << std::endl;

std::cout << "The orbital matrix is: " << std::endl << Cmat << std::endl << std::endl;

std::cout << "The density/bond order matrix is: " << std::endl << P << std::endl << std::endl;

}
