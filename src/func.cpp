#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "HF.hpp"
//Define the function read_xyz.
std::vector<Atom> read_xyz(const std::string& xyz_file,int& atomn){
int j;

//Open the XYZ file.
std::ifstream infile;
infile.open(xyz_file);

std::cout << "Reading from file: " << xyz_file << "\n";

//Read in number of atoms from XYZ file.
infile >> atomn;

//Use 'atoms' to store the atomic number
//and coodinates of each atom.
std::vector<Atom> atoms(atomn);

       for(j = 0; j < atomn; j++) {
       std::string element_symbol;
       double x,y,z;

       infile >> element_symbol >> x >> y >> z;

              int Z;
              if(element_symbol =="H") 
              Z=1;
              else if(element_symbol == "He")
              Z=2;
              else if(element_symbol == "Li")
              Z=3;
              else if(element_symbol == "Be")
              Z=4;
              else if(element_symbol == "B")
              Z=5;
              else if(element_symbol == "C")
              Z=6;
              else if(element_symbol == "N")
              Z=7;
              else if(element_symbol == "O")
              Z=8;
              else if(element_symbol == "F")
              Z=9;
              else if(element_symbol == "Ne")
              Z=10;
              else if(element_symbol == "Na")
              Z=11;
              else if(element_symbol == "Mg")
              Z=12;
              else if(element_symbol == "Al")
              Z=13;
              else if(element_symbol == "Si")
              Z=14;
              else if(element_symbol == "P")
              Z=15;
              else if(element_symbol == "S")
              Z=16;
              else if(element_symbol == "Cl")
              Z=17;
              else if(element_symbol == "Ar")
              Z=18;
              else if(element_symbol == "K")
              Z=19;
              else if(element_symbol == "Ca")
              Z=20;
              else if(element_symbol == "Sc")
              Z=21;
              else if(element_symbol == "Ti")
              Z=22;
              else if(element_symbol == "V")
              Z=23;
              else if(element_symbol == "Cr")
              Z=24;
              else if(element_symbol == "Mn")
              Z=25;
              else if(element_symbol == "Fe")
              Z=26;
              else if(element_symbol == "Co")
              Z=27;
              else if(element_symbol == "Ni")
              Z=28;
              else if(element_symbol == "Cu")
              Z=29;
              else if(element_symbol == "Zn")
              Z=30;
              else if(element_symbol == "Ga")
              Z=31;
              else if(element_symbol == "Ge")
              Z=32;
              else if(element_symbol == "As")
              Z=33;
              else if(element_symbol == "Se")
              Z=34;
              else if(element_symbol == "Br")
              Z=35;
              else if(element_symbol == "Kr")
              Z=36;
              else if(element_symbol == "Rb")
              Z=37;
              else if(element_symbol == "Sr")
              Z=38;
              else if(element_symbol == "Y")
              Z=39;
              else if(element_symbol == "Zr")
              Z=40;
              else if(element_symbol == "Nb")
              Z=41;
              else if(element_symbol == "Mo")
              Z=42;
              else if(element_symbol == "Tc")
              Z=43;
              else if(element_symbol == "Ru")
              Z=44;
              else if(element_symbol == "Rh")
              Z=45;
              else if(element_symbol == "Pd")
              Z=46;
              else if(element_symbol == "Ag")
              Z=47;
              else if(element_symbol == "Cd")
              Z=48;
              else if(element_symbol == "In")
              Z=49;
              else if(element_symbol == "Sn")
              Z=50;
              else if(element_symbol == "Sb")
              Z=51;
              else if(element_symbol == "Te")
              Z=52;
              else if(element_symbol == "I")
              Z=53;
              else if(element_symbol == "Xe")
              Z=54;
              else if(element_symbol == "Cs")
              Z=55;
              else if(element_symbol == "Ba")
              Z=56;
              else if(element_symbol == "La")
              Z=57;
              else if(element_symbol == "Ce")
              Z=58;
              else if(element_symbol == "Pr")
              Z=59;
              else if(element_symbol == "Nd")
              Z=60;
              else if(element_symbol == "Pm")
              Z=61;
              else if(element_symbol == "Sm")
              Z=62;
              else if(element_symbol == "Eu")
              Z=63;
              else if(element_symbol == "Gd")
              Z=64;
              else if(element_symbol == "Tb")
              Z=65;
              else if(element_symbol == "Dy")
              Z=66;
              else if(element_symbol == "Ho")
              Z=67;
              else if(element_symbol == "Er")
              Z=68;
              else if(element_symbol == "Tm")
              Z=69;
              else if(element_symbol == "Yb")
              Z=70;
              else if(element_symbol == "Lu")
              Z=71;
              else if(element_symbol == "Hf")
              Z=72;
              else if(element_symbol == "Ta")
              Z=73;
              else if(element_symbol == "W")
              Z=74;
              else if(element_symbol == "Re")
              Z=75;
              else if(element_symbol == "Os")
              Z=76;
              else if(element_symbol == "Ir")
              Z=77;
              else if(element_symbol == "Pt")
              Z=78;
              else if(element_symbol == "Au")
              Z=79;
              else if(element_symbol == "Hg")
              Z=80;
              else if(element_symbol == "Tl")
              Z=81;
              else if(element_symbol == "Pb")
              Z=82;
              else if(element_symbol == "Bi")
              Z=83;
              else if(element_symbol == "Po")
              Z=84;
              else if(element_symbol == "At")
              Z=85;
              else if(element_symbol == "Rn")
              Z=86;
              else if(element_symbol == "Fr")
              Z=87;
              else if(element_symbol == "Ra")
              Z=88;
              else if(element_symbol == "Ac")
              Z=89;
              else if(element_symbol == "Th")
              Z=90;
              else if(element_symbol == "Pa")
              Z=91;
              else if(element_symbol == "U")
              Z=92;
              else if(element_symbol == "Np")
              Z=93;
              else if(element_symbol == "Pu")
              Z=94;
              else if(element_symbol == "Am")
              Z=95;
              else if(element_symbol == "Cm")
              Z=96;
              else if(element_symbol == "Bk")
              Z=97;
              else if(element_symbol == "Cf")
              Z=98;
              else if(element_symbol == "Es")
              Z=99;
              else if(element_symbol == "Fm")
              Z=100;
              else if(element_symbol == "Md")
              Z=101;
              else if(element_symbol == "No")
              Z=102;
              else if(element_symbol == "Lr")
              Z=103;
              else if(element_symbol == "Rf")
              Z=104;
              else if(element_symbol == "Db")
              Z=105;
              else if(element_symbol == "Sg")
              Z=106;
              else if(element_symbol == "Bh")
              Z=107;
              else if(element_symbol == "Hs")
              Z=108;
              else if(element_symbol == "Mt")
              Z=109;
              else if(element_symbol == "Ds")
              Z=110;
              else if(element_symbol == "Rg")
              Z=111;
              else if(element_symbol == "Cn")
              Z=112;
              else if(element_symbol == "Nh")
              Z=113;
              else if(element_symbol == "Fl")
              Z=114;
              else if(element_symbol == "Mc")
              Z=115;
              else if(element_symbol == "Lv")
              Z=116;
              else if(element_symbol == "Ts")
              Z=117;
              else if(element_symbol == "Og")
              Z=118;
              else{
              std::cout << "Invalid atomic symbol\n";
              }
       //Convert the input coordinates from angstroms to bohr.        
       const double ang_to_bohr = 1/0.52917721092;

       atoms[j].atomic_number = Z;
       atoms[j].x = x * ang_to_bohr;
       atoms[j].y = y * ang_to_bohr;
       atoms[j].z = z * ang_to_bohr;
}
//Close the XYZ file.
infile.close();
      
return atoms;
 
}



//Define the function gauss_product;
std::vector<double> gauss_product(gauss gauss_a, gauss gauss_b){

const double Pi = 3.14159265359;
const double eul=2.71828182846;
double a = gauss_a.alph;
std::vector<double> Ra = {gauss_a.x,gauss_a.y,gauss_a.z}; 
double b = gauss_b.alph;
std::vector<double> Rb = {gauss_b.x,gauss_b.y,gauss_b.z};
double p = a + b;
double diff = ((Rb[0]-Ra[0])*(Rb[0]-Ra[0])) + ((Rb[1]-Ra[1]) * (Rb[1]-Ra[1])) + ((Rb[2]-Ra[2]) * (Rb[2]-Ra[2]));
double N = pow((4*a*b/(Pi*Pi)), 0.75);
double kexp = -a*b/p*diff;
double K = N*pow(eul,kexp);
std::vector<double> Rp = {(a*Ra[0]+b*Rb[0])/p,(a*Ra[1]+b*Rb[1])/p,(a*Ra[2]+b*Rb[2])/p};
std::vector<double> gauss_c = {p, diff, K, Rp[0], Rp[1], Rp[2]};
return gauss_c;
}


//Define the overlap integral.

double overlap(gauss GA, gauss GB){
std::vector<double> volve = gauss_product(GA,GB);
const double Pi = 3.14159265359;
double p = volve[0];
double diff = volve[1];
double K = volve[2];
double prefactor = pow((Pi/p),1.5);
double ans = prefactor*K;
return ans;
}

//Define the kinetic integral.

double kinetic(gauss KA, gauss KB){

const double Pi = 3.14159265359;
std::vector<double> volve = gauss_product(KA,KB);
double p =volve[0];
double diff = volve[1];
double K = volve[2];
std::vector<double> Rp = {volve[3], volve[4], volve[5]};
double prefactor = pow((Pi/p),1.5);
double a = KA.alph;
std::vector<double> Ra = {KA.x, KA.y, KA.z};
double b = KB.alph;
std::vector<double> Rb = {KB.x, KB.y, KB.z};
double reduced_exponent = (a*b)/p;
double ans = reduced_exponent*(3-2*reduced_exponent*diff)*prefactor*K;
return ans;
}


//Define the F0 function.

double fo(double t){
const double Pi = 3.14159265359;
double ans;
if(t ==0) {
return 1;
}
else{
ans=0.5*pow((Pi/t),0.5)*erf(pow(t,0.5));
return ans;
}
}

//Define the Nuclear-electron integral.

double potential(gauss NA, gauss NB, Atom CAtom){
const double Pi =3.14159265359;
std::vector<double> volve = gauss_product(NA,NB);
Atom natom={CAtom.atomic_number, CAtom.x, CAtom.y, CAtom.z};
double p =volve[0];
double diff = volve[1];
double K = volve[2];
std::vector<double> Rp = {volve[3], volve[4], volve[5]};
std::vector<double> Rc = {natom.x, natom.y, natom.z};
int Zc = natom.atomic_number;
double norm = (pow((Rp[0]-Rc[0]),2) + pow((Rp[1]-Rc[1]),2) + pow((Rp[2]-Rc[2]),2));
double ans = (-2*Pi*Zc/p)*K*fo(p*norm);
return ans;
}

//Define the (ab|cd) integral.

double multi(gauss MA, gauss MB, gauss MC, gauss MD){
const double Pi=3.14159265359;
std::vector<double> volve = gauss_product(MA,MB);
double p =volve[0];
double diff_ab = volve[1];
double K_ab = volve[2];
std::vector<double> Rp = {volve[3], volve[4], volve[5]};
std::vector<double> bolve = gauss_product(MC,MD);
double q =bolve[0];
double diff_cd = bolve[1];
double K_cd = bolve[2];
std::vector<double> Rq = {bolve[3], bolve[4], bolve[5]};
double multi_prefactor = 2*pow(Pi,2.5)*pow((p*q*pow((p+q),0.5)),-1);
double norm = (pow((Rp[0]-Rq[0]),2) + pow((Rp[1]-Rq[1]),2) + pow((Rp[2]-Rq[2]),2));
double ans = multi_prefactor*K_ab*K_cd*fo(p*q/(p+q)*norm);
return ans;
}

//Define the succesive density matrix elements.

double succ_dens(Eigen::MatrixXd ptilde, Eigen::MatrixXd p, int B){
double x=0;
for(int xitr=0; xitr < B; xitr++){
       for(int xitr2=0; xitr2 < B; xitr2++){
              x+=pow(B,2)*pow((ptilde(xitr,xitr2)-p(xitr,xitr2)),2);
       }
}
return pow(x,0.5);
}
