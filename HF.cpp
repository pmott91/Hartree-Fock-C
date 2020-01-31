#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>
using std::fixed;
using std::setprecision;

// Build a struct to store atomic number and coordinates.
struct Atom {
       int  atomic_number;
       double x, y, z;
       };

//Declare the function read_xyz which will collect coordinates 
//and number of atoms from an XYZ file.
std::vector<Atom> read_xyz(const std::string& xyz_file,int& atomn);

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
}

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

