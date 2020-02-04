#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "HF.hpp"

int main () {
gauss guess_a,guess_b,guess_c,guess_d;
guess_a.alph=2.5;
guess_a.x=3.5;
guess_a.y=4.5;
guess_a.z=5.5;
guess_b.alph=2.0;
guess_b.x=3.0;
guess_b.y=4.0;
guess_b.z=5.0;
guess_c.alph=6.5;
guess_c.x=7.5;
guess_c.y=8.5;
guess_c.z=9.5;
guess_d.alph=6.0;
guess_d.x=7.0;
guess_d.y=8.0;
guess_d.z=9.0;

double mul_out = multi(guess_a,guess_b,guess_c,guess_d);

std::cout << mul_out << std::endl;

//double t=0;
//double fout = fo(t);
//std::cout << fout;
}
