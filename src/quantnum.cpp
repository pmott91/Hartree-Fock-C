#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "HF.hpp"

std::vector<int> max_quantum_number(std::vector<int>& atomz, int& atomn){

std::vector<int> nlist;

int i;
for (i = 0; i < atomn; i++){
              int qnum;
              if(atomz[i] == 1) 
              qnum=1;
              else if(atomz[i] == 2)
              qnum=2;
              else if(atomz[i] == 3)
              qnum=2;
              else if(atomz[i] == 4)
              qnum=2;
              else if(atomz[i] == 5)
              qnum=2;
              else if(atomz[i] == 6)
              qnum=2;
              else{
              std::cout << "Functionality for this element isn't available yet.\n";
              }
              nlist.push_back(qnum);
}
return nlist;
}
