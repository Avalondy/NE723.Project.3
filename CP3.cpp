/*
Main program of project 3.
Will need to include the header file CP3.h.
Author: Yao Du
Email: ydu9@ncsu.edu
Date: Dec 3, 2019
*/
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>
#include "CP3.h"

int main(int argc, char const *argv[]) {
  std::string args(argv[argc-1]);
  TestClass which_test;
  if (args==("A")) which_test = TestClass::TestA;
  else if(args=="B") which_test = TestClass::TestB;
  else {
    std::cout << "Error! Specify which test to run." << '\n';
    return 0;
  }
  // TestClass which_test = TestClass::TestA; //specify which test to run
  TransportEquation TE_obj(which_test);
  TE_obj.readInputs();
  TE_obj.initializeParameters();
  TE_obj.multiGroupSolve();
  TE_obj.printResults();
  return 0;
}
