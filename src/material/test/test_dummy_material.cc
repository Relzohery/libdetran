/*
 * test_dummy_material.cc
 *
 *  Created on: Jul 25, 2020
 *      Author: rabab
 */



#include "dummy_material.hh"

#define TEST_LIST                     \
        FUNC(test_dummy_material)     \

// Detran headers
#include "TestDriver.hh"
#include "Material.hh"
#include "SoftEquivalence.hh"

// System
#include <iostream>

// Setup

using namespace detran_test;
using namespace detran_material;
using namespace detran_utilities;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_dummy_material(int argc, char *argv[])
{

  std::cout << "begin testing .............." << "\n";
  dummy_material mat(1);
  std::cout << "Done testing .............." << "\n";


  return 0;
}
