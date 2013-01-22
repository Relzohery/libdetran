//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_DLP.cc
 *  @author Jeremy Roberts
 *  @date   Apr 1, 2012
 *  @brief  Test of DLP class
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_DLP)

#include "utilities/TestDriver.hh"
#include "utilities/Definitions.hh"
#include "orthog/DLP.hh"
#include "callow/utils/Initialization.hh"
#include <cmath>

using namespace detran_orthog;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

int test_DLP(int argc, char *argv[])
{

  int N = 2;
  int M = N-1;
  double width = 2.0 / N;

  vec_dbl x(N, 0);
  x[0] = -1.0 + width / 2.0;
  for (int i = 1; i < N; ++i) x[i] = x[i-1] + width;

  DLP P(M, N, true);

//  callow::Vector f(N, 0.0);
//  for (int i = 0; i < N; ++i) f[i] = std::cos(x[i]);
//
//  callow::Vector ft(M + 1, 0.0);
//  callow::Vector f2(N, 0.0);
//
//  P.transform(f, ft);
//  P.inverse(ft, f2);
//
//  f.display();
//  ft.display();
//  f2.display();
//  f2.subtract(f);
//  f2.display();
  P.basis()->display();
  P.coefficients()->display();
  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_MomentIndexer.cc
//---------------------------------------------------------------------------//
