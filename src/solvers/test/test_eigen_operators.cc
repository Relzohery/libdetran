/*
 * eigen_operators.cc
 *
 *  Created on: May 28, 2020
 *      Author: rabab
 */

// LIST OF TEST FUNCTIONS
#define TEST_LIST                              \
        FUNC(test_eigen_operators)                  \

#include "TestDriver.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/test/eigenvalue_fixture.hh"
#include "utilities/MathUtilities.hh"
#include "callow/solver/EigenSolver.hh"
#include "solvers/eigen/EigenArnoldi.hh"
#include "callow/solver/EigenSolverCreator.hh"
#include "EigenArnoldi.hh"
#include "callow/solver/EigenSolverCreator.hh"
#include "utilities/InputDB.hh"
#include <string>
#include "FixedSourceManager.hh"
#include "solvers/eigen/EnergyDependentEigenLHS.hh"



using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

int test_eigen_operators(int argc, char *argv[])
{

EigenvalueData data = get_eigenvalue_data(1, 2);

data.input->put<std::string>("outer_solver", "GMRES");
data.input->put<std::string>("eigen_solver", "GD");
typedef typename Eigensolver<_1D>::Fixed_T              Fixed_T;
typedef typename Eigensolver<_1D>::SP_solver            SP_solver;
typedef typename Fixed_T::SP_manager                  SP_mg_solver;
typedef EnergyIndependentEigenOperator<_1D>         Operator_T;
typedef typename Operator_T::SP_operator          SP_operator;
typedef EnergyDependentEigenLHS<_1D>                LHS_Operator_T;
typename LHS_Operator_T::SP_operator B;

SP_mg_solver mg_solver;
mg_solver = new FixedSourceManager<_1D>(data.input, data.material, data.mesh, true, true);
mg_solver->setup();
mg_solver->set_solver();

SP_operator A;
A = new Operator_T(mg_solver);
B = new LHS_Operator_T(mg_solver);

//std::cout << B->number_columns() << "\n";
//std::cout << B->number_rows() << "\n";

//std::cout << A->number_columns() << "\n";
//std::cout << A->number_rows() << "\n";


int m = A->number_columns();
int n = A->number_rows();

int m1 = B->number_columns();
int n1 = B->number_rows();

callow::Vector x(m, 1.0);
callow::Vector y(n, 0.0);
callow::Vector x1(m1, 1.0);
callow::Vector y1(n1, 0.0);

for (int i=0; i< m; i++)
{
  x[i] = 0.5*i;
}

for (int i=0; i< m1; i++)
{
x1[i] = 0.5*i;
}

std::cout << "matrix one" << "\n";

A->multiply(x, y);

std::cout << "matrix two" << "\n";

B->multiply(x1, y1);

for (int i=0; i < n ; i++)
{
std::cout << y[i] << "\n";
}

std::cout << "************************" << "\n";

for (int i=0; i < n1 ; i++)
{
std::cout << y1[i] << "\n";
}


return 0;
}

