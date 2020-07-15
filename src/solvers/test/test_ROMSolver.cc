/*
 * test_ROMSolver.cc
 *
 *  Created on: Jul 9, 2020
 *      Author: rabab
 */


#define TEST_LIST \
        FUNC(test_ROM_diffusion)\
		FUNC(test_ROM_EnergyIndependent)\
		FUNC(test_ROM_EnergyDependent)


#include "TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "utilities/MathUtilities.hh"
#include "projection_fixture.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "solvers/EigenvalueManager.hh"
#include "callow/vector/Vector.hh"
#include "solvers/rom/ROMSolver.hh"
#include "solvers/rom/ROMBasis.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using namespace std;

typedef callow::MatrixDense::SP_matrix SP_matrix;
typedef callow::Vector::SP_vector      SP_vector;



int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

int test_ROM_diffusion(int argc, char *argv[])
{
 Mesh1D::SP_mesh mesh = get_mesh();
 Material::SP_material mat = get_mat();
 InputDB::SP_input input = get_input();
 input->put<std::string>("bc_west",                    "reflect");
 input->put<std::string>("bc_east",                    "reflect");

 ROMSolver<_1D> ROM(input, mesh, mat, "diffusion");

 SP_matrix U;
 U = new callow::MatrixDense(40, 10);
 ROMBasis::GetBasis("/home/rabab/opt/detran/source/src/solvers/test/flux_basis_assem0_diff", U);
 SP_vector  ROM_flux;
 ROM_flux = new callow::Vector(40, 0.0);
 ROM.Solve(U, ROM_flux);

 // FOM solution
 input->put<std::string>("equation",  "diffusion");

 EigenvalueManager<_1D> manager(input, mat, mesh);
 manager.solve();
 int n = 20;
 callow::Vector phi0_fom(n, 0.0);
 callow::Vector phi1_fom(n, 0.0);

 callow::Vector phi0_rom(n, 0.0);
 callow::Vector phi1_rom(n, 0.0);

 for (int i = 0; i < n; ++i)
{
  phi0_fom[i] = manager.state()->phi(0)[i]/ detran_utilities::norm(manager.state()->phi(0), "L2");
  phi1_fom[i] = manager.state()->phi(1)[i]/detran_utilities::norm(manager.state()->phi(1), "L2");
  phi0_rom[i] = (*ROM_flux)[i];
  phi1_rom[i] = (*ROM_flux)[i + 20];
}
 vec_dbl error1 (n, 0);
 vec_dbl error2 (n, 0);

 for (int i = 0; i < 20; ++i)
 {
   error1[i] = phi0_fom[i]/phi0_fom.norm() - phi0_rom[i]/phi0_rom.norm();

   error2[i] = phi1_fom[i]/phi1_fom.norm() - phi1_rom[i]/phi1_rom.norm();
 }
 return 0;
}


int test_ROM_EnergyIndependent(int argc, char *argv[])
{
  Mesh1D::SP_mesh mesh = get_mesh();
  Material::SP_material mat = get_mat();
  InputDB::SP_input input = get_input();

  ROMSolver<_1D> ROM(input, mesh, mat, "EnergyIndependent");
  SP_matrix U;
  U = new callow::MatrixDense(20, 5);
  ROMBasis::GetBasis("/home/rabab/opt/detran/source/src/solvers/test/fission_density_basis_assem0", U);
  SP_vector  ROM_flux;
  int n=20;
  ROM_flux = new callow::Vector(n, 0.0);
  ROM.Solve(U, ROM_flux);

  EigenvalueManager<_1D> manager(input, mat, mesh);
  manager.solve();

  TEST(soft_equiv(manager.state()->eigenvalue(), ROM.keff(), 1E-4));

 return 0;
}

int test_ROM_EnergyDependent(int argc, char *argv[])
{
 Mesh1D::SP_mesh mesh = get_mesh();
 Material::SP_material mat = get_mat();
 InputDB::SP_input input = get_input();

 //input->put<std::string>("bc_west",                    "reflect");
 //input->put<std::string>("bc_east",                    "reflect");

 ROMSolver<_1D> ROM(input, mesh, mat, "EnergyDependent");

 SP_matrix U;

 int n=40;
 int r=10;
 U = new callow::MatrixDense(n, r);
 ROMBasis::GetBasis("/home/rabab/opt/detran/source/src/solvers/test/flux_basis_assem0_diff", U);
 SP_vector  ROM_flux;
 ROM_flux = new callow::Vector(n, 0.0);

 ROM.Solve(U, ROM_flux);

 return 0;
}


