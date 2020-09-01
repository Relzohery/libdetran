/*
 * test_c5g7.cc
 *
 *  Created on: Jul 28, 2020
 *      Author: rabab
 */




#define TEST_LIST \
        FUNC(test_pin)\
		FUNC(test_assembly)\
        FUNC(test_core)



#include "TestDriver.hh"
#include "utilities/MathUtilities.hh"

#include "callow/utils/Initialization.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/vector/Vector.hh"
#include "solvers/EigenvalueManager.hh"
#include "solvers/rom/ROMSolver.hh"
#include "solvers/rom/ROMBasis.hh"
#include "c5g7_pins.hh"

using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

typedef callow::MatrixDense::SP_matrix SP_matrix;
typedef callow::Vector::SP_vector      SP_vector;



int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

int test_pin(int argc, char *argv[])
{
  Material::SP_material mat = get_mat();

  vec_dbl radii (1, 0.54);
  vec_int  mat_map(2, 0);
  mat_map[0] = 1;
  mat_map[1] = 6 ;
  SP_pincell pin = get_pin(radii, mat_map);
  SP_mesh  mesh = pin->mesh();
  mat->display();

  InputDB::SP_input inp = get_input();
  EigenvalueManager<_2D> manager(inp, mat, mesh);
  manager.solve();
  double keff = manager.state()->eigenvalue();

  return 0;
}

int test_assembly(int argc, char *argv[])
{
  Material::SP_material mat = get_mat();

  Assembly::SP_assembly assm = get_assembly(0);

  SP_mesh  mesh = assm->mesh();
  int n = mesh->number_cells();

  InputDB::SP_input inp = get_input("assembly");

  std::cout << "########################################################" << "\n";
  std::cout << "################### Full Order Model  ##################" << "\n";
  std::cout << "########################################################" << "\n";

  EigenvalueManager<_2D> manager(inp, mat, mesh);
  manager.solve();
  double keff = manager.state()->eigenvalue();

  inp->put<std::string>("equation", "dd");


  std::cout << "###########################################################" << "\n";
  std::cout << "################### Reduced Order Model  ##################" << "\n";
  std::cout << "###########################################################" << "\n";
  ROMSolver<_2D> ROM(inp, mesh, mat);
  int r = 50;

  // get the basis
  SP_matrix U;
  U = new callow::MatrixDense(2601, 20);
  ROMBasis::GetBasis("/home/rabab/opt/detran/source/src/solvers/test/fission_density_assem0_transport_c5g7_r=20", U);

  SP_vector  fd_rom;
  // ROM
  fd_rom = new callow::Vector(n, 0.0);
  ROM.Solve(U, fd_rom);

 double keff_rom = ROM.keff();
 return 0;
}


int test_core(int argc, char *argv[])
{
  Material::SP_material mat = get_mat();
  Core::SP_core core = get_core();
  SP_mesh  mesh = core->mesh();

  InputDB::SP_input inp = get_input("core");

  /*
  std::cout << "########################################################" << "\n";
  std::cout << "################### Full Order Model  ##################" << "\n";
  std::cout << "########################################################" << "\n";

  EigenvalueManager<_2D> manager(inp, mat, mesh);
  manager.solve();
  double keff = manager.state()->eigenvalue();
  std::cout <<  "keff = " << keff << "\n";
  */

  std::cout << "###########################################################" << "\n";
  std::cout << "################### Reduced Order Model  ##################" << "\n";
  std::cout << "###########################################################" << "\n";

  ROMSolver<_2D> ROM(inp, mesh, mat);
  std::cout << mesh->number_cells() << "\n";
  int r = 20;
  int n = 2601*7;
  // get the basis
  SP_matrix U;
  U = new callow::MatrixDense(n, 20);
  ROMBasis::GetBasis("/home/rabab/opt/detran/source/src/solvers/test/flux_basis_core_transport_c5g7_r=20", U);

  SP_vector  fd_rom;
  // ROM
  fd_rom = new callow::Vector(n, 0.0);
  ROM.Solve(U, fd_rom);

  double keff_rom = ROM.keff();
  std::cout <<  "keff rom = " << keff_rom << "\n";

return 0;
}

