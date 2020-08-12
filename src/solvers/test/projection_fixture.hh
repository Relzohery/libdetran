
/*
 * projection_fixture.hh
 *
 *  Created on: Jul 7, 2020
 *      Author: rabab
 */


#include "Mesh1D.hh"
#include "Material.hh"
#include "callow/vector/Vector.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/vector/Vector.hh"
#include "utilities/InputDB.hh"
#include "Definitions.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace detran_material;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace std;
using std::cout;
using std::endl;


Material::SP_material get_mat()
{

  Material::SP_material mat = Material::Create(4, 2, "slabreactor");
  // Material 0: Water
  // Total
  mat->set_sigma_t(0, 0, 0.185343257302);
  mat->set_sigma_t(0, 1, 1.45564669643);

  // Fission
  mat->set_sigma_f(0, 0, 0);
  mat->set_sigma_f(0, 1, 0);
  mat->set_chi(0, 0, 0.0);
  mat->set_chi(0, 1, 0.0);

  // Scattering
  mat->set_sigma_s(0, 0, 0, 0.152722299219);
  mat->set_sigma_s(0, 0, 1, 0.000);
  mat->set_sigma_s(0, 1, 0, 0.0160150759918);
  mat->set_sigma_s(0, 1, 1, 1.43551589635);
  mat->compute_sigma_a();
  mat->compute_diff_coef();


  // Material 1: Fuel I
  // Total
  mat->set_sigma_t(1, 0, 0.229804263342);
  mat->set_sigma_t(1, 1, 1.02891192068);

  // Fission
  mat->set_sigma_f(1, 0, 0.0067060795922);
  mat->set_sigma_f(1, 1, 0.125201394293);
  mat->set_chi(1, 0, 1.0);
  mat->set_chi(1, 1, 0.0);

  // Scattering
  mat->set_sigma_s(1, 0, 0, 0.202916148126);
  mat->set_sigma_s(1, 0, 1, 0.0000);
  mat->set_sigma_s(1, 1, 0, 0.0163393070276);
  mat->set_sigma_s(1, 1, 1, 0.924336913278);
  mat->compute_sigma_a();
  mat->compute_diff_coef();                   ;


  // Material 3: Fuel II

  // Total
  mat->set_sigma_t(2, 0, 0.224822827113);
  mat->set_sigma_t(2, 1, 0.998816575222);

  // Fission
  mat->set_sigma_f(2, 0, 0.00768317162438);
  mat->set_sigma_f(2, 1, 0.125991676867);
  mat->set_chi(2, 0, 1.0);
  mat->set_chi(2, 1, 0.0);

  // Scattering
  mat->set_sigma_s(2, 0, 0, 0.199026443318);
  mat->set_sigma_s(2, 0, 1, 0.0000);
  mat->set_sigma_s(2, 1, 0, 0.0158835171486);
  mat->set_sigma_s(2, 1, 1, 0.887589148437);
  mat->compute_sigma_a();
  mat->compute_diff_coef();                   ;

  // Material 4: Fuel II + Gd

 // Total
  mat->set_sigma_t(3, 0, 0.217344332576);
  mat->set_sigma_t(3, 1, 1.05429714894);

  // Fission
  mat->set_sigma_f(3, 0, 0.00549089320936);
  mat->set_sigma_f(3, 1, 0.0178328973774);
  mat->set_chi(3, 0, 1.0);
  mat->set_chi(3, 1, 0.0);

   // Scattering
  mat->set_sigma_s(3, 0, 0, 0.190555374437);
  mat->set_sigma_s(3, 0, 1, 0.0000);
  mat->set_sigma_s(3, 1, 0, 0.013504085703);
  mat->set_sigma_s(3, 1, 1, 0.562862095564);
  mat->compute_sigma_a();
  mat->compute_diff_coef();                   ;

  mat->finalize();

  return mat;
}

Mesh1D::SP_mesh get_mesh(int fmm = 1, std::string id="assembly")
{
  vec_dbl cm_assembly(7, 0.0);
  cm_assembly[0] = 0.0;
  cm_assembly[1] = 1.1580;
  cm_assembly[2] = 4.4790;
  cm_assembly[3] = 7.8000;
  cm_assembly[4] = 11.1210;
  cm_assembly[5] = 14.4420;
  cm_assembly[6] = 15.6000;

  vec_int fm_assembly(6, 0);

  fm_assembly[0] = 2;
  fm_assembly[1] = 4;
  fm_assembly[2] = 4;
  fm_assembly[3] = 4;
  fm_assembly[4] = 4;
  fm_assembly[5] = 2;

  // a core is composed of 7 adjacent assemblies
  vec_dbl cm_core (43, 0);
  vec_int fm_core(7*6);

  for (int i=0; i <7; i++)
   for (int j=0; j<6; j++)
	 {

      cm_core[j+(6*i) + 1] = cm_assembly[j +1] + 15.6 *i;

	  fm_core[j+(6*i)] = fm_assembly[j];
	 }
  cm_core [0] = 0.0;

  //cm_core.assign(0.0 ,0.0);

  Mesh1D::vec_int mt_0(6);
  mt_0[0] = 0; mt_0[1] = 1;
  mt_0[2] = 2; mt_0[3] = 2;
  mt_0[4] = 1; mt_0[5] = 0;

  Mesh1D::vec_int mt_1(6);
  mt_1[0] = 0; mt_1[1] =1;
  mt_1[2] = 1; mt_1[3] = 1;
  mt_1[4] = 1; mt_1[5] = 0;

  Mesh1D::vec_int mt_2(6);
  mt_2[0] = 0; mt_2[1] =1;
  mt_2[2] = 3; mt_2[3] = 3;
  mt_2[4] = 1; mt_2[5] = 0;

  Mesh1D::vec_int mt_3(6);
  mt_3[0] = 0; mt_3[1] = 3;
  mt_3[2] = 3; mt_3[3] = 3;
  mt_3[4] = 3; mt_3[5] = 0;


  Mesh1D::vec_int mt_core0(42);
  mt_core0[0] = 0; mt_core0[1] = 1;
  mt_core0[2] = 2; mt_core0[3] = 2;
  mt_core0[4] = 1; mt_core0[5] = 0;

  mt_core0[6] = 0; mt_core0[7] = 1;
  mt_core0[8] = 1; mt_core0[9] = 1;
  mt_core0[10] = 1; mt_core0[11] = 0;

  mt_core0[12] = 0; mt_core0[13] = 1;
  mt_core0[14] = 2; mt_core0[15] = 2;
  mt_core0[16] = 1; mt_core0[17] = 0;

  mt_core0[18] = 0; mt_core0[19] = 1;
  mt_core0[20] = 1; mt_core0[21] = 1;
  mt_core0[22] = 1; mt_core0[23] = 0;

  mt_core0[24] = 0; mt_core0[25] = 1;
  mt_core0[26] = 2; mt_core0[27] = 2;
  mt_core0[28] = 1; mt_core0[29] = 0;

  mt_core0[30] = 0; mt_core0[31] = 1;
  mt_core0[32] = 1; mt_core0[33] = 1;
  mt_core0[34] = 1; mt_core0[35] = 0;

  mt_core0[36] = 0; mt_core0[37] = 1;
  mt_core0[38] = 2; mt_core0[39] = 2;
  mt_core0[40] = 1; mt_core0[41] = 0;

  if (id == "assembly")
  {
  Mesh1D::SP_mesh mesh = Mesh1D::Create(fm_assembly, cm_assembly, mt_0);
  return mesh;
  }
  else if (id == "core")
  {
   Mesh1D::SP_mesh mesh = Mesh1D::Create(fm_core, cm_core, mt_core0);

   return mesh;
  }

}

InputDB::SP_input get_input()
{
  InputDB::SP_input inp(new InputDB("Slab Reactor"));
  inp->put<std::string>("problem_type", "eigenvalue");
  inp->put<int>("number_groups",                  2);
  inp->put<std::string>("inner_solver",           "SI");
  inp->put<int>("inner_max_iters",            1000);
  inp->put<int>("inner_max_iters",            1000);
  inp->put<double>("inner_tolerance",            1e-7);
  inp->put<int>("inner_print_level",          0);
  inp->put<std::string>("outer_solver",              "GS");
  inp->put<int>("outer_max_iters",            1000);
  inp->put<double>("outer_tolerance",            1e-7);
  inp->put<int>("inner_print_level",          0);
  inp->put<int>("outer_print_level",          1);
 // inp->put<std::string>("eigen_solver",       "arnoldi");
  inp->put<int>("eigen_max_iters",            200);
  inp->put<double>("eigen_tolerance",            1e-7);
  inp->put<std::string>("bc_west",                    "vacuum");
  inp->put<std::string>("bc_east",                    "vacuum");
  inp->put<int>("quad_number_polar_octant",   16);

  InputDB::SP_input db(new InputDB("callow dp"));
  db->put<double>("linear_solver_atol",              1e-12);
  db->put<double>("linear_solver_rtol",              1e-12);
  db->put<std::string>("linear_solver_type",              "petsc");
  db->put<int>("linear_solver_maxit",             5000);
  db->put<int>("linear_solver_gmres_restart",     30);
  db->put<int>("linear_solver_monitor_level",     0);
  db->put<std::string>("pc_type",                         "petsc_pc");
  db->put<std::string>("petsc_pc_type",                   "lu");
  db->put<std::string>("eigen_solver_type",               "slepc");
  db->put<int>("eigen_solver_monitor_level",      2);
  db->put<InputDB::SP_input>("inner_solver_db",               db);
  db->put<InputDB::SP_input>("inner_pc_db",                   db);
  db->put<InputDB::SP_input>("outer_solver_db",               db);
  db->put<InputDB::SP_input>("eigen_solver_db",               db);


  return inp;
}


