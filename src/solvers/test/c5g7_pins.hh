/*
 * c5g7_pins.hh
 *
 *  Created on: Jul 27, 2020
 *      Author: rabab
 */

#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"
#include "material/Material.hh"
#include "utilities/InputDB.hh"
#include "material/test/material_fixture.hh"
#include "geometry/RegionFactory.hh"
#include "callow/vector/Vector.hh"
#include "utilities/InputDB.hh"
#include "PinCell.hh"
#include "Assembly.hh"
#include "Core.hh"
#include "c5g7_material.hh"



using namespace detran;
using namespace detran_test;
using namespace detran_utilities;
using namespace std;

using namespace detran_material;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace callow;
using namespace std;
using std::cout;
using std::endl;

typedef Mesh::SP_mesh SP_mesh;
typedef PinCell::SP_pincell SP_pincell;
typedef PinCell::vec_dbl vec_dbl;
typedef PinCell::vec_int vec_int;


//typedef detran_geometry::Assembly::SP_assembly SP_assembly;
//typedef detran_geometry::Core::SP_core SP_core;

 SP_pincell get_pin(vec_dbl radii, vec_int  mat_map)
{
  // PinCell(double pitch, vec_dbl radii,
  //         vec_int mat_map, bool fuel_flag = true);
  Point pitch(1.26, 1.26);
  bool fuel_flag = true; // this is a fuel pin (i.e. not a moderator box)
//  PinCell(const Point      &pitch,
//          const vec_int    &mat_map,
//          const vec_dbl    &radii = vec_dbl(0),
//          const size_t      division = DIVISION_NONE,
//          const Point      &pincenter = Point(0));
  SP_pincell pin(new detran_geometry::PinCell(pitch, mat_map, radii));
  pin->meshify(3, true);
  return pin;
}

 InputDB::SP_input get_input(std::string id = "pin")
 {
   InputDB::SP_input inp(new InputDB("c5g7 pin"));
   //inp->put<std::string>("equation",                       "dd");
   inp->put<std::string>("problem_type",  "eigenvalue");
   inp->put<int>("number_groups",                  7);
   inp->put<std::string>("inner_solver",           "SI");
   inp->put<int>("inner_max_iters",            1000);
   inp->put<std::string>("outer_solver",              "GS");
   inp->put<int>("outer_max_iters",            1000);
   inp->put<double>("outer_tolerance",            1e-9);
   inp->put<double>("inner_tolerance",                1e-9);
   inp->put<int>("inner_print_level",          0);
   inp->put<int>("outer_print_level",          1);
   inp->put<int>("quad_number_polar_octant",       5);
   inp->put<int>("quad_number_azimuth_octant",     10);
   inp->put<std::string>("eigen_solver",       "PI");
   inp->put<std::string>("bc_west",                    "reflect");
   inp->put<std::string>("bc_east",                    "reflect");
   inp->put<std::string>("bc_north",                    "reflect");
   inp->put<std::string>("bc_south",                    "reflect");
   if (id == "pin" | id == "core")
   {
    InputDB::SP_input db(new InputDB("callow dp"));
    db->put<double>("linear_solver_atol",              1e-9);
    db->put<double>("linear_solver_rtol",              1e-9);
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
   }

   if (id == "core")
   {
    inp->put<std::string>("bc_west",                    "reflect");
    inp->put<std::string>("bc_east",                    "vacuum");
	inp->put<std::string>("bc_north",                    "reflect");
    inp->put<std::string>("bc_south",                    "vacuum");
   }
/*
   else if (id == "assembly")
   {
     InputDB::SP_input solver_db(new InputDB("solver dp"));
     solver_db->put<double>("linear_solver_atol",              0.0);
     solver_db->put<double>("linear_solver_rtol",              1e-8);
     solver_db->put<std::string>("linear_solver_type",              "petsc");
     solver_db->put<int>("linear_solver_maxit",             50000);
     solver_db->put<int>("linear_solver_gmres_restart",     30);
     solver_db->put<int>("linear_solver_monitor_level",     1);

     InputDB::SP_input preconditioner_db(new InputDB("preconditioner_db"));
     preconditioner_db->put<double>("linear_solver_atol",              0.0);
     preconditioner_db->put<double>("linear_solver_rtol",              0.1);
     preconditioner_db->put<std::string>("linear_solver_type",              "petsc");
     preconditioner_db->put<int>("linear_solver_maxit",             5000);
     preconditioner_db->put<int>("linear_solver_gmres_restart",     30);
     preconditioner_db->put<int>("linear_solver_monitor_level",     0);
     preconditioner_db->put<std::string>("pc_type",                         "petsc_pc");
     preconditioner_db->put<std::string>("petsc_pc_type",                   "ilu");
     preconditioner_db->put<int>("petsc_pc_factor_levels",          2);

     preconditioner_db->put<InputDB::SP_input>("inner_solver_db", solver_db);
     preconditioner_db->put<InputDB::SP_input>("inner_pc_db",     preconditioner_db);
     preconditioner_db->put<InputDB::SP_input>("outer_solver_db", solver_db);
     preconditioner_db->put<InputDB::SP_input>("outer_pc_db",     preconditioner_db);

   }
   */

   return inp;
 }

 Assembly::SP_assembly get_assembly(int id)
 {
   typedef RegionFactory RF;
   int n = 17;
   Assembly::SP_assembly a = Assembly::Create(n, n);
   int G = 4;
   int F = 5;
   int pin_map_a_0[] = {
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	   0,0,0,0,0,G,0,0,G,0,0,G,0,0,0,0,0,
	   0,0,0,G,0,0,0,0,0,0,0,0,0,G,0,0,0,
	   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	   0,0,G,0,0,G,0,0,G,0,0,G,0,0,G,0,0,
	   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	   0,0,G,0,0,G,0,0,F,0,0,G,0,0,G,0,0,
	   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	   0,0,G,0,0,G,0,0,G,0,0,G,0,0,G,0,0,
	   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	   0,0,0,G,0,0,0,0,0,0,0,0,0,G,0,0,0,
	   0,0,0,0,0,G,0,0,G,0,0,G,0,0,0,0,0,
	   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
                     };

   int pin_map_a_1[] = {
	   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
	   1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,
	   1,2,2,2,2,G,2,2,G,2,2,G,2,2,2,2,1,
	   1,2,2,G,2,3,3,3,3,3,3,3,2,G,2,2,1,
	   1,2,2,2,3,3,3,3,3,3,3,3,3,2,2,2,1,
	   1,2,G,3,3,G,3,3,G,3,3,G,3,3,G,2,1,
	   1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
	   1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
	   1,2,G,3,3,G,3,3,F,3,3,G,3,3,G,2,1,
	   1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
	   1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
	   1,2,G,3,3,G,3,3,G,3,3,G,3,3,G,2,1,
	   1,2,2,2,3,3,3,3,3,3,3,3,3,2,2,2,1,
	   1,2,2,G,2,3,3,3,3,3,3,3,2,G,2,2,1,
       1,2,2,2,2,G,2,2,G,2,2,G,2,2,2,2,1,
	   1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,
	   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

   int pin_map_a_2[] = {
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
	   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6
	 };

	 Assembly::vec_int pin_map(n * n);
	 for (int i = 0; i < pin_map.size(); ++i)
	 {

	 if (id == 0) pin_map[i] = pin_map_a_0[i];
	 else if (id == 1) pin_map[i] = pin_map_a_1[i];
	 else if (id == 2) pin_map[i] = pin_map_a_2[i];

	 };

	 vec_dbl radii (1, 0.54);
	 vec_int  mat_map(2, 0);
	 mat_map[0] = 0;
	 mat_map[1] = 6 ;
	 SP_pincell pin0 = get_pin(radii, mat_map);

	 mat_map[0] = 1;
	 SP_pincell pin1 = get_pin(radii, mat_map);

	 mat_map[0] = 2;
	 SP_pincell pin2 = get_pin(radii, mat_map);

	 mat_map[0] = 3;
	 SP_pincell pin3 = get_pin(radii, mat_map);

	 mat_map[0] = 4;
	 SP_pincell pin4 = get_pin(radii, mat_map);

	 mat_map[0] = 5;
	 SP_pincell pin5 = get_pin(radii, mat_map);

	 mat_map[0] = 6;
	 SP_pincell pin6 = get_pin(radii, mat_map);

	 a->add_pincell(pin0);
	 a->add_pincell(pin1);
	 a->add_pincell(pin2);
	 a->add_pincell(pin3);
	 a->add_pincell(pin4);
	 a->add_pincell(pin5);
	 a->add_pincell(pin6);

	 a->finalize(pin_map);
	 /*
     Geometry::SP_geometry geo = Geometry::Create(n*1.26, n*1.26, 0.0);
	   for (size_t i = 0; i < regions.size(); ++i)
	   {
	     geo->add_region(regions[i]);
	   }

	   return geo;

*/

  return a;
}

Core::SP_core get_core()

{
  typedef Assembly::SP_assembly         SP_assembly;
  typedef std::vector<SP_assembly>      vec_assembly;
  typedef Mesh::vec_int                 vec_int;

  Assembly::SP_assembly assem0 =  get_assembly(0);
  Assembly::SP_assembly assem1 =  get_assembly(1);
  Assembly::SP_assembly assem2 =  get_assembly(2);

  vec_assembly assemblies (3, assem0);
  assemblies[1] = assem1;
  assemblies[2] = assem2;

  int core_map_a [] = {0,1,2,
		               1,0,2,
		               2,2,2};

  vec_int core_map(3 *3);
 	 for (int i = 0; i < core_map.size(); ++i)
 	 {
       core_map[i] = core_map_a [i];
 	 };

  Core::SP_core c5g7_core = Core::Create(3);
  c5g7_core->add_assembly( assem0 );
  c5g7_core->add_assembly( assem0 );
  c5g7_core->add_assembly( assem0 );


  c5g7_core->finalize(core_map);


return c5g7_core;


}

