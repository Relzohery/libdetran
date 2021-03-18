/*
 * offline_stage.hh
 *
 *  Created on: Jan 1, 2021
 *      Author: rabab
 */

//#ifndef SOURCE_SRC_SOLVERS_ROM_OFFLINE_STAGE_HH_
//#define SOURCE_SRC_SOLVERS_ROM_OFFLINE_STAGE_HH_

#include "utilities/InputDB.hh"
#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "solvers/rom/DEIM.hh"

class offline_stage
{
  public:
	typedef detran_utilities::InputDB::SP_input       SP_input;
	typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
	typedef detran_material::Material::SP_material    SP_material;
	typedef Matrix::SP_matrix                         SP_matrix;
	typedef MatrixDense::SP_matrix                    SP_matrixDense;
	typedef std::vector<SP_matrixDense>               vec_matrix;

	offline_stage(SP_matrix A, SP_matrixDense Uf, SP_matrixDense UD, unsigned int r);

	vec_matrix Decompositon();

	double* Interpolation_indices;

	SP_matrixDense Ur_deim;

  private:
	int d_nnz;
	int* d_cols;
	int* d_rows;
	int d_m;
	int d_n;
	int d_r;
    vec_matrix Aq;
    //SP_matrixDense UD;
    SP_matrix d_operator;
    SP_matrixDense d_U;
    SP_matrixDense d_Uf;

};



//#endif /* SOURCE_SRC_SOLVERS_ROM_OFFLINE_STAGE_HH_ */