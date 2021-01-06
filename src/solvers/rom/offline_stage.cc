/*
 * offline_stage.cc
 *
 *  Created on: Jan 1, 2021
 *      Author: rabab
 */



#include "offline_stage.hh"
typedef MatrixDense::SP_matrix                    SP_matrixDense;
typedef std::vector<SP_matrixDense>               vec_matrix;

offline_stage::offline_stage(SP_matrix A, SP_matrixDense Uf, SP_matrixDense UD, unsigned int r)
:
 d_cols(A->columns()),
 d_rows(A->rows()),
 d_nnz(A->number_nonzeros()),
 d_r(r),
 d_operator(A),
 d_n(A->number_columns()),
 d_m(A->number_rows()),
 d_U(UD),
 d_Uf(Uf)

{
   //Ensure(A->number_nonzeros() == UD->number_rows());
   int d_n = A->number_columns();
   int d_m = A->number_rows();
}

vec_matrix offline_stage::VectorToMatrix()
{

  detran::OperatorProjection projector(1);
  for (int r=0; r< d_r; r++)
  { int nz = 0;
	callow::Matrix::SP_matrix d_A_prime;
	d_A_prime = new callow::Matrix(d_m, d_n, d_operator->number_nonzeros());
    Matrix Ai(*d_operator);
    for (int i=0; i<d_m; i++)
    {
      for (int c= d_operator->rows()[i]; c< d_operator->rows()[i+1]; c++)
      {
        double v = (*d_U)(nz, r);
        d_A_prime->insert(i, d_operator->columns()[c], v);
        nz += 1;
      }
    }
    d_A_prime->assemble();
    d_A_prime->print_matlab("A1.txt");

    projector.SetOperators(d_A_prime, d_Uf);

    d_Uf->print_matlab("flux_basis.txt");

    int rank = d_Uf->number_columns();
    SP_matrixDense Ar;
    Ar = new callow::MatrixDense(rank, rank);

    projector.Project(Ar);
    Ar->print_matlab("expanded_matrices.txt");
    Aq.push_back(Ar);
  }
return Aq;
}

void offline_stage::project()
{


}
