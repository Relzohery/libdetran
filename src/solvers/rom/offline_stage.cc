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

 //DEIM D(d_U, d_r);
//   D.Search();
//   Interpolation_indices = D.interpolation_indices();
//   Ur_deim = new callow::MatrixDense(d_r, d_r);
//   Ur_deim = D.ReducedBasis();
}

vec_matrix offline_stage::Decompositon()
{
  detran::OperatorProjection projector(1);

  int rank = d_Uf->number_columns();

  Vector y(d_n, 0.0);
  Vector y2(rank, 0.0);
  // convert each vector to matrix and project
  for (int r=0; r< d_r; r++)
  {
    SP_matrixDense Ar;
	Ar = new callow::MatrixDense(rank, rank);
	// loop over each column of Uf
    for (int rb=0; rb<rank ; rb++)
    {
      int nz = 0;
      for (int i = 0; i < d_m; ++i)
      {
        double temp = 0.0;
        // for all columns
        for (int p = d_rows[i]; p < d_rows[i + 1]; ++p)
        {
          int j = d_operator->columns()[p];

          temp += (*d_Uf)(j, rb )*(*d_U)(nz, r);
          nz += 1;
        }
        y[i] = temp;
      }

     d_Uf->multiply_transpose(y, y2);
     double *y_ = &y2[0];
     Ar->insert_col(rb, y_, 0);
   }
    Aq.push_back(Ar);

  } // columns of U-DEIM

return Aq;
}



void offline_stage::map_indices()
{
  DEIM D(d_U, d_r);
  D.Search();
  l = D.interpolation_indices();

  d_target_rows = new int[d_r];
  d_target_cols = new int[d_r];

  int nz = 0;
  for (int i=0; i < d_r; i++)
  {
    d_target_cols[i] = d_cols[l[i]];
    for (int j=0; j < d_m; j++)
    {
      //std::cout << l[i] << "  " <<  d_rows[j] <<  "  " << d_rows[j+1] << " &&&&& \n";
      if (d_rows[j] <= l[i] && l[i] <= d_rows[j+1])
      {
        d_target_rows[i] = j;
      }
    }
  }
}
