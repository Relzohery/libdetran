/*
 * DEIM.cc
 *
 *  Created on: Jan 1, 2021
 *      Author: rabab
 */

#include "DEIM.hh"

namespace detran
{

DEIM::DEIM(SP_matrixDense U, int r)
:d_r(r),
 d_U(U),
 d_n(U->number_rows())
{

}


void DEIM::Search()
{
  db = get_db();
  d_solver = LinearSolverCreator::Create(db);
  d_l = new int[d_r];
  int tmp = 0;
  int l = 0;
  // get the first interpolation point
  for (int j=0; j<d_n-1; j++)
  {
	double v1 = (*d_U)(j+1, 0);
	double v2 = (*d_U)(tmp, 0);

    if (std::abs(v1) > std::abs(v2))
    {
      tmp = j+1;
      l = j+1;
    }
  }
  d_l[0] = l;

  double c = (*d_U)(l, 1)/(*d_U)(l, 0);
  l = 0;
  double tmp_r = 0.0;
  double r0 = (*d_U)(0, 1) - c*(*d_U)(0, 0);
  for (int i=1; i<d_n-1; i++)
  {
   double r = (*d_U)(i, 1) - c*(*d_U)(i, 0);
   if (std::abs(r) > std::abs(r0))
   {
     l = i;
     r0 = r;
   }
  }
  d_l[1] = l;

 // now, need linear solver
  for (int i=2; i<d_r+1; i++)
  {
	d_solver = LinearSolverCreator::Create(db);
	int M = i;
    SP_matrixDense Ur;
    SP_matrixDense Urm;
    ///
    Ur = new MatrixDense(M, M);
    Urm = new MatrixDense(d_n, M);

    SP_vector b;
    b = new callow::Vector(M, 0.0);
    SP_vector x;
    x = new callow::Vector(M, 0.0);

    for (int k=0; k<M; k++)
    {
	 (*b)[k] = (*d_U)(d_l[k], i);
     for (int m=0; m<M; m++)
     {
       (*Ur)(k, m) = (*d_U)(d_l[k], m);
     }
    }

   for (int n=0; n<d_n; n++)
   {
    for (int m=0; m<M; m++)
    {
      (*Urm)(n, m) = (*d_U)(n, m);
    }
   }
   d_Um = Ur;
   d_solver->set_operators(Ur, db);
   d_solver->solve(*b, *x);

   Vector y(d_n, 0.0);
   Urm->multiply(*x, y);
   l = 0;
   double r0 = (*d_U)(0, i) - y[0];
   for (int ii=1; ii<d_n; ii++)
   {
     double r = (*d_U)(ii, i) - y[ii];
     if (std::abs(r) > std::abs(r0))
     {
       l = ii;
       r0 = r;
     }
   }
   d_l[i] = l;
 }
}
}
