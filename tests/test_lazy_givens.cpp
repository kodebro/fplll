/* Copyright (C) 2015 Martin Albrecht

   This file is part of fplll. fplll is free software: you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation,
   either version 2.1 of the License, or (at your option) any later version.

   fplll is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with fplll. If not, see <http://www.gnu.org/licenses/>. */

#include <cstring>
#include <fplll.h>
#include <gso.h>
#include <gso_gram.h>
#include <gso_interface.h>
#include <gso_givens.h>
#include <lll.h>
#include <test_utils.h>

using namespace std;
using namespace fplll;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif




template <class ZT, class FT> void lazy_givens(MatGSOGivens<Z_NR<ZT>,FP_NR<FT>> m, int dim)
{
  int cols = m.b.get_cols();
  int rows = m.b.get_rows();
  //Matrix<ZN_R<ZT>> b,U,UT; 
  ZZ_mat<ZT> b,U,UT;

        MatGSOGivens<Z_NR<ZT>, FP_NR<FT>> M(b, U, UT, 1 << 8 | GSO_ROW_EXPO);
        Mlazy.update_gso();
        LLLReduction<Z_NR<ZT>, FP_NR<FT>> LLLObjgram(Mlazy, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_VERBOSE);
        LLLObjgram.lll();


  //int startrow = 0;
  for(int startrow = 0; startrow < rows - dim; startrow+=1) {
        cerr << startrow << endl;
        b.resize(min(dim, rows - startrow),cols);
        //U.resize(min(dim, rows - startrow),cols);
        //U.gen_identity(min(dim, rows - startrow));

        for(int i = startrow; i < min(dim + startrow, rows); i++)
          for(int j = 0; j < cols; j++)
            b(i-startrow,j) = m.b(i,j);

       // cerr << b << endl << endl;
      

        MatGSOGivens<Z_NR<ZT>, FP_NR<FT>> Mlazy(b, U, UT, 1 << 8 | GSO_ROW_EXPO);
        Mlazy.update_gso();
        LLLReduction<Z_NR<ZT>, FP_NR<FT>> LLLObjgram(Mlazy, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_VERBOSE);
        LLLObjgram.lll();
        //cerr << U << endl;

        cerr << b[0] << endl << endl;

        for(int i = startrow; i < min(dim + startrow, rows); i++)
          for(int j = 0; j < cols; j++)
            m.b(i,j) = b(i-startrow,j);

    }

}




template <class ZT, class FT> int test_lll(ZZ_mat<ZT> &A)
{

  ZZ_mat<ZT> U;
  ZZ_mat<ZT> UT;

  // _______________________________________________
  // -----------------------------------------------
  // Create copy of A.

  ZZ_mat<ZT> A2;
  int r = A.r;
  int c = A.c;
  A2.resize(r, c);
  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < c; j++)
    {
      A2(i,j) = A(i,j);
    }
  }

  MatGSO<Z_NR<ZT>, FP_NR<FT>> M(A, U, UT, 0);
  MatGSOGivens<Z_NR<ZT>, FP_NR<FT>> MGivens(A2, U, UT, (1 << 8) || GSO_ROW_EXPO);
  // LAZY GIVENS

  lazy_givens<ZT, FT>(MGivens,40);
  //lazy_givens<Z_NR<ZT>, FP_NR<FT>>(MGivens,100);
  //lazy_givens<Z_NR<ZT>, FP_NR<FT>>(MGivens,200);  
  //lazy_givens<Z_NR<ZT>, FP_NR<FT>>(MGivens,60);
  //lazy_givens<Z_NR<ZT>, FP_NR<FT>>(MGivens,80);    
  return 1;

}

template <class ZT, class FT> int test_filename(const char *input_filename)
{
  ZZ_mat<ZT> A;
  int status = 0;
  status |= read_matrix(A, input_filename);
  status |= test_lll<ZT, FT>(A);
  return status;
}

/**
   @brief Construct d × (d+1) integer relations matrix with bit size b and test LLL.

   @param d                dimension
   @param b                bit size

   @return zero on success
*/

template <class ZT, class FT> int test_int_rel(int d, int b)
{
  ZZ_mat<ZT> A;
  A.resize(d, d + 1);
  A.gen_intrel(b);
  return test_lll<ZT, FT>(A);
}

int main(int /*argc*/, char ** /*argv*/)
{

  int status = 0;


  status |= test_int_rel<mpz_t, double>(200,1000);


  if (status == 0)
  {
    cerr << "All tests passed." << endl;
    return 0;
  }
  else
  {
    return -1;
  }

  return 0;
}
