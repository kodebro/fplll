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
#include <stdlib.h>     /* srand, rand */
#include <gso.h>
#include <gso_gram.h>
#include <gso_interface.h>
#include <gso_givens.h>
#include <nr/matrix.h>
#include <test_utils.h>
//#include <random>

using namespace std;
using namespace fplll;

#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif


/*
void test_gso_stability(int rows, int cols, int max_entry)
{ 
  mpfr_set_default_prec (200);

  ZZ_mat<mpz_t> A;
  A.resize(rows, cols);

  ZZ_mat<mpz_t> U_double;
  ZZ_mat<mpz_t> UT_double;
  ZZ_mat<mpz_t> U_mpfr;
  ZZ_mat<mpz_t> UT_mpfr;

  FP_NR<mpfr_t> ftmp1, ftmp2, ftmp3;
  FP_NR<mpfr_t> max_diff_gso, max_diff_givens, max_diff_mpfr;

  srand (1);
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < rows; j++)
      A(i, j) = (rand()%(max_entry*2)) - max_entry;

  MatGSO<Z_NR<mpz_t>, FP_NR<double>> M_double(A, U_double, UT_double, 0);
  M_double.update_gso();


  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> M_mpfr(A, U_mpfr, UT_mpfr, 0);
  M_mpfr.update_gso();

  max_diff_gso = 0.0;
  max_diff_givens = 0.0;
  max_diff_mpfr = 0.0;
  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < i; j++)
    {
      ftmp1 = M_double.r(i,j).get_data();
      ftmp2.sub(ftmp1, M_mpfr.r(i,j));
      ftmp3.abs(ftmp2);
      
      if ( ftmp3 > max_diff_gso )
        max_diff_gso = ftmp3;

    }
  }

  cerr << "Comparing the accuracy of gso and givens for a matrix of " << rows << " rows and " << cols << " columns." << endl;
  cerr << "max_diff_gso    = " << max_diff_gso << endl;

}
*/
template<class ZT, class FT> int test_gso_stability(ZZ_mat<ZT> A)
{ 
  mpfr_set_default_prec (200);


  //ZZ_mat<mpz_t> U_givens;
  //ZZ_mat<mpz_t> UT_givens;
  //ZZ_mat<mpz_t> U_double;
  //ZZ_mat<mpz_t> UT_double;  
  //ZZ_mat<mpz_t> U_mpfr;
  //ZZ_mat<mpz_t> UT_mpfr;

  //ZZ_mat<mpz_t> B;
  //B.resize(A.r,A.c);
  Matrix<Z_NR<ZT>> B, UB, UTB, U, UT;
  B.resize(A.r,A.c);

  Matrix<FP_NR<FT>> displayMatrix;
  displayMatrix.resize(A.r,A.c);


  for(int i = 0; i < A.r; i++) {
    for(int j = 0; j < A.c; j++) {
        B(i,j) = A(i,j); 
    }
  }

  FP_NR<FT> ftmp1, ftmp2, ftmp3;
  FP_NR<FT> max_diff_gso, max_diff_givens, max_diff_mpfr;

  //MatGSO<Z_NR<mpz_t>, FP_NR<double>> M_double(A, U_double, UT_double, 0);
  //M_double.update_gso();

  MatGSOGivens<Z_NR<ZT>, FP_NR<double>> M_givens(A, UB, UTB, 0);
  MatGSO<Z_NR<ZT>, FP_NR<FT>> M_mpfr(B, U, UT, 0);
  M_givens.recomputation_count = 1000; //00000;


  M_mpfr.update_gso();
  M_givens.update_gso();

  max_diff_gso = 0.0;
  max_diff_givens = 0.0;
  max_diff_mpfr = 0.0;
  int rows = A.r;

  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < i; j++)
    {
      //ftmp1 = M_double.r(i,j).get_data();
      //ftmp2.sub(ftmp1, M_mpfr.r(i,j));
      //ftmp3.abs(ftmp2);

      //if ( ftmp3 > max_diff_gso )
      //  max_diff_gso = ftmp3;

      ftmp1 = M_givens.mu(i,j).get_data();
      //(M_givens.r(i,j)).get_mpfr(ftmp1);
      ftmp2.sub(ftmp1, M_mpfr.mu(i,j));
      ftmp3.abs(ftmp2);
      displayMatrix(i,j) = ftmp3;

      if ( ftmp3 > max_diff_givens )
        max_diff_givens = ftmp3;


    }
  }


  
  cerr << "Comparing the GSO- accuracy of gso and givens for a matrix of dimension " << A.c << endl;
  cerr << "max_diff_gso    = " << max_diff_gso << " and max_diff_givens = " << max_diff_givens << endl;
  //cerr << displayMatrix << endl;

  //if (max_diff_givens > 1) { return 1; }





  LLLReduction<Z_NR<ZT>, FP_NR<FT>> LLLObj(M_mpfr, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_VERBOSE);
  LLLReduction<Z_NR<ZT>, FP_NR<double>> LLLObj_givens(M_givens, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_VERBOSE);


  // and LLL reduce both objects
  LLLObj_givens.lll();

  LLLObj.lll();


  // ------------------------------------------------
  // ************************************************

  // _________________________________________________
  // -------------------------------------------------
  // Check whether M and MGivens are really reduced after LLL reduction
  int is_reduced  = is_lll_reduced<Z_NR<ZT>, FP_NR<FT>>(M_mpfr, LLL_DEF_DELTA, LLL_DEF_ETA);

  MatGSO<Z_NR<ZT>, FP_NR<FT>> M2(M_givens.b, U, UT, 0);
  M2.update_gso();
  int is_greduced = is_lll_reduced<Z_NR<ZT>, FP_NR<FT>>(M2, LLL_DEF_DELTA, LLL_DEF_ETA);

  if (is_reduced != 1 || is_greduced != 1)
  {
    if (is_reduced != 1)
    {
      cerr << "The basis GSO-object is not LLL-reduced after calling LLL\n";
    }
    if (is_greduced != 1)
    {
      cerr << "The givens GSO-object is not LLL-reduced after calling LLL\n";
    }

  }

//cerr << M_mpfr.b << endl << endl;
//cerr << M_givens.b << endl;
  

  max_diff_gso = 0.0;
  max_diff_givens = 0.0;
  max_diff_mpfr = 0.0;
  rows = A.r;

  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < i; j++)
    {
      //ftmp1 = M_double.r(i,j).get_data();
      //ftmp2.sub(ftmp1, M_mpfr.r(i,j));
      //ftmp3.abs(ftmp2);

      //if ( ftmp3 > max_diff_gso )
      //  max_diff_gso = ftmp3;

      ftmp1 = M_givens.mu(i,j).get_data();
      //(M_givens.r(i,j)).get_mpfr(ftmp1);
      ftmp2.sub(ftmp1, M2.mu(i,j));
      ftmp3.abs(ftmp2);
      displayMatrix(i,j) = ftmp3;

      if ( ftmp3 > max_diff_givens )
        max_diff_givens = ftmp3;


    }
  }
  
  cerr << "Comparing the LLL-accuracy of gso and givens for a matrix of dimension " << A.c << endl;
  cerr << "max_diff_gso    = " << max_diff_gso << " and max_diff_givens = " << max_diff_givens << endl;
    if (is_greduced !=1) {     
    for(int i = 0; i < rows; i++)
    {
      for(int j = 0; j < i; j++)
      {
        if (displayMatrix(i,j) < .1) { cerr << "0"; } else  {cerr << "x";}
      }
      cerr << endl;
    }
    cerr << endl;
    cerr << displayMatrix << endl << endl;
    return 1; 

    }

  if (max_diff_givens > 1) { 
    for(int i = 0; i < rows; i++)
    {
      for(int j = 0; j < i; j++)
      {
        if (displayMatrix(i,j) < .001) { cerr << "0"; } else  {cerr << "x";}
      }
      cerr << endl;
    }
    cerr << endl;
    return 0; }
  //if (max_diff_givens  > 1) { return 1; }

  return 0;
}


template <class ZT, class FT> int test_filename(const char *input_filename)
{
  ZZ_mat<ZT> A;
  int status = read_matrix(A, input_filename);
  // if status == 1, read_matrix fails.
  if (status == 1)
  {
    return 1;
  }
  int retvalue = test_gso_stability<ZT, FT>(A);
  return retvalue;
}

template <class ZT, class FT> int test_int_rel(int d, int b)
{
  ZZ_mat<ZT> A;
  A.resize(d, d + 1);
  A.gen_intrel(b);
  return test_gso_stability<ZT, FT>(A);
}

int main(int /*argc*/, char ** /*argv*/)
{
  ZZ_mat<mpz_t> A;
  int status = 0;

/*
  srand(1);
  int max_entry = 1000000000;
  for (int k = 200; k < 300; k+=20) { 
    A.resize(k,k);
    for (int i = 0; i < k; i++)
      for (int j = 0; j < k; j++)
        A(i, j) = (rand()%(max_entry*2)) - max_entry;
    status += test_gso_stability<mpz_t, mpfr_t>(A);
  }
  */
 // status |= test_int_rel<mpz_t, double>(90, 10);

/*
  srand(1);
  int max_entry = 10;
  for (int k = 100; k < 300; k+=50) { 
    A.resize(k,k);
    for (int i = 0; i < k; i++)
      for (int j = 0; j < k; j++)
        A(i, j) = (rand()%(max_entry*2)) - max_entry;
    status += test_gso_stability<mpz_t, mpfr_t>(A);
  }

*/


  /*for (int i = 4; i < 10; i++) { 
    A.resize(i,i);
    A.hilbert_matrix();
    test_gso_stability(A);
  } */ 

  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t, double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, double>(50, 20);
  status |= test_int_rel<mpz_t, double>(40, 10);

  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t, mpfr_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, mpfr_t>(50, 20);
  status |= test_int_rel<mpz_t, mpfr_t>(40, 10);
  


#ifdef FPLLL_WITH_LONG_DOUBLE
  status |= test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |=
      test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |=
      test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |=
      test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |=
      test_filename<mpz_t, long double>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, long double>(50, 20);
  status |= test_int_rel<mpz_t, long double>(40, 10);
#endif
#ifdef FPLLL_WITH_QD
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t, dd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, dd_real>(50, 20);
  status |= test_int_rel<mpz_t, dd_real>(40, 10);

  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t, qd_real>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, qd_real>(50, 20);
  status |= test_int_rel<mpz_t, qd_real>(40, 10);
#endif
#ifdef FPLLL_WITH_DPE
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example2_in");
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice");
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice2");
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice3");
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice4");
  status |= test_filename<mpz_t, dpe_t>(TESTDATADIR "/tests/lattices/example_cvp_in_lattice5");
  status |= test_int_rel<mpz_t, dpe_t>(50, 20);
  status |= test_int_rel<mpz_t, dpe_t>(40, 10);
#endif

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
