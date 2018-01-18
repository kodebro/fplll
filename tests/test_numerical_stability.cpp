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
#include <nr/matrix.h>
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
void test_gso_stability(ZZ_mat<mpz_t> A)
{ 
  mpfr_set_default_prec (200);


  ZZ_mat<mpz_t> U_double;
  ZZ_mat<mpz_t> UT_double;
  ZZ_mat<mpz_t> U_mpfr;
  ZZ_mat<mpz_t> UT_mpfr;

  FP_NR<mpfr_t> ftmp1, ftmp2, ftmp3;
  FP_NR<mpfr_t> max_diff_gso, max_diff_givens, max_diff_mpfr;

  MatGSO<Z_NR<mpz_t>, FP_NR<double>> M_double(A, U_double, UT_double, 0);
  M_double.update_gso();


  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> M_mpfr(A, U_mpfr, UT_mpfr, 0);
  M_mpfr.update_gso();

  max_diff_gso = 0.0;
  max_diff_givens = 0.0;
  max_diff_mpfr = 0.0;

  int rows = A.r;

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

  cerr << "Comparing the accuracy of gso and givens for a matrix of dimension " << A.c << endl;
  cerr << "max_diff_gso    = " << max_diff_gso << endl;

}


int main(int /*argc*/, char ** /*argv*/)
{
  ZZ_mat<mpz_t> A;
  for (int i = 4; i < 10; i++) { 
    A.resize(i,i);
    A.hilbert_matrix();
    test_gso_stability(A);
  }

  //for (int i = 4; i < 10; i++)
  //  test_gso_stability(1<<i, 1<<i, 1000);
  //
  int status = 0;


/*
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
*/
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
