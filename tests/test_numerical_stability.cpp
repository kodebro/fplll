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
int test_gso_stability(ZZ_mat<mpz_t> A)
{ 
  mpfr_set_default_prec (200);


  ZZ_mat<mpz_t> U_double;
  ZZ_mat<mpz_t> UT_double;
  ZZ_mat<mpz_t> U_mpfr;
  ZZ_mat<mpz_t> UT_mpfr;

  ZZ_mat<mpz_t> B,C,A2,B2,C2,A3,B3,C3,A4,B4,C4;
  B.resize(A.r,A.c);
  C.resize(A.r,A.c); 
    B2.resize(A.r,A.c);
  C2.resize(A.r,A.c); 
    B3.resize(A.r,A.c);
  C3.resize(A.r,A.c); 
    A2.resize(A.r,A.c);
  A3.resize(A.r,A.c);  
  A4.resize(A.r,A.c); 
  B4.resize(A.r,A.c);
  C4.resize(A.r,A.c);  

  for(int i = 0; i < A.r; i++) {
    for(int j = 0; j < A.c; j++) {
        B(i,j) = A(i,j); C(i,j) = A(i,j); B2(i,j) = A(i,j); B3(i,j) = A(i,j); C2(i,j) = A(i,j); C3(i,j) = A(i,j); A2(i,j) = A(i,j); A3(i,j) = A(i,j); B4(i,j) = A(i,j); C4(i,j) = A(i,j); A4(i,j) = A(i,j);
    }
  }

  FP_NR<mpfr_t> ftmp1, ftmp2, ftmp3;
  FP_NR<mpfr_t> max_diff_gso, max_diff_givens, max_diff_mpfr;

  MatGSO<Z_NR<mpz_t>, FP_NR<double>> M_double(A, U_double, UT_double, 0);
  M_double.update_gso();

  MatGSOGivens<Z_NR<mpz_t>, FP_NR<double>> M_givens_double(B, U_double, UT_double, 0);
  M_givens_double.update_gso();
  

  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> M_mpfr(C, U_mpfr, UT_mpfr, 0);
  M_mpfr.update_gso();

  int row_begin = 1;
  int row_end = 5;

  M_mpfr.row_op_begin(row_begin,row_end+1);
  //M_mpfr.row_swap(row_begin,row_end);
  M_mpfr.row_swap(row_end,row_begin);
  M_mpfr.row_op_end(row_begin,row_end+1);

  M_double.row_op_begin(row_begin,row_end+1);
  //M_double.row_swap(row_begin,row_end);
  M_double.row_swap(row_end,row_begin);
  M_double.row_op_end(row_begin,row_end+1);

  //M_givens_double.row_swap(row_begin,row_end);
  M_givens_double.row_swap(row_end,row_begin);

  M_mpfr.update_gso();
  M_double.update_gso();


  

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

      ftmp1 = M_givens_double.r(i,j).get_data();
      ftmp2.sub(ftmp1, M_mpfr.r(i,j));
      ftmp3.abs(ftmp2);

      if ( ftmp3 > max_diff_givens )
        max_diff_givens = ftmp3;


    }
  }
  
  cerr << "Comparing the swap-accuracy of gso and givens for a matrix of dimension " << A.c << endl;
  cerr << "max_diff_gso    = " << max_diff_gso << " and max_diff_givens = " << max_diff_givens << endl;
  if (max_diff_givens  > 1) { return 1; }

  MatGSO<Z_NR<mpz_t>, FP_NR<double>> M_double2(A2, U_double, UT_double, 0);
  M_double2.update_gso();

  MatGSOGivens<Z_NR<mpz_t>, FP_NR<double>> M_givens_double2(B2, U_double, UT_double, 0);
  M_givens_double2.update_gso();
  

  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> M_mpfr2(C2, U_mpfr, UT_mpfr, 0);
  M_mpfr2.update_gso();


  M_mpfr2.row_op_begin(row_begin,row_end+1);
  //M_mpfr2.row_sub(row_end,row_begin);
  M_mpfr2.row_sub(row_begin,row_end);
  M_mpfr2.row_op_end(row_begin,row_end+1);

  M_double2.row_op_begin(row_begin,row_end+1);
  //M_double2.row_sub(row_end,row_begin);
  M_double2.row_sub(row_begin,row_end);
  M_double2.row_op_end(row_begin,row_end+1);

  M_givens_double2.row_sub(row_begin,row_end);
  //M_givens_double2.row_sub(row_end,row_begin);

  M_mpfr2.update_gso();
  M_double2.update_gso();




   /*
   cerr << "Normal GSO" << endl;
   M_double.print_mu_r_g(cerr);

   cerr << "Givens GSO" << endl;
   M_givens_double.print_mu_r_g(cerr);
   */
  
  max_diff_gso = 0.0;
  max_diff_givens = 0.0;
  max_diff_mpfr = 0.0;

  rows = A.r;

  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < i; j++)
    {
      ftmp1 = M_double2.r(i,j).get_data();
      ftmp2.sub(ftmp1, M_mpfr2.r(i,j));
      ftmp3.abs(ftmp2);

      if ( ftmp3 > max_diff_gso )
        max_diff_gso = ftmp3;

      ftmp1 = M_givens_double2.r(i,j).get_data();
      ftmp2.sub(ftmp1, M_mpfr2.r(i,j));
      ftmp3.abs(ftmp2);

      if ( ftmp3 > max_diff_givens )
        max_diff_givens = ftmp3;


    }
  }
  
  cerr << "Comparing the addmul accuracy of gso and givens for a matrix of dimension " << A.c << endl;
  cerr << "max_diff_gso    = " << max_diff_gso << " and max_diff_givens = " << max_diff_givens << endl;
  if (max_diff_givens  > 1) { return 2; }


  row_begin = 1; row_end = 5;


  MatGSO<Z_NR<mpz_t>, FP_NR<double>> M_double3(A3, U_double, UT_double, 0);
  M_double3.update_gso();

  MatGSOGivens<Z_NR<mpz_t>, FP_NR<double>> M_givens_double3(B3, U_double, UT_double, 0);
  M_givens_double3.update_gso();
  

  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> M_mpfr3(C3, U_mpfr, UT_mpfr, 0);
  M_mpfr3.update_gso();

  //cerr << M_double3.b << endl << endl;
  //cerr << M_givens_double3.b << endl << endl;


  double mult = 10.0;
  long exp = 0;

  M_mpfr3.row_op_begin(row_begin,row_end+1);
  M_mpfr3.row_addmul_we(row_begin,row_end,mult,exp);
  //M_mpfr3.row_addmul_we(row_end,row_begin,mult,exp);
  M_mpfr3.row_op_end(row_begin,row_end+1);

  M_double3.row_op_begin(row_begin,row_end+1);
  M_double3.row_addmul_we(row_begin,row_end,mult,exp);  
  //M_double3.row_addmul_we(row_end,row_begin,mult,exp);
  M_double3.row_op_end(row_begin,row_end+1);

   M_givens_double3.row_addmul_we(row_begin,row_end,mult,exp);
   //M_givens_double3.row_addmul_we(row_end,row_begin,mult,exp);

  M_mpfr3.update_gso();
  M_double3.update_gso();





  //M_mpfr3.print_mu_r_g(cerr);
  //M_double3.print_mu_r_g(cerr);
  //M_givens_double3.print_mu_r_g(cerr);

   /*
   cerr << "Normal GSO" << endl;
   M_double.print_mu_r_g(cerr);

   cerr << "Givens GSO" << endl;
   M_givens_double.print_mu_r_g(cerr);
   */
  
  max_diff_gso = 0.0;
  max_diff_givens = 0.0;
  max_diff_mpfr = 0.0;

  rows = A.r;

  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < i; j++)
    {
      ftmp1 = M_double3.r(i,j).get_data();
      ftmp2.sub(ftmp1, M_mpfr3.r(i,j));
      ftmp3.abs(ftmp2);

      if ( ftmp3 > max_diff_gso )
        max_diff_gso = ftmp3;

      ftmp1 = M_givens_double3.r(i,j).get_data();
      ftmp2.sub(ftmp1, M_mpfr3.r(i,j));
      ftmp3.abs(ftmp2);

      if ( ftmp3 > max_diff_givens )
        max_diff_givens = ftmp3;


    }
  }
  
  cerr << "Comparing the addmul_we accuracy of gso and givens for a matrix of dimension " << A.c << endl;
  cerr << "max_diff_gso    = " << max_diff_gso << " and max_diff_givens = " << max_diff_givens << endl;
  if (max_diff_givens > 1) { return 3; }



  MatGSO<Z_NR<mpz_t>, FP_NR<double>> M_double4(A4, U_double, UT_double, 0);
  M_double4.update_gso();

  MatGSOGivens<Z_NR<mpz_t>, FP_NR<double>> M_givens_double4(B4, U_double, UT_double, 0);
  M_givens_double4.update_gso();
  

  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> M_mpfr4(C4, U_mpfr, UT_mpfr, 0);
  M_mpfr4.update_gso();

  //cerr << M_double3.b << endl << endl;
  //cerr << M_givens_double3.b << endl << endl;


  M_mpfr4.row_op_begin(row_begin,row_end+1);
  M_mpfr4.move_row(row_begin,row_end);  
  //M_mpfr4.move_row(row_end,row_begin);
  M_mpfr4.row_op_end(row_begin,row_end+1);

  M_double4.row_op_begin(row_begin,row_end+1);
  M_double4.move_row(row_begin,row_end);    
  //M_double4.move_row(row_end,row_begin);
  M_double4.row_op_end(row_begin,row_end+1);


  M_givens_double4.move_row(row_begin,row_end);  
  //M_givens_double4.move_row(row_end,row_begin);

  M_mpfr4.update_gso();
  M_double4.update_gso();






  //M_mpfr3.print_mu_r_g(cerr);
  //M_double3.print_mu_r_g(cerr);
  //M_givens_double3.print_mu_r_g(cerr);

   /*
   cerr << "Normal GSO" << endl;
   M_double.print_mu_r_g(cerr);

   cerr << "Givens GSO" << endl;
   M_givens_double.print_mu_r_g(cerr);
   */
  
  max_diff_gso = 0.0;
  max_diff_givens = 0.0;
  max_diff_mpfr = 0.0;

  rows = A.r;

  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < i; j++)
    {
      ftmp1 = M_double4.r(i,j).get_data();
      ftmp2.sub(ftmp1, M_mpfr4.r(i,j));
      ftmp3.abs(ftmp2);

      if ( ftmp3 > max_diff_gso )
        max_diff_gso = ftmp3;

      ftmp1 = M_givens_double4.r(i,j).get_data();
      ftmp2.sub(ftmp1, M_mpfr4.r(i,j));
      ftmp3.abs(ftmp2);

      if ( ftmp3 > max_diff_givens )
        max_diff_givens = ftmp3;


    }
  }
  
  cerr << "Comparing the move_row accuracy of gso and givens for a matrix of dimension " << A.c << endl;
  cerr << "max_diff_gso    = " << max_diff_gso << " and max_diff_givens = " << max_diff_givens << endl;
  if (max_diff_givens > 1) { return 4; }


  return 0;

}


int main(int /*argc*/, char ** /*argv*/)
{
  ZZ_mat<mpz_t> A;
  int status = 0;

/*
  srand(1);
  int max_entry = 10;
  for (int k = 20; k < 30; k+=5) { 
    A.resize(k,k);
    for (int i = 0; i < k; i++)
      for (int j = 0; j < k; j++)
        A(i, j) = (rand()%(max_entry*2)) - max_entry;
    status += test_gso_stability(A);
  }
*/

  srand(1);
  int max_entry = 10;
  for (int k = 100; k < 300; k+=50) { 
    A.resize(k,k);
    for (int i = 0; i < k; i++)
      for (int j = 0; j < k; j++)
        A(i, j) = (rand()%(max_entry*2)) - max_entry;
    status += test_gso_stability(A);
  }




  /*for (int i = 4; i < 10; i++) { 
    A.resize(i,i);
    A.hilbert_matrix();
    test_gso_stability(A);
  } */ 


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
