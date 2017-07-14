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
template <class ZT, class FT> Matrix<FT> matrix_relative_difference(Matrix<FT> r1, Matrix<FT> r2)
{
  Matrix<FT> diff_matrix = Matrix<FT>(r1.get_rows(), r1.get_cols());
  diff_matrix.fill(0.0);
  FT relativation_factor = 0.0;
  for (int i = 0; i < r1.get_rows(); i++)
  {
    for (int j = 0; j < i; j++)
    {  // j < i, because r is lower-triangular, and has only 1 on the diagonal.
      relativation_factor = abs(r1[i][j]) + abs(r2[i][j]);
      if (relativation_factor.is_zero())
      {
        diff_matrix[i][j] = abs(r1[i][j] - r2[i][j]);
      }
      else
      {
        diff_matrix[i][j] = abs(r1[i][j] - r2[i][j]) / relativation_factor;
      }
    }
  }
  return diff_matrix;
}

// Returns true when the r-matrices of M1 and M2 are entry-wise equal, up to an error 'error'.
template <class ZT, class FT> bool rs_are_equal(MatGSO<ZT, FT> M1, MatGSOGram<ZT, FT> M2, FT error)
{
  Matrix<FT> r1   = M1.get_r_matrix();
  Matrix<FT> r2   = M2.get_r_matrix();
  Matrix<FT> diff = matrix_relative_difference<ZT, FT>(r1, r2);

  FT max_entry = 0.0;
  max_entry    = diff.get_max();
  if (max_entry > error)
  {
    diff.print(cerr);
    cerr << endl << endl;
    return false;
  }
  return true;
}
*/

template <class ZT, class FT> Matrix<FT> matrix_difference(Matrix<FT> mu1, Matrix<FT> mu2)
{
  Matrix<FT> diff_matrix = Matrix<FT>(mu1.get_rows(), mu1.get_cols());
  diff_matrix.fill(0.0);
  FT relativation_factor = 0.0;
  for (int i = 0; i < mu1.get_rows(); i++)
  {
    for (int j = 0; j < i; j++)
    { 
       diff_matrix(i,j)= abs(mu1(i,j) - mu2(i,j));
    }
  }
  return diff_matrix;
}

template <class ZT, class FT> bool mus_are_equal(MatGSO<ZT, FT> M1, MatGSOGram<ZT, FT> M2, FT error)
{
  Matrix<FT> mu1  = M1.get_mu_matrix();
  Matrix<FT> mu2 = M2.get_mu_matrix();
  Matrix<FT> diff = matrix_difference<ZT, FT>(mu1, mu2);

  FT max_entry = 0.0;
  max_entry    = diff.get_max();
  if (max_entry > error)
  {
    cerr << "Difference is too big:" << endl;
    //diff.print(cerr);
    /*cerr << endl << endl;
    mu1.print(cerr);
    cerr << endl << endl;
    mu2.print(cerr);
    cerr << endl << endl;*/

    return false;
  }
  cerr << "Maximum difference: " << max_entry << endl;
  return true;
}



template <class ZT> void read_matrix(ZZ_mat<ZT> &A, const char *input_filename)
{
  ifstream is(input_filename);
  if (!is)
  {
    cerr << "Could not open file!" << endl;
  }  // throw std::runtime_error("could not open input file");
  is >> A;
}

template <class ZT, class FT> int test_gso(ZZ_mat<ZT> &A)
{
  // TEST A
  // Method:
  // (1) Apply 'normal' MatGSO to A, with floating point gram matrix
  // -   Extract mu-matrix of A
  // (2) Apply 'normal' MatGSO to A, with exact gram matrix
  //     Extract mu-matrix of A

  // - Compute G = A^T A.
  // - Apply GramMatGSO to G.
  // - Extract mu-matrix of this 


  // -> The r-matrices should be equal.

  // TEST B
  // Apply some 'random' elementary operation on A and on G
  // (of course, the operations on A and G are only 'abstractly' the same)
  // check if their r-matrices are still equal.

  // TEST A
  // ----------------------------
  int r = A.r;

  ZZ_mat<ZT> U;
  ZZ_mat<ZT> UT;

  MatGSO<Z_NR<ZT>, FP_NR<FT>> Mbuf(A, U, UT, GSO_INT_GRAM);
  Mbuf.discover_all_rows();
  Matrix<Z_NR<ZT>> G = Mbuf.get_g_matrix();

  MatGSO<Z_NR<ZT>, FP_NR<FT>> M(A, U, UT, GSO_DEFAULT); // with floating point gram matrix
  M.update_gso();
  ZZ_mat<ZT> A1(A);
  MatGSO<Z_NR<ZT>, FP_NR<FT>> M2(A1, U, UT, GSO_INT_GRAM); //  with exact gram matrix
  M2.update_gso();
  MatGSOGram<Z_NR<ZT>, FP_NR<FT>> M3(G, U, UT, GSO_INT_GRAM); // with only a gram matrix, and no basis (GramGSO object)
  M3.update_gso();

  FP_NR<FT> err  = .001;
  bool retvalue1 = mus_are_equal(M, M3, err) && mus_are_equal(M2, M3, err);

  // TEST B
  // ------------------------

  for (int i = 0; i < rand() % 10 + 1; i++)
  {
    int k = rand() % r;
    int j = rand() % r;
    M.move_row(k, j);
    M2.move_row(k, j);
    M3.move_row(k, j);
  }
  M.update_gso();
  M2.update_gso();
  M3.update_gso();
  bool retvalue2 = mus_are_equal(M, M3, err) && mus_are_equal(M2, M3, err) ;

  M.row_op_begin(0, r);
  M2.row_op_begin(0, r);
  M3.row_op_begin(0, r);
  for (int i = 0; i < rand() % 10 + 1; i++)
  {
    int k = rand() % r;
    int j = rand() % r;
    M.row_add(k, j);
    M2.row_add(k, j);
    M3.row_add(k, j);
  }
  M.row_op_end(0, r);
  M2.row_op_end(0, r);
  M3.row_op_end(0, r);

  M.update_gso();
  M2.update_gso();
  M3.update_gso();
  bool retvalue3 = mus_are_equal(M, M3, err) && mus_are_equal(M2, M3, err);

  return (!retvalue1) * 1 + (!retvalue2) * 2 + (!retvalue3) * 4;
}

template <class ZT, class FT> int test_filename(const char *input_filename)
{
  ZZ_mat<ZT> A;
  read_matrix(A, input_filename);
  int retvalue = test_gso<ZT, FT>(A);
  if (retvalue & 1)
  {
    cerr
        << input_filename
        << " shows different GSO-outputs for grammatrix representation, basis representation (fp gram) or basis representation (exact gram).\n";
  }
  if (retvalue & 2)
  {
    cerr << input_filename << " shows different GSO-outputs for grammatrix representation and "
                              " basis representation (fp gram) or basis representation (exact gram) after moving rows.\n";
  }
  if (retvalue & 4)
  {
    cerr << input_filename << " shows different GSO-outputs for grammatrix representation and "
                              " basis representation (fp gram) or basis representation (exact gram) after adding rows.\n";
  }
  if (retvalue > 0)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

/**
   @brief Construct d × (d+1) integer relations matrix with bit size b and test LLL.

   @param d                dimension
   @param b                bit size
   @param method           LLL method to test
   @param float_type       floating point type to test
   @param flags            flags to use
   @param prec             precision used for is_lll_reduced

   @return zero on success
*/

template <class ZT, class FT>
int test_int_rel(int d, int b, FloatType float_type = FT_DEFAULT, int prec = 0)
{
  ZZ_mat<ZT> A;
  A.resize(d, d + 1);
  A.gen_intrel(b);
  int retvalue = test_gso<ZT, FT>(A);
  if (retvalue & 1)
  {
    cerr
        << "Integer relation matrix with parameters " << d << " and " << b
        << " shows different GSO-outputs for grammatrix representation, basis representation (fp gram) or basis representation (exact gram).\n";
  }
  if (retvalue & 2)
  {
    cerr
        << "Integer relation matrix with parameters " << d << " and " << b
        << " shows different GSO-outputs for grammatrix representation and "
                              " basis representation (fp gram) or basis representation (exact gram) after moving rows.\n";
  }
  if (retvalue & 4)
  {
    cerr
        << "Integer relation matrix with parameters " << d << " and " << b
         << " shows different GSO-outputs for grammatrix representation and "
                              " basis representation (fp gram) or basis representation (exact gram) after adding rows.\n";
  }
  if (retvalue > 0)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int main(int /*argc*/, char ** /*argv*/)
{

  int status = 0;

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
