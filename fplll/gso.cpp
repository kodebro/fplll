/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2011 Xavier Pujol.

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

/* Template source file */

#include "gso.h"

FPLLL_BEGIN_NAMESPACE

template <class ZT, class FT> void MatGSO<ZT, FT>::initialize_r_givens_matrix()
{
  mu_givens.fill(0.0);

  if (enable_row_expo)
  {
    // throw std::runtime_error("Error: givens rotations are not yet implemented for enable_row_expo
    // is true.");
  }
  else
  {
    if (r_givens.get_rows() != d)
    {
      throw std::runtime_error("Error: r_givens does not have good dimensions.");
    }
    if (r_givens.get_cols() != b.get_cols())
    {
      throw std::runtime_error("Error: r_givens does not have good dimensions.");
    }
    for (int i = 0; i < d; i++)
    {
      for (int j = 0; j < r_givens.get_cols(); j++)
      {
        r_givens(i, j).set_z(b(i, j));
      }
    }
  }
}

template <class ZT, class FT> void MatGSO<ZT, FT>::givens_rotation(int col_i, int col_j, int row_k)
{
  // TODO: This maybe can be sped up if we manage to remove some FT 's
  // we might want to have some more 'tmp1';


  // If the place where we want to introduce
  // a zero is already zero, we do nothing.
  if (r_givens(row_k,col_j).is_zero()) {
    return;
  }

  // Computes the 'c' and 's' of the givens, 
  // by means of the hypot-function,
  // which hopefully avoids overflow.
  FT c, s;
  ftmp1.hypot(r_givens(row_k, col_i), r_givens(row_k, col_j));
  c.div(r_givens(row_k, col_i), ftmp1);
  s.div(r_givens(row_k, col_j), ftmp1);


  // Here we take multiples of columns
  // -  Remark that we leave alone the first
  //    row_k - 1 entries of the columns.
  //    For effectivity this is obviously better,
  //    but for numerical stability I don't know.
  // -  Also, we could force the entry that should be
  //    zero as 'zero'. Now know whether this is good or not.
  //
  // -> Answer to both -> Probably no effect on stability.

  for (int k = row_k; k < r_givens.get_rows(); k++)
  {
    ftmp1 = r_givens(k, col_i); 
    ftmp2 = r_givens(k, col_j); 
    r_givens(k, col_i).mul(ftmp1, c); // r_(k,col_i) = c*r_(k,col_i) + s*r_(k,col_j)
    r_givens(k, col_i).addmul(ftmp2, s);

    r_givens(k, col_j).neg(s);
    r_givens(k, col_j).mul(ftmp1, r_givens(k, col_j));
    r_givens(k, col_j).addmul(ftmp2, c); // r_(k,col_j) = -s*r_(k,col_i) + c*r_(k,col_j)
  }
  // "Forcing" zero 
  r_givens(row_k, col_j) = 0.0;
}

template <class ZT, class FT> void MatGSO<ZT, FT>::givens_row_reduction(int row_k, int rightmost_nonzero_entry)
{
  for (int i = rightmost_nonzero_entry - 1; i > row_k; i--)
    givens_rotation(i - 1, i, row_k);
  for (int i = row_k; i < mu_givens.get_rows(); i++)
    mu_givens(i,row_k).div(r_givens(i,row_k),r_givens(row_k,row_k));
}

template <class ZT, class FT> void MatGSO<ZT, FT>::update_bf(int i)
{
  int n = max(n_known_cols, init_row_size[i]);
  if (enable_row_expo)
  {
    long max_expo = LONG_MIN;
    for (int j = 0; j < n; j++)
    {
      b(i, j).get_f_exp(bf(i, j), tmp_col_expo[j]);
      max_expo = max(max_expo, tmp_col_expo[j]);
    }
    for (int j = 0; j < n; j++)
    {
      bf(i, j).mul_2si(bf(i, j), tmp_col_expo[j] - max_expo);
    }
    row_expo[i] = max_expo;
  }
  else
  {
    for (int j = 0; j < n; j++)
    {
      bf(i, j).set_z(b(i, j));
    }
  }
}

template <class ZT, class FT> bool MatGSO<ZT, FT>::update_gso_row(int i, int last_j)
{
  givens_row_reduction(i, r_givens.get_cols());
  // FPLLL_TRACE_IN("Updating GSO up to (" << i << ", " << last_j << ")");
  // FPLLL_TRACE("n_known_rows=" << n_known_rows << " n_source_rows=" << n_source_rows);
  if (i >= n_known_rows)
  {
    discover_row();
  }
  FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && last_j >= 0 && last_j < n_source_rows);

  int j = max(0, gso_valid_cols[i]);

  for (; j <= last_j; j++)
  {
    get_gram(ftmp1, i, j);
    FPLLL_DEBUG_CHECK(j == i || gso_valid_cols[j] >= j);
    for (int k = 0; k < j; k++)
    {
      ftmp2.mul(mu(j, k), r(i, k));
      ftmp1.sub(ftmp1, ftmp2);
    }
    r(i, j) = ftmp1;
    if (i > j)
    {
      mu(i, j).div(ftmp1, r(j, j));
      if (!mu(i, j).is_finite())
        return false;
    }
  }

  gso_valid_cols[i] = j;  // = max(0, gso_valid_cols[i], last_j + 1)
  // FPLLL_TRACE_OUT("End of GSO update");


  return true;
}

template <class ZT, class FT> void MatGSO<ZT, FT>::invalidate_gram_row(int i)
{
  for (int j = 0; j <= i; j++)
    gf(i, j).set_nan();
}

template <class ZT, class FT> void MatGSO<ZT, FT>::discover_row()
{
  FPLLL_DEBUG_CHECK(n_known_rows < d);
  /* Early reduction (cols_locked=true) is not allowed when enable_int_gram=true,
     since n_known_cols might be too small to compute all the g(i,j). */
  FPLLL_DEBUG_CHECK(!(cols_locked && enable_int_gram));
  int i = n_known_rows;
  n_known_rows++;
  if (!cols_locked)
  {
    n_source_rows = n_known_rows;
    n_known_cols  = max(n_known_cols, init_row_size[i]);
  }

  if (enable_int_gram)
  {
    for (int j = 0; j <= i; j++)
    {
      dot_product(g(i, j), b[i], b[j], n_known_cols);
    }
  }
  else
  {
    invalidate_gram_row(i);
  }
  gso_valid_cols[i] = 0;
}

template <class ZT, class FT> void MatGSO<ZT, FT>::row_add(int i, int j)
{
  b[i].add(b[j], n_known_cols);
  if (enable_transform)
  {
    u[i].add(u[j]);
    if (enable_inverse_transform)
      u_inv_t[j].sub(u_inv_t[i]);
  }

  if (enable_int_gram)
  {
    // g(i, i) += 2 * g(i, j) + g(j, j)
    ztmp1.mul_2si(g(i, j), 1);
    ztmp1.add(ztmp1, g(j, j));
    g(i, i).add(g(i, i), ztmp1);

    for (int k = 0; k < n_known_rows; k++)
      if (k != i)
        sym_g(i, k).add(sym_g(i, k), sym_g(j, k));
  }
}

template <class ZT, class FT> void MatGSO<ZT, FT>::row_sub(int i, int j)
{
  b[i].sub(b[j], n_known_cols);
  if (enable_transform)
  {
    u[i].sub(u[j]);
    if (enable_inverse_transform)
      u_inv_t[j].add(u_inv_t[i]);
  }

  if (enable_int_gram)
  {
    // g(i, i) += g(j, j) - 2 * g(i, j)
    ztmp1.mul_2si(g(i, j), 1);
    ztmp1.sub(g(j, j), ztmp1);
    g(i, i).add(g(i, i), ztmp1);

    for (int k = 0; k < n_known_rows; k++)
      if (k != i)
        sym_g(i, k).sub(sym_g(i, k), sym_g(j, k));
  }
}

template <class ZT, class FT> void MatGSO<ZT, FT>::row_addmul_si(int i, int j, long x)
{
  b[i].addmul_si(b[j], x, n_known_cols);
  if (enable_transform)
  {
    u[i].addmul_si(u[j], x);
    if (enable_inverse_transform)
      u_inv_t[j].addmul_si(u_inv_t[i], -x);
  }

  if (enable_int_gram)
  {
    /* g(i, i) += 2 * x * g(i, j) +  x^2 * g(j, j)
      (must be done before updating g(i, j)) */
    ztmp1.mul_si(g(i, j), x);
    ztmp1.mul_2si(ztmp1, 1);
    g(i, i).add(g(i, i), ztmp1);
    ztmp1.mul_si(g(j, j), x);
    ztmp1.mul_si(ztmp1, x);
    g(i, i).add(g(i, i), ztmp1);

    // g(i, k) += g(j, k) * x  for k != i
    for (int k = 0; k < n_known_rows; k++)
    {
      if (k == i)
        continue;
      ztmp1.mul_si(sym_g(j, k), x);
      sym_g(i, k).add(sym_g(i, k), ztmp1);
    }
  }
}

template <class ZT, class FT>
void MatGSO<ZT, FT>::row_addmul_si_2exp(int i, int j, long x, long expo)
{
  b[i].addmul_si_2exp(b[j], x, expo, n_known_cols, ztmp1);
  if (enable_transform)
  {
    u[i].addmul_si_2exp(u[j], x, expo, ztmp1);
    if (enable_inverse_transform)
      u_inv_t[j].addmul_si_2exp(u_inv_t[i], -x, expo, ztmp1);
  }

  if (enable_int_gram)
  {
    /* g(i, i) += 2 * (2^e * x) * g(i, j) + 2^(2*e) * x^2 * g(j, j)
      (must be done before updating g(i, j)) */
    ztmp1.mul_si(g(i, j), x);
    ztmp1.mul_2si(ztmp1, expo + 1);
    g(i, i).add(g(i, i), ztmp1);
    ztmp1.mul_si(g(j, j), x);
    ztmp1.mul_si(ztmp1, x);
    ztmp1.mul_2si(ztmp1, 2 * expo);
    g(i, i).add(g(i, i), ztmp1);

    // g(i, k) += g(j, k) * (2^e * x) for k != i
    for (int k = 0; k < n_known_rows; k++)
    {
      if (k == i)
        continue;
      ztmp1.mul_si(sym_g(j, k), x);
      ztmp1.mul_2si(ztmp1, expo);
      sym_g(i, k).add(sym_g(i, k), ztmp1);
    }
  }
}

template <class ZT, class FT>
void MatGSO<ZT, FT>::row_addmul_2exp(int i, int j, const ZT &x, long expo)
{
  b[i].addmul_2exp(b[j], x, expo, n_known_cols, ztmp1);
  if (enable_transform)
  {
    u[i].addmul_2exp(u[j], x, expo, ztmp1);
    if (enable_inverse_transform)
    {
      ZT minus_x;
      minus_x.neg(x);
      u_inv_t[j].addmul_2exp(u_inv_t[i], minus_x, expo, ztmp1);
    }
  }

  if (enable_int_gram)
  {
    /* g(i, i) += 2 * (2^e * x) * g(i, j) + 2^(2*e) * x^2 * g(j, j)
      (must be done before updating g(i, j)) */
    ztmp1.mul(g(i, j), x);
    ztmp1.mul_2si(ztmp1, expo + 1);
    g(i, i).add(g(i, i), ztmp1);
    ztmp1.mul(g(j, j), x);
    ztmp1.mul(ztmp1, x);
    ztmp1.mul_2si(ztmp1, 2 * expo);
    g(i, i).add(g(i, i), ztmp1);

    // g(i, k) += g(j, k) * (2^e * x) for k != i
    for (int k = 0; k < n_known_rows; k++)
    {
      if (k == i)
        continue;
      ztmp1.mul(sym_g(j, k), x);
      ztmp1.mul_2si(ztmp1, expo);
      sym_g(i, k).add(sym_g(i, k), ztmp1);
    }
  }
}

// in row_addmul_we we must have i > j
template <class ZT, class FT>
void MatGSO<ZT, FT>::row_addmul_we(int i, int j, const FT &x, long expo_add)
{
  FPLLL_DEBUG_CHECK(j >= 0 && /* i > j &&*/ i < n_known_rows && j < n_source_rows);

  long expo;
  long lx = x.get_si_exp_we(expo, expo_add);

  if (expo == 0)
  {
    if (lx == 1)
      row_add(i, j);
    else if (lx == -1)
      row_sub(i, j);
    else if (lx != 0)
      row_addmul_si(i, j, lx);
  }
  else if (row_op_force_long)
  {
    row_addmul_si_2exp(i, j, lx, expo);
  }
  else
  {
    x.get_z_exp_we(ztmp2, expo, expo_add);
    row_addmul_2exp(i, j, ztmp2, expo);
  }
}
// In row_swap, i < j
template <class ZT, class FT> void MatGSO<ZT, FT>::row_swap(int i, int j)
{
  FPLLL_DEBUG_CHECK(!enable_inverse_transform);
  b.swap_rows(i, j);
  if (enable_transform)
  {
    u.swap_rows(i, j);
  }

  if (enable_int_gram)
  {
    if (j < i)  // Leaving this check here?
    {
      throw std::runtime_error("Error: in row_swap, i > j, causing errors in the grammatrix.");
    }
    for (int k = 0; k < i; k++)
      g(i, k).swap(g(j, k));
    for (int k = i + 1; k < j; k++)
      g(k, i).swap(g(j, k));
    for (int k = j + 1; k < n_known_rows; k++)
      g(k, i).swap(g(k, j));
    g(i, i).swap(g(j, j));
  }
}

template <class ZT, class FT> void MatGSO<ZT, FT>::move_row(int old_r, int new_r)
{
  FPLLL_DEBUG_CHECK(!cols_locked);
  if (new_r < old_r)
  {
    FPLLL_DEBUG_CHECK(old_r < n_known_rows && !cols_locked);
    for (int i = new_r; i < n_known_rows; i++)
    {
      invalidate_gso_row(i, new_r);
    }
    rotate(gso_valid_cols.begin() + new_r, gso_valid_cols.begin() + old_r,
           gso_valid_cols.begin() + old_r + 1);
    mu.rotate_right(new_r, old_r);
    r.rotate_right(new_r, old_r);
    b.rotate_right(new_r, old_r);
    if (enable_transform)
    {
      u.rotate_right(new_r, old_r);
      if (enable_inverse_transform)
        u_inv_t.rotate_right(new_r, old_r);
    }
    if (enable_int_gram)
    {
      g.rotate_gram_right(new_r, old_r, n_known_rows);
    }
    else
    {
      gf.rotate_gram_right(new_r, old_r, n_known_rows);
      bf.rotate_right(new_r, old_r);
    }
    if (enable_row_expo)
      rotate(row_expo.begin() + new_r, row_expo.begin() + old_r, row_expo.begin() + old_r + 1);
  }
  else if (new_r > old_r)
  {
    for (int i = old_r; i < n_known_rows; i++)
    {
      invalidate_gso_row(i, old_r);
    }
    rotate(gso_valid_cols.begin() + old_r, gso_valid_cols.begin() + old_r + 1,
           gso_valid_cols.begin() + new_r + 1);
    mu.rotate_left(old_r, new_r);
    r.rotate_left(old_r, new_r);
    b.rotate_left(old_r, new_r);
    if (enable_transform)
    {
      u.rotate_left(old_r, new_r);
      if (enable_inverse_transform)
        u_inv_t.rotate_left(old_r, new_r);
    }
    if (enable_int_gram)
    {
      if (old_r < n_known_rows - 1)
      {
        g.rotate_gram_left(old_r, min(new_r, n_known_rows - 1), n_known_rows);
      }
    }
    else
    {
      if (old_r < n_known_rows - 1)
        gf.rotate_gram_left(old_r, min(new_r, n_known_rows - 1), n_known_rows);
      bf.rotate_left(old_r, new_r);
    }
    if (enable_row_expo)
      rotate(row_expo.begin() + old_r, row_expo.begin() + old_r + 1, row_expo.begin() + new_r + 1);
    if (new_r >= n_known_rows)
    {
      rotate(init_row_size.begin() + old_r, init_row_size.begin() + old_r + 1,
             init_row_size.begin() + new_r + 1);
      if (old_r < n_known_rows)
      {
        n_known_rows--;
        n_source_rows        = n_known_rows;
        init_row_size[new_r] = max(b[new_r].size_nz(), 1);
      }
    }
  }
}

template <class ZT, class FT> void MatGSO<ZT, FT>::size_increased()
{
  int old_d = mu.get_rows();

  if (d > alloc_dim)
  {
    if (enable_int_gram)
    {
      g.resize(d, d);
    }
    else
    {
      bf.resize(d, b.get_cols());
      gf.resize(d, d);
    }
    mu.resize(d, d);
    r.resize(d, d);
    r_givens.resize(d, b.get_cols());
    mu_givens.resize(d, b.get_cols());
    gso_valid_cols.resize(d);
    init_row_size.resize(d);
    if (enable_row_expo)
    {
      row_expo.resize(d);
    }
    alloc_dim = d;
  }

  for (int i = old_d; i < d; i++)
  {
    init_row_size[i] = max(b[i].size_nz(), 1);
    if (!enable_int_gram)
    {
      bf[i].fill(0);  // update_bf might not copy all the zeros of b[i]
      update_bf(i);
    }
  }
}

template class MatGSO<Z_NR<long>, FP_NR<double>>;
template class MatGSO<Z_NR<double>, FP_NR<double>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<double>>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class MatGSO<Z_NR<long>, FP_NR<long double>>;
template class MatGSO<Z_NR<double>, FP_NR<long double>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<long double>>;

#endif

#ifdef FPLLL_WITH_QD
template class MatGSO<Z_NR<long>, FP_NR<dd_real>>;
template class MatGSO<Z_NR<double>, FP_NR<dd_real>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<dd_real>>;

template class MatGSO<Z_NR<long>, FP_NR<qd_real>>;
template class MatGSO<Z_NR<double>, FP_NR<qd_real>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<qd_real>>;
#endif

#ifdef FPLLL_WITH_DPE
template class MatGSO<Z_NR<long>, FP_NR<dpe_t>>;
template class MatGSO<Z_NR<double>, FP_NR<dpe_t>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<dpe_t>>;
#endif

template class MatGSO<Z_NR<long>, FP_NR<mpfr_t>>;
template class MatGSO<Z_NR<double>, FP_NR<mpfr_t>>;
template class MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>>;

FPLLL_END_NAMESPACE
