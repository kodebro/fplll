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

#include "gso_givens.h"

FPLLL_BEGIN_NAMESPACE

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::initialize_l_givens_matrix()
{
  mu.fill(0.0); // Empty mu matrix

  if (enable_row_expo)
  {
    throw std::runtime_error("Error: givens rotations are not yet implemented for enable_row_expo is true.");
  }
  else
  {
    /*if (l_givens.get_rows() != d)
    {
      throw std::runtime_error("Error: l_givens does not have good dimensions.");
    }
    if (l_givens.get_cols() != b.get_cols())
    {
      throw std::runtime_error("Error: l_givens does not have good dimensions.");
    }*/
    for (int i = 0; i < d; i++)
    {
      for (int j = 0; j < l_givens.get_cols(); j++)
      {
        l_givens(i, j).set_z(b(i, j));
      }
    }
  }
}


// The givens rotation introduces a zero in place  (k,j) in the matrix, and puts the residue
// on place (k,i). Most of the cases, j = i-1
//              j|i
//   _________________________
//  |						  |    
//  |_________________________|
//k |__________|x|y|__________|
//  |						  |
//  |						  |
//  |_________________________|
//
//        TRAMSFORMS TO
//              
//
//              j|i
//   _________________________
//  |		    * *	     	  |    
//  |___________*_*___________|
//k |__________|z|0|__________|
//  |			* *			  |
//  |		    * *		      |
//  |-------------------------|
//
//  The stars indicate that the values
//  there can be changed [in linear combinations]. 
// The rest of  the matrix is fixed.
// 


template <class ZT, class FT> void MatGSOGivens<ZT, FT>::givens_rotation(int col_i, int col_j, int row_k)
{
  // TODO: This maybe can be sped up if we manage to remove some FT 's
  // we might want to have some more 'tmp1';


  // If the place where we want to introduce
  // a zero is already zero, we do nothing.
  if (l_givens(row_k,col_j).is_zero()) {
    return;
  }

  // Computes the 'c' and 's' of the givens, 
  // by means of the hypot-function,
  // which hopefully avoids overflow.
  FT c, s;
  ftmp1.hypot(l_givens(row_k, col_i), l_givens(row_k, col_j));
  c.div(l_givens(row_k, col_i), ftmp1);
  s.div(l_givens(row_k, col_j), ftmp1);


  // Here we take multiples of columns
  // -  Remark that we leave alone the first
  //    row_k - 1 entries of the columns.
  //    For effectivity this is obviously better,
  //    but for numerical stability I don't know.
  // -  Also, we could force the entry that should be
  //    zero as 'zero'. Now know whether this is good or not.
  //
  // -> Answer to both -> Probably no effect on stability.

  for (int k = row_k; k < l_givens.get_rows(); k++)
  {
    ftmp1 = l_givens(k, col_i); 
    ftmp2 = l_givens(k, col_j); 
    l_givens(k, col_i).mul(ftmp1, c); // r_(k,col_i) = c*r_(k,col_i) + s*r_(k,col_j)
    l_givens(k, col_i).addmul(ftmp2, s);

    l_givens(k, col_j).neg(s);
    l_givens(k, col_j).mul(ftmp1, l_givens(k, col_j));
    l_givens(k, col_j).addmul(ftmp2, c); // r_(k,col_j) = -s*r_(k,col_i) + c*r_(k,col_j)
  }
  // "Forcing" zero 
  l_givens(row_k, col_j) = 0.0;
}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::compute_mu_and_r_columns(int starting_column, int last_column)
{
  for(int col = starting_column; col <= last_column; col++) {

    for (int i = col; i < mu.get_rows(); i++)
      mu(i,col).div(l_givens(i,col),l_givens(col,col));

    ftmp1 = l_givens(col,col);
    for (int i = col; i < l_givens.get_rows(); i++)
      r(i, col).mul(ftmp1, l_givens(i, col));  

  }
}

// The givens_row_reduction introduces a zero-sequence from the diagonal of row_k
// to the 'rightmost_nonzero_entry' column.
// 
//              k|        r
//   _________________________
//  |						  |    
//  |_________________________|
//k |__________|a|b|c|d|e|f|__|
//  |						  |
//  |						  |
//  |_________________________|
//
//        TRAMSFORMS TO
//              
//
//              k|        r
//   _________________________
//  |			* * * * * *   |    
//  |_________________________|
//k |__________|a|0|0|0|0|0|__|
//  |			* * * * * *   |
//  |			* * * * * *   |
//  |-------------------------|
//
// The stars indicate that the entries
// there might be changed into linear
// combinations of their neighbour-entries.
// The rest of the matrix is fixed.

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::givens_row_reduction(int row_k, int rightmost_nonzero_entry)
{
  // - Could be improved by pushing the zero to the rightmost non-zero entry.
  // - Or, just pushing the zero to the diagonal.
  for (int i = rightmost_nonzero_entry; i > row_k; i--)
    givens_rotation(i - 1, i, row_k);


  compute_mu_and_r_column(row_k);
}


template <class ZT, class FT> void MatGSOGivens<ZT, FT>::clean_mu()
{

  for(int i = 0; i < mu.get_rows(); i++) {
    mu(i,i) = 1.0;
    for(int j = i + 1; j < mu.get_cols(); j++) {
      mu(i,j) = 0.0;  
    }
  }

}


template <class ZT, class FT> void MatGSOGivens<ZT, FT>::update_bf(int i)
{
   // Does nothing
   // Maybe recompute the L-matrix of the givens?
}


template <class ZT, class FT> bool MatGSOGivens<ZT, FT>::update_gso_row(int i, int last_j)
{
  givens_row_reduction(i, l_givens.get_cols()-1);
  
  // FPLLL_TRACE_IN("Updating GSO up to (" << i << ", " << last_j << ")");
  // FPLLL_TRACE("n_known_rows=" << n_known_rows << " n_source_rows=" << n_source_rows);
  if (i >= n_known_rows)
  {
    discover_row();
  }
  //FPLLL_DEBUG_CHECK(i >= 0 && i < n_known_rows && last_j >= 0 && last_j < n_source_rows);

  // TODO, Maybe add some functionality here.

  return true;
}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::invalidate_gram_row(int i)
{
	// TODO maybe some functionality later.
}

// TODO Can discover_row be deleted?
template <class ZT, class FT> void MatGSOGivens<ZT, FT>::discover_row()
{
  FPLLL_DEBUG_CHECK(n_known_rows < d);
  FPLLL_DEBUG_CHECK(!(cols_locked));

}

// Givens ready.
template <class ZT, class FT> void MatGSOGivens<ZT, FT>::row_add(int i, int j)
{
 
	  b[i].add(b[j], n_known_cols);
	  if (enable_transform
	)  {
	    u[i].add(u[j]);
	    if (enable_inverse_transform)
	      u_inv_t[j].sub(u_inv_t[i]);
	  }


      // IS THIS STABLE?
    if (j < i) {
      // Doing b_i <- b_i + c b_j doesn't affect the triangularity

      l_givens[i].add(l_givens[j],j+1);


      // TODO: making this lazy

      for(int k = 0; k <= j; k++) {
        mu[i][k].div(l_givens[i][k],l_givens(k,k));
        r[i][k].mul(l_givens[i][k],l_givens(k,k));

      }


      // No other recomputation of mu or r needed.

     } else {
      // i < j, so affects triangularity 
      l_givens[i].add(l_givens[j],j+1);

      givens_row_reduction(i,j);
      for(int k = i+1; k <= j; k++) {
        givens_rotation(k,k+1,k);
      }

      // TODO making this lazy
      for(int k = 0; k < i; k++) {
        mu[i][k].div(l_givens[i][k],l_givens(k,k));
        r[i][k].mul(l_givens[i][k],l_givens(k,k));
      }
      compute_mu_and_r_columns(i,j);
    }
  
}

// Givens ready
template <class ZT, class FT> void MatGSOGivens<ZT, FT>::row_sub(int i, int j)
{

  b[i].sub(b[j], n_known_cols);
  if (enable_transform)
  {
    u[i].sub(u[j]);
    if (enable_inverse_transform)
      u_inv_t[j].add(u_inv_t[i]);
  }

      // IS THIS STABLE?
    if (j < i) {
      // Doing b_i <- b_i + c b_j doesn't affect the triangularity

      l_givens[i].sub(l_givens[j],j+1);


      // TODO: making this lazy

      for(int k = 0; k <= j; k++) {
        mu[i][k].div(l_givens[i][k],l_givens(k,k));
        r[i][k].mul(l_givens[i][k],l_givens(k,k));

      }


      // No other recomputation of mu or r needed.

     } else {
      // i < j, so affects triangularity 
      l_givens[i].sub(l_givens[j],j+1);

      givens_row_reduction(i,j);
      for(int k = i+1; k <= j; k++) {
        givens_rotation(k,k+1,k);
      }

      // TODO making this lazy
      for(int k = 0; k < i; k++) {
        mu[i][k].div(l_givens[i][k],l_givens(k,k));
        r[i][k].mul(l_givens[i][k],l_givens(k,k));
      }
      compute_mu_and_r_columns(i,j);
    }
}

// Givens-ready
template <class ZT, class FT> void MatGSOGivens<ZT, FT>::row_addmul_si(int i, int j, long x)
{


  // TODO addmul_si not possible, because
  // l_givens is of type NumVect
  // and not of MatrixRow
  // Is this well-resolved??

  b[i].addmul_si(b[j], x, n_known_cols);

  if (enable_transform)
  {
    u[i].addmul_si(u[j], x);
    if (enable_inverse_transform)
      u_inv_t[j].addmul_si(u_inv_t[i], -x);
  }

      // IS THIS STABLE?
    if (j < i) {
      // Doing b_i <- b_i + c b_j doesn't affect the 

      // TODO Changing to addmul okay?
      l_givens[i].addmul(l_givens[j],x,j+1);
      // TODO: making this lazy

      for(int k = 0; k <= j; k++) {
        mu[i][k].div(l_givens[i][k],l_givens(k,k));
        r[i][k].mul(l_givens[i][k],l_givens(k,k));

      }


      // No other recomputation of mu or r needed.

     } else {
      // i < j, so affects triangularity 
      l_givens[i].addmul(l_givens[j],x,j+1);

      givens_row_reduction(i,j);
      for(int k = i+1; k <= j; k++) {
        givens_rotation(k,k+1,k);
      }

      // TODO making this lazy
      for(int k = 0; k < i; k++) {
        mu[i][k].div(l_givens[i][k],l_givens(k,k));
        r[i][k].mul(l_givens[i][k],l_givens(k,k));
      }
      compute_mu_and_r_columns(i,j);
    }




}


// TODO Still needs to be implemented!
template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::row_addmul_si_2exp(int i, int j, long x, long expo)
{
  throw std::runtime_error("Error: exponents are not yet implemented for givens rotations");
  b[i].addmul_si_2exp(b[j], x, expo, n_known_cols, ztmp1);
  if (enable_transform)
  {
    u[i].addmul_si_2exp(u[j], x, expo, ztmp1);
    if (enable_inverse_transform)
      u_inv_t[j].addmul_si_2exp(u_inv_t[i], -x, expo, ztmp1);
  }

// Exponents to be implemented later
/*
      // IS THIS STABLE?
    if (j < i) {
      // Doing b_i <- b_i + c b_j doesn't affect the triangularity

      l_givens[i].addmul_2exp(l_givens[j],x,expo, ftmp1);


      // TODO: making this lazy

      for(int k = 0; k <= j; k++) {
        mu[i][k].div(l_givens[i][k],l_givens(k,k));
        r[i][k].mul(l_givens[i][k],l_givens(k,k));

      }

      // No other recomputation of mu or r needed.

     } else {
      // i < j, so affects triangularity 
      l_givens[i].addmul_2exp(l_givens[j],x,expo, ftmp1);

      givens_row_reduction(i,j);
      for(int k = i+1; k <= j; k++) {
        givens_rotation(k,k+1,k);
      }

      // TODO making this lazy
      for(int k = 0; k < i; k++) {
        mu[i][k].div(l_givens[i][k],l_givens(k,k));
        r[i][k].mul(l_givens[i][k],l_givens(k,k));
      }
      compute_mu_and_r_columns(i,j);
    }

*/
}


// Still needs to be implemented!!
template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::row_addmul_2exp(int i, int j, const ZT &x, long expo)
{
  throw std::runtime_error("Error: exponents are not yet implemented for givens rotations");
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

// Exponents to be implemented later
/*
  FT tmpx;
  tmpx.set_z(x);

      // IS THIS STABLE?
    if (j < i) {
      // Doing b_i <- b_i + c b_j doesn't affect the triangularity

        l_givens[i].addmul_2exp(l_givens[j], tmpx, expo, j+1, ftmp1);


      // TODO: making this lazy

      for(int k = 0; k <= j; k++) {
        mu[i][k].div(l_givens[i][k],l_givens(k,k));
        r[i][k].mul(l_givens[i][k],l_givens(k,k));

      }

      // No other recomputation of mu or r needed.

     } else {
      // i < j, so affects triangularity 
      l_givens[i].addmul_2exp(l_givens[j],x,expo, ftmp1);

      givens_row_reduction(i,j);
      for(int k = i+1; k <= j; k++) {
        givens_rotation(k,k+1,k);
      }

      // TODO making this lazy
      for(int k = 0; k < i; k++) {
        mu[i][k].div(l_givens[i][k],l_givens(k,k));
        r[i][k].mul(l_givens[i][k],l_givens(k,k));
      }
      compute_mu_and_r_columns(i,j);
    }
*/


}

// in row_addmul_we we must have i > j
template <class ZT, class FT>
void MatGSOGivens<ZT, FT>::row_addmul_we(int i, int j, const FT &x, long expo_add)
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
// Is Givens-ready.
template <class ZT, class FT> void MatGSOGivens<ZT, FT>::row_swap(int i, int j)
{

  FPLLL_DEBUG_CHECK(!enable_inverse_transform);
   if (j < i)  
   {
   	int k = i;
   	i = j;
   	j = k;
   }

  // *******************
  // Swap the rows of b
  //********************
  b.swap_rows(i, j);
  if (enable_transform)
  {
    u.swap_rows(i, j);
  }


  // *****************
  // Givens equivalent
  // *****************
  l_givens.swap_rows(i,j); 
  mu.swap_rows(i,j);
  r.swap_rows(i,j);

  givens_row_reduction(i,j);
  for(int k = i+1; k <= j; k++) {
  	givens_rotation(k,k+1,k);
  }
  compute_mu_and_r_columns(i,j);


}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::move_row(int old_r, int new_r)
{
  FPLLL_DEBUG_CHECK(!cols_locked);
  if (new_r < old_r)
  {

    b.rotate_right(new_r, old_r);

    if (enable_transform)
    {
      u.rotate_right(new_r, old_r);
      if (enable_inverse_transform)
        u_inv_t.rotate_right(new_r, old_r);
    }

    // Disturbs triangularity a bit like swap.
    // So, we correct the givens a bit.
    l_givens.rotate_right(new_r,old_r); 
    mu.rotate_right(new_r,old_r);
    r.rotate_right(new_r,old_r);

    givens_row_reduction(new_r,old_r+1);
    for(int k = new_r+1; k <= old_r; k++) {
  	  givens_rotation(k,k+1,k);
    }
    compute_mu_and_r_columns(new_r,old_r);

  }
  else if (new_r > old_r)
  {


    b.rotate_left(old_r, new_r);
    if (enable_transform)
    {
      u.rotate_left(old_r, new_r);
      if (enable_inverse_transform)
        u_inv_t.rotate_left(old_r, new_r);
    }



    // Disturbs triangularity a bit like swap.
    // So, we correct the givens a bit.
    l_givens.rotate_left(old_r, new_r);
    mu.rotate_left(old_r, new_r);
    r.rotate_left(old_r, new_r);


    for(int k = old_r; k <= new_r-1; k++) { // Maybe change new_r - 1 to new_r 
  	  givens_rotation(k,k+1,k);
    }
    compute_mu_and_r_columns(old_r,new_r);

  }
}

template <class ZT, class FT> void MatGSOGivens<ZT, FT>::size_increased()
{

  if (d > alloc_dim)
  {

      mu.resize(d, b.get_cols());
      l_givens.resize(d, b.get_cols());
   	  r.resize(d, b.get_cols());   		
  }
  
}

template class MatGSOGivens<Z_NR<long>, FP_NR<double>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<double>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<double>>;

#ifdef FPLLL_WITH_LONG_DOUBLE
template class MatGSOGivens<Z_NR<long>, FP_NR<long double>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<long double>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<long double>>;

#endif

#ifdef FPLLL_WITH_QD
template class MatGSOGivens<Z_NR<long>, FP_NR<dd_real>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<dd_real>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<dd_real>>;

template class MatGSOGivens<Z_NR<long>, FP_NR<qd_real>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<qd_real>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<qd_real>>;
#endif

#ifdef FPLLL_WITH_DPE
template class MatGSOGivens<Z_NR<long>, FP_NR<dpe_t>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<dpe_t>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<dpe_t>>;
#endif

template class MatGSOGivens<Z_NR<long>, FP_NR<mpfr_t>>;
template class MatGSOGivens<Z_NR<double>, FP_NR<mpfr_t>>;
template class MatGSOGivens<Z_NR<mpz_t>, FP_NR<mpfr_t>>;

FPLLL_END_NAMESPACE
