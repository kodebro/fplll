#ifndef BKZ_PARAMS_H
#define BKZ_PARAMS_H

/* (C) 2014-2016 Martin Albrecht.

   This file is part of fplll. fplll is free software: you can
   redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software
   Foundation, either version 2.1 of the License, or (at your option)
   any later version.

   fplll is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with fplll. If not, see
   <http://www.gnu.org/licenses/>.

*/

#include <string>
#include <vector>
#include "defs.h"

FPLLL_BEGIN_NAMESPACE

class Pruning {

public:

  double radius_factor;              //< radius/Gaussian heuristic
  double probability;               //< success probability
  std::vector<double> coefficients;  //< pruning coefficients

  Pruning() : radius_factor(1.), probability(1.) {};

  /** Sets all pruning coefficients to 1, except the last <level>
      coefficients, these will be linearly with slope -1 / blockSize.

      @param level number of levels in linear descent
  */

  static Pruning LinearPruning(int block_size, int level) {

    Pruning pruning = Pruning();
    int start_descent = block_size - level;

    if (start_descent > block_size)
      start_descent = block_size;

    if (start_descent < 1)
      start_descent = 1;

    pruning.coefficients.resize(block_size);
    for (int k = 0; k < start_descent; k++) {
      pruning.coefficients[k] = 1.0;
    }
    for (int k = 0; k < block_size - start_descent; k++) {
      pruning.coefficients[start_descent + k] = ((double)(block_size - k - 1)) / block_size;
    }
    // TODO: need to adapt probability
    pruning.radius_factor = 1.0;
    pruning.probability = 1.0;

    return pruning;
  }
};


class Strategy {
public:
  vector<Pruning> pruning_parameters;
  vector<int> preprocessing_blocksizes;

  static Strategy EmptyStrategy() {
    Strategy strat;
    strat.pruning_parameters.emplace_back(Pruning());
    return strat;
  };

  const Pruning &get_pruning(double radius, double gh) const;
};


class BKZParam {
public:
  BKZParam(int block_size, vector<Strategy> &strategies,
           double delta = LLL_DEF_DELTA, int flags = BKZ_DEFAULT,
           int max_loops = 0, double max_time = 0, double auto_abort_scale = 1.0,
           int auto_abort_max_no_dec = 5, double gh_factor = 1.1,
           double min_success_probability = 0.5,
           int rerandomization_density = 3)
    : block_size(block_size), strategies(strategies), delta(delta), flags(flags), max_loops(max_loops), max_time(max_time),
      auto_abort_scale(auto_abort_scale), auto_abort_max_no_dec(auto_abort_max_no_dec),
      gh_factor(gh_factor), dump_gso_filename("gso.log"),
      min_success_probability(min_success_probability),
      rerandomization_density(rerandomization_density) {
    if (strategies.empty()) {
      strategies = vector<Strategy>();
      for (long b = 0; b <= block_size; ++b) {
        strategies.emplace_back(std::move(Strategy::EmptyStrategy()));
      }
    }
  };

  /** Block size used for enumeration **/
  int block_size;

  /** Strategies (pruning coefficients, preprocessing)  */

  vector<Strategy> &strategies;

  /** LLL parameter delta **/
  double delta;

  /** See BKZFlags **/
  int flags;

  /** Maximum number of loops to execute **/
  int max_loops;

  /** Maximum time to spend **/
  double max_time;

  /** If BKZ_AUTOABORT is set, We abort if newSlope < autoAbort_scale * oldSlope
      is true for autoAbort_maxNoDec loops.
   */
  double auto_abort_scale;
  int    auto_abort_max_no_dec;

  /** If BKZ_GH_BND is set, the enumeration bound will be set to ghFactor times
      the Gaussian Heuristic
  */

  double gh_factor;

  /** If BKZ_DUMP_GSO is set, the norms of the GSO matrix are written to this
      file after each complete round.
  */

  string dump_gso_filename;

  /** minimum success probability when using extreme pruning */

  double min_success_probability;

  /** density of rerandomization operation when using extreme pruning **/

  int rerandomization_density;

};

FPLLL_END_NAMESPACE
#endif /* BKZ_PARAMS_H */
