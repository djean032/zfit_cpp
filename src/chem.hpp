#ifndef CHEM_HPP
#define CHEM_HPP

#include "../utils.hpp"
int solve_rates(CVRhsFn rhs, double tin, double tout,
               std::vector<double> &initial_conditions, sunrealtype params[10]);
int chem_func(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data);

#endif // !CHEM_HPP
