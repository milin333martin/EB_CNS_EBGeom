#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>
#include <cmath>
#include <algorithm>
#include <cns_prob.H>

// Forward declaration - the actual implementation is in STLLoader.cpp
void initialize_EB2(const amrex::Geometry& geom, const int required_level, const int max_coarsening_level);
