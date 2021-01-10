/// libgp.hpp
/// Marcio Gameiro
/// 
/// Note: include this file for header-only use

#pragma once

#define INLINE_IF_HEADER_ONLY inline 

#include "gp_utils.cc"
#include "cov.cc"
#include "sampleset.cc"
#include "cg.cc"
#include "cov_linear_ard.cc"
#include "cov_linear_one.cc"
#include "cov_matern3_iso.cc"
#include "cov_matern5_iso.cc"
#include "cov_noise.cc"
#include "cov_periodic_matern3_iso.cc"
#include "cov_periodic.cc"
#include "cov_prod.cc"
#include "cov_rq_iso.cc"
#include "cov_se_ard.cc"
#include "cov_se_iso.cc"
#include "cov_sum.cc"
#include "rprop.cc"
#include "input_dim_filter.cc"
#include "cov_factory.cc"
#include "gp.cc"
