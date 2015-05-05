// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dae_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2015-05-05 10:11:29 d3g096
 * 
 * @brief  
 * 
 * 
 */

#include "dae_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class DAESolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// DAESolverImplementation:: constructors / destructor
// -------------------------------------------------------------
DAESolverImplementation::DAESolverImplementation(const parallel::Communicator& comm, 
                                                 const int local_size,
                                                 DAESolverImplementation::JacobianBuilder& jbuilder,
                                                 DAESolverImplementation::FunctionBuilder& fbuilder)
  : parallel::Distributed(comm),
    utility::Configurable("DAESolver"),
    utility::Uncopyable(),
    p_J(comm, local_size, local_size),
    p_Fbuilder(fbuilder), p_Jbuilder(jbuilder)
{
  
}

DAESolverImplementation::~DAESolverImplementation(void)
{
}

// -------------------------------------------------------------
// DAESolverImplementation::p_configure
// -------------------------------------------------------------
void
DAESolverImplementation::p_configure(utility::Configuration::CursorPtr props)
{
  // empty
}


} // namespace math
} // namespace gridpack
