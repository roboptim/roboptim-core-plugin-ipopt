// Copyright (C) 2009 by Thomas Moulard, AIST, CNRS, INRIA.
//
// This file is part of the roboptim.
//
// roboptim is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// roboptim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with roboptim.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/variant/get.hpp>

#include <coin/IpIpoptApplication.hpp>

#include <roboptim/core/plugin/ipopt-td.hh>
#include "ipopt-common.hxx"

#include "shared-tests/common.hh"
#include "shared-tests/hs071.hh"

int run_test ()
{
  F f;

  IpoptSolverTd::problem_t pb (f);
  initialize_problem<IpoptSolverTd::problem_t,
    roboptim::TwiceDerivableFunction> (pb);

  // Initialize solver
  IpoptSolverTd solver (pb);

  Ipopt::SmartPtr<Ipopt::Journal> stdout_jrnl =
    solver.getIpoptApplication ()->Jnlst ()->AddFileJournal
    ("console", "stdout", Ipopt::J_ITERSUMMARY);

  // Compute the minimum and retrieve the result.
  IpoptSolverTd::result_t res = solver.minimum ();

  // Display solver information.
  std::cout << solver << std::endl;

  // Check if the minimization has succeed.
  if (res.which () != IpoptSolverTd::SOLVER_VALUE)
    {
      std::cout << "A solution should have been found. Failing..."
                << std::endl
                << boost::get<SolverError> (res).what ()
                << std::endl;
      return 1;
    }

  // Get the result.
  Result& result = boost::get<Result> (res);

  // Display the result.
  std::cout << "A solution has been found: " << std::endl;
  std::cout << result << std::endl;
  return 0;
}


GENERATE_TEST ()
