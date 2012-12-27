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
#include <boost/mpl/vector.hpp>

#include <roboptim/core/solver-factory.hh>

#include "shared-tests/common.hh"
#include "shared-tests/hs071_bis.hh"

using namespace roboptim;
typedef Solver<TwiceDerivableFunction,
	       boost::mpl::vector<TwiceDerivableFunction> > solver_t;

int run_test ()
{
  F f;

  solver_t::problem_t pb (f);
  initialize_problem<solver_t::problem_t,
    roboptim::TwiceDerivableFunction> (pb);

  std::cout << " ############## Displaying the problem ###########" << std::endl;
  std::cout << pb << std::endl;
  std::cout << " ############## End of displaying the problem ###########" << std::endl;

  // Initialize solver
  SolverFactory<solver_t> factory ("ipopt-td", pb);
  solver_t& solver = factory ();

//  solver.parameters()["ipopt.linear_solver"].value = "ma27";
//  solver.parameters()["ipopt.file_print_level"].value = 12;
//  solver.parameters()["ipopt.print_level"].value = 12;

  // Compute the minimum and retrieve the result.
  solver_t::result_t res = solver.minimum ();

  // Display solver information.
  std::cout << solver << std::endl;

  // Check if the minimization has succeed.
  if (res.which () != solver_t::SOLVER_VALUE)
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
