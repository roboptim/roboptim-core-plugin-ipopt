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

#include <cassert>
#include <cstring>

#include <boost/foreach.hpp>
#include <boost/variant/apply_visitor.hpp>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>

#include <roboptim/core/util.hh>

#include "roboptim/core/plugin/ipopt.hh"
#include "roboptim/core/plugin/ipopt-parameters-updater.hh"
#include "ipopt-common.hxx"
#include "tnlp.hh"

namespace roboptim
{
  typedef Solver<DifferentiableFunction,
		 boost::mpl::vector<LinearFunction, DifferentiableFunction> >
  ipopt_solver_t;

  template class IpoptSolverCommon<ipopt_solver_t>;

  IpoptSolver::IpoptSolver (const problem_t& pb) throw ()
    : parent_t (pb, Ipopt::SmartPtr<Ipopt::TNLP>
		(new detail::Tnlp<IpoptSolver> (pb, *this)))

  {
    parameters ()["ipopt.hessian_approximation"].value = "limited-memory";
  }
} // end of namespace roboptim

extern "C"
{
  using namespace roboptim;
  typedef IpoptSolver::solver_t solver_t;

  ROBOPTIM_DLLEXPORT unsigned getSizeOfProblem ();
  ROBOPTIM_DLLEXPORT solver_t* create (const IpoptSolver::problem_t& pb);
  ROBOPTIM_DLLEXPORT void destroy (solver_t* p);

  ROBOPTIM_DLLEXPORT unsigned getSizeOfProblem ()
  {
    return sizeof (IpoptSolver::problem_t);
  }

  ROBOPTIM_DLLEXPORT solver_t* create (const IpoptSolver::problem_t& pb)
  {
    return new IpoptSolver (pb);
  }

  ROBOPTIM_DLLEXPORT void destroy (solver_t* p)
  {
    delete p;
  }
}


// Local Variables:
// compile-command: "make -k -C ../_build"
// End:
