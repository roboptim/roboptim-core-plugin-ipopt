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
#include <typeinfo>

#include <boost/foreach.hpp>
#include <boost/variant/apply_visitor.hpp>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>

#include <roboptim/core/util.hh>

#include "roboptim/core/plugin/ipopt/ipopt-td.hh"
#include "roboptim/core/plugin/ipopt/ipopt-parameters-updater.hh"
#include "ipopt-common.hxx"
#include "tnlp.hh"

namespace roboptim
{
  typedef Solver<TwiceDifferentiableFunction,
		 boost::mpl::vector<LinearFunction,
				    TwiceDifferentiableFunction> >
  solver_ipopt_td_t;

  template class IpoptSolverCommon<solver_ipopt_td_t>;

  // On Microsoft Windows, working with pre-built Ipopt
  // binaries requires the use of the IpoptApplicationFactory.
  IpoptSolverTd::IpoptSolverTd (const problem_t& pb) throw ()
    : parent_t (pb, Ipopt::SmartPtr<Ipopt::TNLP>
		(new detail::Tnlp<IpoptSolverTd> (pb, *this)))
  {
    parameters ()["ipopt.hessian_approximation"].value = "exact";

#ifdef ROBOPTIM_CORE_PLUGIN_IPOPT_VERBOSE
    Ipopt::SmartPtr<Ipopt::Journal> stdout_jrnl =
      getIpoptApplication ()->Jnlst ()->AddFileJournal
      ("console", "stdout", Ipopt::J_ITERSUMMARY);
#endif // ROBPOTIM_CORE_PLUGIN_IPOPT_VERBOSE
  }
} // end of namespace roboptim

extern "C"
{
  using namespace roboptim;
  typedef IpoptSolverTd::solver_t solver_t;

  ROBOPTIM_DLLEXPORT unsigned getSizeOfProblem ();
  ROBOPTIM_DLLEXPORT const char* getTypeIdOfConstraintsList ();
  ROBOPTIM_DLLEXPORT solver_t* create (const IpoptSolverTd::problem_t& pb);
  ROBOPTIM_DLLEXPORT void destroy (solver_t* p);

  ROBOPTIM_DLLEXPORT unsigned getSizeOfProblem ()
  {
    return sizeof (IpoptSolverTd::problem_t);
  }

  ROBOPTIM_DLLEXPORT const char* getTypeIdOfConstraintsList ()
  {
    return typeid (IpoptSolverTd::problem_t::constraintsList_t).name ();
  }

  ROBOPTIM_DLLEXPORT solver_t* create (const IpoptSolverTd::problem_t& pb)
  {
    return new IpoptSolverTd (pb);
  }

  ROBOPTIM_DLLEXPORT void destroy (solver_t* p)
  {
    delete p;
  }
}


// Local Variables:
// compile-command: "cd ../_build && make -k"
// End:
