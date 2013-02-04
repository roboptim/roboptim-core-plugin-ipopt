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

#ifndef ROBOPTIM_CORE_IPOPT_TD_HH
# define ROBOPTIM_CORE_IPOPT_TD_HH
# include <roboptim/core/sys.hh>
# include <roboptim/core/portability.hh>

# include <boost/mpl/vector.hpp>

# include <coin/IpSmartPtr.hpp>

# include <roboptim/core/solver.hh>
# include <roboptim/core/twice-derivable-function.hh>
# include <roboptim/core/plugin/ipopt-common.hh>

/// \brief Ipopt classes.
namespace Ipopt
{
  class IpoptApplication;
} // end of namespace Ipopt


namespace roboptim
{
  namespace detail
  {
    /// \internal
    class TnlpTd;
  }

  /// \addtogroup roboptim_problem
  /// @{

  /// \brief Ipopt based solver.
  ///
  /// Instantiate this solver to solve problems with Ipopt.
  ///
  /// \warning Ipopt needs twice derivable functions, so be sure
  /// to provide hessians in your function's problems.
  class ROBOPTIM_DLLEXPORT IpoptSolverTd
    : public IpoptSolverCommon<
    Solver <TwiceDifferentiableFunction,
	    boost::mpl::vector<LinearFunction, TwiceDifferentiableFunction> > >
  {
  public:
    friend class detail::TnlpTd;

    /// \brief RobOptim solver type.
    typedef Solver<TwiceDifferentiableFunction,
      boost::mpl::vector<LinearFunction,
			 TwiceDifferentiableFunction> > solver_t;

    /// \brief Parent type.
    typedef IpoptSolverCommon<solver_t> parent_t;

    /// \brief Instantiate the solver from a problem.
    ///
    /// \param problem problem that will be solved
    explicit IpoptSolverTd (const problem_t& problem) throw ();

    virtual ~IpoptSolverTd () throw () {}
  };

  /// @}
} // end of namespace roboptim
#endif //! ROBOPTIM_CORE_IPOPT_TD_HH
