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

# include <boost/mpl/vector.hpp>

# include <coin/IpSmartPtr.hpp>

# include <roboptim/core/plugin/ipopt/config.hh>
# include <roboptim/core/solver.hh>
# include <roboptim/core/twice-derivable-function.hh>
# include <roboptim/core/plugin/ipopt/ipopt-common.hh>

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
    template <typename T>
    class Tnlp;
  }

  /// \addtogroup roboptim_problem
  /// @{

  /// \brief Ipopt based solver.
  ///
  /// Instantiate this solver to solve problems with Ipopt.
  ///
  /// \warning Ipopt needs twice derivable functions, so be sure
  /// to provide hessians in your function's problems.
  class ROBOPTIM_CORE_PLUGIN_IPOPT_DLLAPI IpoptSolverTd
    : public IpoptSolverCommon<Solver <EigenMatrixDense> >
  {
  public:
    /// \brief RobOptim solver type.
    typedef Solver<EigenMatrixDense> solver_t;

    /// \brief Parent type.
    typedef IpoptSolverCommon<solver_t> parent_t;

    /// \brief Common function type.
    ///
    /// Fuction type which can contain any kind of constraint.
    typedef TwiceDifferentiableFunction commonConstraintFunction_t;

    /// \brief Instantiate the solver from a problem.
    ///
    /// \param problem problem that will be solved
    explicit IpoptSolverTd (const problem_t& problem);

    virtual ~IpoptSolverTd () {}

    template <typename T>
      friend class ::roboptim::detail::Tnlp;
  };

  /// @}
} // end of namespace roboptim
#endif //! ROBOPTIM_CORE_IPOPT_TD_HH
