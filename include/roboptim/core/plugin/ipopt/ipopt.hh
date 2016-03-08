// Copyright (C) 2010 by Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_CORE_IPOPT_HH
# define ROBOPTIM_CORE_IPOPT_HH

# include <boost/mpl/vector.hpp>

# include <coin/IpSmartPtr.hpp>
# include <coin/IpReturnCodes.hpp> // for AlgorithmMode

# include <roboptim/core/plugin/ipopt/config.hh>
# include <roboptim/core/solver.hh>
# include <roboptim/core/derivable-function.hh>
# include <roboptim/core/plugin/ipopt/ipopt-common.hh>

/// \brief Ipopt classes.
namespace Ipopt
{
  class IpoptApplication;
  class IpoptData;
  class IpoptCalculatedQuantities;
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
  class ROBOPTIM_CORE_PLUGIN_IPOPT_DLLAPI IpoptSolver
    : public IpoptSolverCommon<Solver<EigenMatrixDense> >
  {
  public:
    /// \brief RobOptim solver type.
    typedef Solver<EigenMatrixDense> solver_t;

    /// \brief Parent type.
    typedef IpoptSolverCommon<solver_t> parent_t;

    /// \brief Common function type.
    ///
    /// Fuction type which can contain any kind of constraint.
    typedef DifferentiableFunction commonConstraintFunction_t;

    /// \brief Instantiate the solver from a problem.
    ///
    /// \param problem problem that will be solved
    explicit IpoptSolver (const problem_t& problem);

    virtual ~IpoptSolver () {}

    template <typename T>
      friend class ::roboptim::detail::Tnlp;
  };

  /// @}
} // end of namespace roboptim

#endif //! ROBOPTIM_CORE_IPOPT_HH
