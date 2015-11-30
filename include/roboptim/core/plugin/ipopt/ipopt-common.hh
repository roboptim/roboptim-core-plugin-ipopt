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

#ifndef ROBOPTIM_CORE_IPOPT_COMMON_HH
# define ROBOPTIM_CORE_IPOPT_COMMON_HH
# include <roboptim/core/sys.hh>
# include <roboptim/core/portability.hh>

# include <boost/mpl/vector.hpp>
# include <boost/optional.hpp>

# include <roboptim/core/solver.hh>
# include <roboptim/core/linear-function.hh>
# include <roboptim/core/differentiable-function.hh>
# include <roboptim/core/twice-differentiable-function.hh>

/// \brief Ipopt classes.
namespace Ipopt
{
  class IpoptApplication;
  class TNLP;
  class IpoptData;
  class IpoptCalculatedQuantities;
  template <typename T> class SmartPtr;
} // end of namespace Ipopt


namespace roboptim
{
  /// \addtogroup roboptim_problem
  /// @{

  template <typename T>
  class IpoptSolverCommon;

  /// \brief Ipopt common solver.
  ///
  /// This solver shares common piece of code of the two solvers.
  template <typename T>
  class IpoptSolverCommon : public T
  {
  public:
    /// \brief Categorize constraints.
    ///
    /// Used with the which method of the Boost.Variant, it
    /// allows to check for a constraint's real type.
    ///
    /// \warning Make sure to keep enum values in the
    /// same order than the MPL vector used to specify CLIST.
    enum ConstraintType
      {
	/// \brief Constraint is a linear function.
	LINEAR = 0,
	/// \brief Constraint is a differentiable or a twice
	/// differentiable function depending on the solve.
	NONLINEAR = 1
      };

    /// \brief Parent type.
    typedef T parent_t;

    typedef typename T::callback_t callback_t;
    typedef typename T::problem_t problem_t;

    /// \brief Instantiate the solver from a problem.
    ///
    /// \param pb problem that will be solved.
    /// \param tnlp pointer to Ipopt solver.
    explicit IpoptSolverCommon (const problem_t& pb,
				Ipopt::SmartPtr<Ipopt::TNLP> tnlp);

    virtual ~IpoptSolverCommon ();

    /// \brief Solve the problem.
    void solve ();

    /// \brief Get Ipopt Application object for Ipopt specific tuning.
    ///
    /// Consult Ipopt documentation for information regarding
    /// IpoptApplication class.
    virtual Ipopt::SmartPtr<Ipopt::IpoptApplication> getIpoptApplication ();

    virtual void
    setIterationCallback (callback_t callback)
    {
      callback_ = callback;
    }

    const callback_t& callback () const
    {
      return callback_;
    }
  private:
    /// \brief Initialize parameters.
    ///
    /// Add solver parameters. Called during construction.
    void initializeParameters ();

    /// \brief Read parameters and update associated options in Ipopt.
    ///
    /// Called before solving problem.
    void updateParameters ();

    /// \brief Smart pointer to the Ipopt non linear problem description.
    Ipopt::SmartPtr<Ipopt::TNLP> nlp_;
    /// \brief Smart pointer to the Ipopt application instance.
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app_;

    /// \brief Intermediate callback (called at each end
    /// of iteration).
    callback_t callback_;

    /// \brief Number of times the solve method was called.
    int solveCounter_;
  };

  /// @}
} // end of namespace roboptim

#endif //! ROBOPTIM_CORE_IPOPT_COMMON_HH
