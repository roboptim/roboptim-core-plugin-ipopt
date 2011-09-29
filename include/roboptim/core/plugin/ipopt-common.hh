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

# include <roboptim/core/solver.hh>
# include <roboptim/core/derivable-function.hh>
# include <roboptim/core/twice-derivable-function.hh>

/// \brief Ipopt classes.
namespace Ipopt
{
  class IpoptApplication;
  class TNLP;
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
    /// \brief Parent type.
    typedef T parent_t;

    typedef typename T::problem_t problem_t;

    /// \brief Instantiate the solver from a problem.
    ///
    /// \param problem problem that will be solved
    explicit IpoptSolverCommon (const problem_t& pb,
				Ipopt::SmartPtr<Ipopt::TNLP> tnlp) throw ();

    virtual ~IpoptSolverCommon () throw ();

    /// \brief Solve the problem.
    void solve () throw ();

    /// \brief Get Ipopt Application object for Ipopt specific tuning.
    ///
    /// Consult Ipopt documentation for information regarding
    /// IpoptApplication class.
    virtual Ipopt::SmartPtr<Ipopt::IpoptApplication> getIpoptApplication ()
      throw ();
  private:

    /// \brief Initialize parameters.
    ///
    /// Add solver parameters. Called during construction.
    void initializeParameters () throw ();

    /// \brief Read parameters and update associated options in Ipopt.
    ///
    /// Called before solving problem.
    void updateParameters () throw ();

    /// \brief Smart pointer to the Ipopt non linear problem description.
    Ipopt::SmartPtr<Ipopt::TNLP> nlp_;
    /// \brief Smart pointer to the Ipopt application instance.
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app_;
  };

extern template class IpoptSolverCommon<
  Solver<DifferentiableFunction,
	 boost::mpl::vector<DifferentiableFunction> > >;

extern template class IpoptSolverCommon<
  Solver<TwiceDifferentiableFunction,
	 boost::mpl::vector<TwiceDifferentiableFunction> > >;

  /// @}
} // end of namespace roboptim

#endif //! ROBOPTIM_CORE_IPOPT_COMMON_HH
