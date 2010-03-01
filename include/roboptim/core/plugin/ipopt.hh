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

#ifndef ROBOPTIM_CORE_IPOPT_HH
# define ROBOPTIM_CORE_IPOPT_HH
# include <roboptim/core/sys.hh>
# include <roboptim/core/portability.hh>

# include <boost/mpl/vector.hpp>

# include <coin/IpSmartPtr.hpp>

# include <roboptim/core/solver.hh>
# include <roboptim/core/twice-derivable-function.hh>


/// \brief Ipopt classes.
namespace Ipopt
{
  class TNLP;
  class IpoptApplication;
} // end of namespace Ipopt


namespace roboptim
{
  namespace detail
  {
    /// \internal
    class MyTNLP;
  };

  /// \addtogroup roboptim_problem
  /// @{

  /// \brief Ipopt based solver.
  ///
  /// Instantiate this solver to solve problems with Ipopt.
  ///
  /// \warning Ipopt needs twice derivable functions, so be sure
  /// to provide hessians in your function's problems.
  class ROBOPTIM_DLLEXPORT IpoptSolver
    : public Solver<TwiceDerivableFunction,
		    boost::mpl::vector<TwiceDerivableFunction> >
  {
  public:
    friend class detail::MyTNLP;

    /// \brief Parent type.
    typedef Solver<TwiceDerivableFunction,
		   boost::mpl::vector<TwiceDerivableFunction> > parent_t;

    /// \brief Instantiate the solver from a problem.
    ///
    /// \param problem problem that will be solved
    explicit IpoptSolver (const problem_t& problem) throw ();

    virtual ~IpoptSolver () throw ();

    /// \brief Solve the problem.
    virtual void solve () throw ();

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

  /// @}
} // end of namespace roboptim
#endif //! ROBOPTIM_CORE_IPOPT_HH
