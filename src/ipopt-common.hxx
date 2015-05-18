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

#ifndef ROBOPTIM_CORE_IPOPT_COMMON_HXX
# define ROBOPTIM_CORE_IPOPT_COMMON_HXX

# include <roboptim/core/sys.hh>
# include <roboptim/core/portability.hh>

# include <stdexcept>

# include <boost/mpl/vector.hpp>

# include <coin/IpSmartPtr.hpp>
# include <coin/IpIpoptApplication.hpp>

# include "roboptim/core/plugin/ipopt/ipopt-parameters-updater.hh"
# include "roboptim/core/plugin/ipopt/ipopt-common.hh"

# ifndef IPOPT_DEFAULT_LINEAR_SOLVER
   // Enable by default MUMPS which is the only open-source
   // solver provided by Ipopt.
#  define IPOPT_DEFAULT_LINEAR_SOLVER "mumps"
# endif //! IPOPT_DEFAULT_LINEAR_SOLVER

namespace roboptim
{
  template<typename T>
  IpoptSolverCommon<T>::
  IpoptSolverCommon (const problem_t& pb,
		     Ipopt::SmartPtr<Ipopt::TNLP> tnlp)
    : parent_t (pb),
      nlp_ (tnlp),
      app_ (IpoptApplicationFactory ()),
      callback_ ()
  {
    app_->Jnlst()->DeleteAllJournals();

    // Re-throw non-Ipopt exceptions
    app_->RethrowNonIpoptException (true);

    // Initialize parameters.
    initializeParameters ();
  }

  template<typename T>
  IpoptSolverCommon<T>::
  ~IpoptSolverCommon ()
  {}

  // /!\ this->result_ is filled by tnlp.hxx, do not overwrite!
#define SWITCH_ERROR(NAME, ERROR)		\
  case NAME:					\
  break

#define SWITCH_FATAL(NAME, ERROR)		\
  case NAME:					\
  throw std::runtime_error (ERROR);		\
  break

#define SWITCH_OK(NAME, CASES)			\
  case NAME:					\
  {						\
    int status = app_->OptimizeTNLP (nlp_);	\
    switch (status)				\
      {						\
	CASES;					\
      }						\
  }						\
  break

#define MAP_IPOPT_ERRORS(MACRO)						\
  MACRO (Ipopt::Infeasible_Problem_Detected,				\
	 "Infeasible problem detected");				\
  MACRO (Ipopt::Search_Direction_Becomes_Too_Small,			\
	 "Search direction too small");					\
  MACRO (Ipopt::Diverging_Iterates, "Diverging iterates");		\
  MACRO (Ipopt::Maximum_Iterations_Exceeded,				\
	 "Maximum iterations exceeded");				\
  MACRO (Ipopt::Restoration_Failed, "Restoration failed");		\
  MACRO (Ipopt::Error_In_Step_Computation, "Error in step computation"); \
  MACRO (Ipopt::Not_Enough_Degrees_Of_Freedom,				\
	 "Not enough degrees of freedom");				\
  MACRO (Ipopt::Invalid_Problem_Definition,				\
	 "Invalid problem definition");					\
  MACRO (Ipopt::Invalid_Option, "Invalid option");			\
  MACRO (Ipopt::Invalid_Number_Detected, "Invalid number detected");	\
  MACRO (Ipopt::Unrecoverable_Exception, "Unrecoverable exception");	\
  MACRO (Ipopt::Insufficient_Memory, "Insufficient memory");		\
  MACRO (Ipopt::Internal_Error, "Internal error");			\
  MACRO (Ipopt::Maximum_CpuTime_Exceeded, "Maximum CPU time exceeded")

#define MAP_IPOPT_FATALS(MACRO)						\
  MACRO(Ipopt::NonIpopt_Exception_Thrown, "Non-Ipopt exception thrown")

#define MAP_IPOPT_OKS(MACRO)						\
  MACRO (Ipopt::Solve_Succeeded, MAP_IPOPT_ERRORS(SWITCH_ERROR);	\
	 MAP_IPOPT_FATALS(SWITCH_FATAL));				\
  MACRO (Ipopt::Solved_To_Acceptable_Level, MAP_IPOPT_ERRORS(SWITCH_ERROR); \
	 MAP_IPOPT_FATALS(SWITCH_FATAL));				\
  MACRO (Ipopt::Feasible_Point_Found, MAP_IPOPT_ERRORS(SWITCH_ERROR);	\
	 MAP_IPOPT_FATALS(SWITCH_FATAL));				\
  MACRO (Ipopt::User_Requested_Stop, MAP_IPOPT_ERRORS(SWITCH_ERROR);	\
	 MAP_IPOPT_FATALS(SWITCH_FATAL))

  template<typename T>
  void IpoptSolverCommon<T>::
  solve ()
  {
    // Read parameters and forward them to Ipopt.
    updateParameters ();
    Ipopt::ApplicationReturnStatus status = app_->Initialize ("");

    switch (status)
      {
	MAP_IPOPT_OKS (SWITCH_OK);
	MAP_IPOPT_ERRORS (SWITCH_ERROR);
	MAP_IPOPT_FATALS (SWITCH_FATAL);
      }
    assert (this->result_.which () != T::SOLVER_NO_SOLUTION);
  }

#undef SWITCH_ERROR
#undef SWITCH_FATAL
#undef SWITCH_OK
#undef MAP_IPOPT_ERRORS
#undef MAP_IPOPT_FATALS
#undef MAP_IPOPT_OKS

  template<typename T>
  Ipopt::SmartPtr<Ipopt::IpoptApplication> IpoptSolverCommon<T>::
  getIpoptApplication ()
  {
    return app_;
  }

#define DEFINE_PARAMETER(KEY, DESCRIPTION, VALUE)	\
  do {							\
    this->parameters_[KEY].description = DESCRIPTION;	\
    this->parameters_[KEY].value = VALUE;		\
  } while (0)

  template<typename T>
  void IpoptSolverCommon<T>::
  initializeParameters ()
  {
    this->parameters_.clear ();

    // Shared parameters.
    DEFINE_PARAMETER ("max-iterations", "number of iterations", 3000);

    // IPOPT specific.
    // Much more options are available for Ipopt, see ``Options reference'' of
    // Ipopt documentation.

    //  Output
    DEFINE_PARAMETER ("ipopt.print_level", "output verbosity level", 5);
    DEFINE_PARAMETER ("ipopt.print_user_options",
		      "print all options set by the user", std::string ("no"));
    DEFINE_PARAMETER ("ipopt.print_options_documentation",
		      "switch to print all algorithmic options", std::string ("no"));
    DEFINE_PARAMETER
      ("ipopt.output_file",
       "file name of desired output file (leave unset for no file output)",
       std::string (""));
    DEFINE_PARAMETER ("ipopt.file_print_level",
		      "verbosity level for output file", 5);
    DEFINE_PARAMETER ("ipopt.option_file_name",
		      "file name of options file (to overwrite default)", std::string (""));

    //  Termination
    DEFINE_PARAMETER ("ipopt.tol",
		      "desired convergence tolerance (relative)", 1e-8);
    DEFINE_PARAMETER ("ipopt.dual_inf_tol",
		      "desired threshold for the dual infeasibility", 1.);
    DEFINE_PARAMETER ("ipopt.constr_viol_tol",
		      "desired threshold for the constraint violation", 1e-4);
    DEFINE_PARAMETER
      ("ipopt.compl_inf_tol",
       "desired threshold for the complementarity conditions", 1e-4);

    DEFINE_PARAMETER ("ipopt.acceptable_tol",
		      "\"acceptable\" convergence tolerance (relative)", 1e-6);
    DEFINE_PARAMETER
      ("ipopt.acceptable_iter",
       "number of \"acceptable\" iterates before triggering termination", 15);
    DEFINE_PARAMETER
      ("ipopt.acceptable_constr_viol_tol",
       "\"acceptance\" threshold for the constraint violation", 1e-2);
    DEFINE_PARAMETER
      ("ipopt.acceptable_dual_inf_tol",
       "\"acceptance\" threshold for the dual infeasibility", 1e10);
    DEFINE_PARAMETER
      ("ipopt.acceptable_compl_inf_tol",
       "\"acceptance\" threshold for the complementarity conditions", 1e-2);

    // NLP scaling
    DEFINE_PARAMETER
      ("ipopt.nlp_scaling_method",
       "technique used for scaling the problem internally before it is solved",
       std::string ("gradient-based"));
    DEFINE_PARAMETER
      ("ipopt.nlp_scaling_max_gradient",
       "maximum gradient after NLP scaling", 1e2);

    //  Barrier parameter
    DEFINE_PARAMETER ("ipopt.mu_strategy",
		      "update strategy for barrier parameter", std::string ("adaptive"));

    //  Restoration Phase
    DEFINE_PARAMETER
      ("ipopt.expect_infeasible_problem",
       "enable heuristics to quickly detect an infeasible problem",
       std::string ("no"));
    DEFINE_PARAMETER
      ("ipopt.start_with_resto",
       "tells algorithm to switch to restoration phase in first iteration",
       std::string ("no"));

    // Linear solver choice.
    DEFINE_PARAMETER ("ipopt.linear_solver", "linear solver",
                      std::string (IPOPT_DEFAULT_LINEAR_SOLVER));

    // Derivative test.
    DEFINE_PARAMETER ("ipopt.derivative_test", "enable derivative checker",
                      std::string ("none"));
  }

#undef DEFINE_PARAMETER

  template<typename T>
  void IpoptSolverCommon<T>::
  updateParameters ()
  {
    const std::string prefix = "ipopt.";
    typedef const std::pair<const std::string, Parameter> const_iterator_t;
    BOOST_FOREACH (const_iterator_t& it, this->parameters_)
      {
	if (it.first.substr (0, prefix.size ()) == prefix)
	  {
	    boost::apply_visitor
	      (IpoptParametersUpdater
	       (app_, it.first.substr (prefix.size ())), it.second.value);
	  }
      }

    // Remap standardized parameters.
    boost::apply_visitor
      (IpoptParametersUpdater
       (app_, "max_iter"), this->parameters_["max-iterations"].value);
  }
} // end of namespace roboptim

#endif //! ROBOPTIM_CORE_IPOPT_COMMON_HXX
