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

#ifndef ROBOPTIM_CORE_PLUGING_IPOPT_TNLP_HXX
# define ROBOPTIM_CORE_PLUGING_IPOPT_TNLP_HXX

# include <boost/shared_ptr.hpp>

# include <boost/mpl/assert.hpp>
# include <boost/mpl/at.hpp>
# include <boost/mpl/size.hpp>
# include <boost/mpl/vector.hpp>

# include <roboptim/core/plugin/ipopt/ipopt-td.hh>
# include <roboptim/core/plugin/ipopt/ipopt-sparse.hh>
# include <roboptim/core/debug.hh>

# include <coin/IpIpoptCalculatedQuantities.hpp>
# include <coin/IpIpoptData.hpp>
# include <coin/IpOrigIpoptNLP.hpp>
# include <coin/IpTNLPAdapter.hpp>

namespace roboptim
{
  using namespace Ipopt;

  namespace detail
  {
    template <typename T>
    log4cxx::LoggerPtr Tnlp<T>::logger
    (log4cxx::Logger::getLogger ("roboptim.ipopt"));

    /// \internal
#ifdef ROBOPTIM_CORE_IPOPT_PLUGIN_CHECK_GRADIENT
    template <typename T, typename F>
    void IpoptCheckGradient (const F& function,
			     unsigned functionId,
			     const Eigen::Map<const Function::vector_t>& x,
			     int constraintId,
			     T& solver)
    {
      using boost::format;
      try
	{
	  checkGradientAndThrow (function, functionId, x, 1.);
	}
      catch (BadGradient& bg)
	{
	  solver.invalidateGradient ();
	  std::cerr
	    << ((functionId < 0)
		? "Invalid cost function gradient:"
		: (format ("Invalid constraint function gradient (id = %1%):")
		   % constraintId).str ())
	    << std::endl
	    << function.getName ()
	    << std::endl
	    << bg
	    << std::endl;
	}
    }
#else
    template <typename T, typename F>
    void IpoptCheckGradient (const F&,
			     unsigned,
			     const Eigen::Map<const Function::vector_t>&,
			     int,
			     T&)
    {}
#endif //!ROBOPTIM_CORE_IPOPT_PLUGIN_CHECK_GRADIENT

    template <typename T>
    void fillMultipliers (Function::vector_t& multipliers,
                          const Function::value_type* z_L,
                          const Function::value_type* z_U,
                          const Function::value_type* lambda,
                          Function::size_type n,
                          Function::size_type m,
                          const Solver<T>& solver)
    {
      const Eigen::Map<const Function::argument_t> map_zL (z_L, n);
      const Eigen::Map<const Function::argument_t> map_zU (z_U, n);
      const Eigen::Map<const Function::argument_t> map_lambda (lambda, m);

      multipliers.resize (n + m + 1);
      multipliers.setZero ();

      // First, argument bounds multipliers
      Function::size_type i = 0;
      typedef IpoptSolver::problem_t::intervals_t::const_iterator citer_t;
      for (citer_t it = solver.problem ().argumentBounds ().begin ();
           it != solver.problem ().argumentBounds ().end (); ++it)
      {
        // Active lower bound
        // Note: we expect λ < 0 if active lower bound constraint, and λ > 0 if
        // active upper bound constraint
        if ((*it).first != -Function::infinity () && map_zL[i] > map_zU[i])
          multipliers[i] = -map_zL[i];
        // Active upper bound
        else if ((*it).second != Function::infinity () && map_zU[i] >= map_zL[i])
          multipliers[i] = map_zU[i];

        i++;
      }

      // Then constraint multipliers
      multipliers.segment (n, m) = map_lambda;

      // Finally, for 1D objective function: 1
      multipliers[n+m] = 1.;
    }

    void
    jacobianFromGradients
    (DerivableFunction::matrix_t& jac,
     const IpoptSolver::problem_t::constraints_t& c,
     const DerivableFunction::vector_t& x);

    template  <typename T>
    typename Tnlp<T>::size_type computeConstraintsOutputSize (const T& solver)
    {
      using namespace boost;

      typename Tnlp<T>::size_type result = 0;
      typedef typename T::problem_t::constraints_t::const_iterator citer_t;
      for (citer_t it = solver.problem ().constraints ().begin ();
	   it != solver.problem ().constraints ().end (); ++it)
	{
	  result += (*it)->outputSize ();
	}
      return result;
    }

    template <typename T>
    Tnlp<T>::Tnlp (const typename solver_t::problem_t& pb, solver_t& solver)
      : solver_ (solver),
        solverState_ (pb),
	costFunction_ (&pb.function()),
	differentiableCostFunction_ (),
	twiceDifferentiableCostFunction_ (),
	constraints_ (),
	differentiableConstraints_ (),
	twiceDifferentiableConstraints_ (),
	costBuf_ (),
	costGradientBuf_ (),
	constraintsBuf_ (),
	jacobianBuf_ (),
	constraintJacobians_ ()
    {
      if (pb.function().template asType<differentiableFunction_t>())
      {
        if (pb.function().template asType<twiceDifferentiableFunction_t>())
        {
          twiceDifferentiableCostFunction_ = pb.function().template castInto<twiceDifferentiableFunction_t>();
        }
        differentiableCostFunction_ = pb.function().template castInto<differentiableFunction_t>();
      }
      typedef typename T::problem_t::constraints_t::const_iterator citer_t;
      for (citer_t it = pb.constraints ().begin ();
           it != pb.constraints ().end (); ++it)
      {
        if ((*(*it)).template asType<differentiableFunction_t>())
        {
          if ((*(*it)).template asType<twiceDifferentiableFunction_t>())
          {
            twiceDifferentiableConstraints_.push_back(boost::static_pointer_cast<twiceDifferentiableFunction_t>(*it));
          }
            differentiableConstraints_.push_back(boost::static_pointer_cast<differentiableFunction_t>(*it));
        }
        constraints_.push_back(*it);
      }
    }

    template <typename T>
    typename Tnlp<T>::size_type Tnlp<T>::constraintsOutputSize ()
    {
      return computeConstraintsOutputSize (solver_);
    }

    template <>
    bool
    Tnlp<IpoptSolverSparse>::get_nlp_info (Index& n, Index& m, Index& nnz_jac_g,
					   Index& nnz_h_lag,
					   TNLP::IndexStyleEnum& index_style);

    template <typename T>
    bool
    Tnlp<T>::get_nlp_info (Index& n, Index& m, Index& nnz_jac_g,
                           Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
    {
      n = static_cast<Index> ((*costFunction_).inputSize ());
      m = static_cast<Index> (constraintsOutputSize ());
      nnz_jac_g = n * m;
      nnz_h_lag = n * (n + 1) / 2;
      index_style = TNLP::C_STYLE;
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::get_bounds_info (Index ROBOPTIM_DEBUG_ONLY(n), Number* x_l,
                              Number* x_u, Index ROBOPTIM_DEBUG_ONLY(m), Number* g_l,
                              Number* g_u)
    {
      assert ((*costFunction_).inputSize () - n == 0);
      assert (constraintsOutputSize () - m == 0);

      typedef IpoptSolver::problem_t::intervals_t::const_iterator citer_t;
      for (citer_t it = solver_.problem ().argumentBounds ().begin ();
	   it != solver_.problem ().argumentBounds ().end (); ++it)
	*(x_l++) = (*it).first, *(x_u++) = (*it).second;

      typedef IpoptSolver::problem_t::intervalsVect_t::const_iterator
	citerVect_t;

      for (citerVect_t it = solver_.problem ().boundsVector ().begin ();
	   it != solver_.problem ().boundsVector ().end (); ++it)
	for (citer_t it2 = it->begin (); it2 != it->end (); ++it2)
	  *(g_l++) = it2->first, *(g_u++) = it2->second;
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::get_scaling_parameters (Number& obj_scaling,
                                     bool& use_x_scaling,
                                     Index ROBOPTIM_DEBUG_ONLY(n),
                                     Number* x_scaling,
                                     bool& use_g_scaling, Index m,
                                     Number* g_scaling)
    {
      typedef typename solver_t::problem_t::scalingVect_t scalingVect_t;
      ROBOPTIM_DEBUG_ONLY(std::size_t n_ = static_cast<std::size_t> (n));
      ROBOPTIM_DEBUG_ONLY(std::size_t m_ = static_cast<std::size_t> (m));

      // TODO: add support for objective scaling in roboptim-core
      obj_scaling = 1.;

      assert (solver_.problem ().argumentScaling ().size () == n_);

      use_x_scaling = true;
      std::copy (solver_.problem ().argumentScaling ().begin (),
		 solver_.problem ().argumentScaling ().end (),
		 x_scaling);

      std::size_t i = 0;
      use_g_scaling = (m > 0);
      for (typename scalingVect_t::const_iterator
	     c  = solver_.problem ().scalingVector ().begin ();
           c != solver_.problem ().scalingVector ().end (); ++c)
	{
	  for (std::size_t j = 0; j < c->size (); ++j)
	    g_scaling[i++] = (*c)[j];
	}

      assert (i == m_);

      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::get_variables_linearity (Index n, LinearityType* var_types)
    {
      assert ((*costFunction_).inputSize () - n == 0);
      //FIXME: detect from problem.
      for (Index i = 0; i < n; ++i)
	var_types[i] = TNLP::NON_LINEAR;
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::get_constraints_linearity (Index ROBOPTIM_DEBUG_ONLY(m),
                                     LinearityType* const_types)
    {
      using namespace boost;

      assert (constraintsOutputSize () - m == 0);

      typedef typename constraints_t::const_iterator citer_t;

      unsigned idx = 0;
      for (citer_t it = constraints_.begin();
	   it != constraints_.end (); ++it)

	{
	  LinearityType type =
	    ((*it)->template asType<GenericLinearFunction<traits_t> >()) ? TNLP::LINEAR : TNLP::NON_LINEAR;

	  for (size_type j = 0; j < (*it)->outputSize (); ++j)
	    const_types[idx++] = type;
	}
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::get_starting_point (Index n, bool init_x, Number* x, bool init_z,
                                 Number* z_L, Number* z_U, Index ROBOPTIM_DEBUG_ONLY(m),
                                 bool ROBOPTIM_DEBUG_ONLY(init_lambda), Number*)
    {
      assert ((*costFunction_).inputSize () - n == 0);
      assert (constraintsOutputSize () - m == 0);

      //FIXME: handle all modes.
      assert(init_lambda == false);

      // Set bound multipliers.
      if (init_z)
	{
	  //FIXME: for now, if required, scale is one.
	  //When do we need something else?
	  for (Index i = 0; i < n; ++i)
	    z_L[i] = 1., z_U[i] = 1.;
	}

      // Set the starting point.
      if (!solver_.problem ().startingPoint () && init_x)
	{
	  solver_.result_ =
	    SolverError ("Ipopt method needs a starting point.");
	  return false;
	}
      if (!solver_.problem ().startingPoint ())
	return true;

      Eigen::Map<Function::argument_t> x_ (x, n);
      x_ = *solver_.problem ().startingPoint ();
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::get_warm_start_iterate (IteratesVector&)
    {
      //FIXME: implement this.
      //IteratesVector is defined in src/Algorithm/IteratesVector.hpp
      //and is not distributed (not a distributed header).
      //Hence, seems difficult to manipulate this type.
      //Idea 1: offer the possibility either to retrive this data after
      //solving a problem.
      //Idea 2: creating this type manually from problem + rough guess
      //(or previous solution).
      return false;
    }

    template <typename T>
    bool
    Tnlp<T>::eval_f (Index n, const Number* x, bool, Number& obj_value)
    {
      assert ((*costFunction_).inputSize () - n == 0);

      costBuf_ = typename function_t::vector_t (1);
      Eigen::Map<const typename function_t::argument_t> x_ (x, n);
      (*costFunction_) (*costBuf_, x_);

      obj_value = (*costBuf_)[0];
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::eval_grad_f (Index n, const Number* x, bool, Number* grad_f)
    {
      assert ((*costFunction_).inputSize () - n == 0);

      if (!costGradientBuf_)
	costGradientBuf_ = typename differentiableFunction_t::gradient_t
	  ((*costFunction_).inputSize ());

      const Eigen::Map<const typename function_t::argument_t> x_ (x, n);
      differentiableCostFunction_->gradient (*costGradientBuf_, x_, 0);

      IpoptCheckGradient
        (differentiableCostFunction_, 0, x_, -1, solver_);

      Eigen::Map<typename function_t::vector_t> grad_f_ (grad_f, n);
      grad_f_ =  *costGradientBuf_;
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::eval_g (Index n, const Number* x, bool,
		     Index m, Number* g)
    {
      using namespace boost;
      ROBOPTIM_DEBUG_ONLY(size_type n_ = static_cast<size_type> (n));

      assert ((*costFunction_).inputSize () == n_);
      assert (constraintsOutputSize () == m);

      if (!constraintsBuf_)
	constraintsBuf_ =
	  typename function_t::result_t (constraintsOutputSize ());

      Eigen::Map<const typename function_t::argument_t> x_ (x, n);

      typedef typename constraints_t::const_iterator citer_t;

      size_type idx = 0;
      for (citer_t it = constraints_.begin ();
	   it != constraints_.end (); ++it)
	{
	  constraintsBuf_->segment (idx, (*it)->outputSize ()) = (*(*it)) (x_);
	  idx += (*it)->outputSize ();
        }

      Eigen::Map<typename function_t::result_t> g_ (g, m, 1);
      g_ =  *constraintsBuf_;
      return true;
    }


    template <>
    bool
    Tnlp<IpoptSolverSparse>::eval_jac_g(Index n, const Number* x, bool,
					Index m,
					Index nele_jac,
					Index* iRow, Index *jCol,
					Number* values);

    template <typename T>
    bool
    Tnlp<T>::eval_jac_g(Index n, const Number* x, bool,
			Index m, Index, Index* iRow,
			Index *jCol, Number* values)
    {
      using namespace boost;

      ROBOPTIM_DEBUG_ONLY(size_type n_ = static_cast<size_type> (n));
      assert ((*costFunction_).inputSize () == n_);
      assert (constraintsOutputSize () == m);

      if (!values)
	{
	  int idx = 0;
	  // If dense RobOptim jacobian matrices are column-major:
	  if (GenericFunctionTraits<traits_t>::
	      StorageOrder == Eigen::ColMajor)
	    {
	      for (int j = 0; j < n; ++j)
		for (int i = 0; i < m; ++i)
		  {
		    iRow[idx] = i, jCol[idx] = j;
		    ++idx;
		  }
	    }
	  else
	    {
	      for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
		  {
		    iRow[idx] = i, jCol[idx] = j;
		    ++idx;
		  }
	    }
	}
      else
	{
	  if (!jacobianBuf_)
	    {
	      jacobianBuf_ = typename differentiableFunction_t::matrix_t
		(constraintsOutputSize (),
		 (*costFunction_).inputSize ());
	      jacobianBuf_->setZero ();
	    }

	  Eigen::Map<const typename function_t::vector_t> x_ (x, n);

	  typedef typename differentiableConstraints_t::const_iterator citer_t;

	  size_type idx = 0;
	  int constraintId = 0;
	  for (citer_t it = differentiableConstraints_.begin ();
	       it != differentiableConstraints_.end (); ++it)
	    {
	      (*it)->jacobian (jacobianBuf_->block (idx, 0, (*it)->outputSize (), n), x_);
	      idx += (*it)->outputSize ();

	      IpoptCheckGradient
		(*it, 0, x_,
		 constraintId++, solver_);
	    }

	  Eigen::Map<typename differentiableFunction_t::jacobian_t> values_ (values, m, n);
	  values_ =  *jacobianBuf_;
	}

      return true;
    }

    /// Compute Ipopt hessian from several hessians.
    template <>
    inline void
    Tnlp<IpoptSolverTd>::compute_hessian
    (TwiceDifferentiableFunction::hessian_t& h,
     const solver_t::vector_t& x,
     Number obj_factor,
     const Number* lambda)
    {
      using namespace boost;

      typedef twiceDifferentiableConstraints_t::const_iterator citer_t;

      TwiceDifferentiableFunction::hessian_t fct_h =
        twiceDifferentiableCostFunction_->hessian (x, 0);
      h = obj_factor * fct_h;

      int i = 0;
      for (citer_t it = twiceDifferentiableConstraints_.begin();
           it != twiceDifferentiableConstraints_.end (); ++it)
        {
          h += lambda[i++] * (*it)->hessian (x, 0);
        }
    }

    template <>
    inline bool
    Tnlp<IpoptSolverTd>::eval_h
    (Index n, const Number* x, bool,
     Number obj_factor, Index ROBOPTIM_DEBUG_ONLY(m), const Number* lambda,
     bool, Index ROBOPTIM_DEBUG_ONLY(nele_hess), Index* iRow,
     Index* jCol, Number* values)
    {
      ROBOPTIM_DEBUG_ONLY(size_type n_ = static_cast<size_type> (n));

      assert ((*costFunction_).inputSize () == n_);
      assert (constraintsOutputSize () == m);

      //FIXME: check if a hessian is provided.

      if (!values)
	{
	  //FIXME: always dense for now.
	  int idx = 0;
	  for (int i = 0; i < n; ++i)
	    for (int j = 0; j < i + 1; ++j)
	      {
		iRow[idx] = i, jCol[idx] = j;
		++idx;
	      }
	  assert (idx == nele_hess);
	}
      else
	{
	  Eigen::Map<const function_t::argument_t> x_ (x, n);

	  TwiceDifferentiableFunction::hessian_t h
	    ((*costFunction_).inputSize (),
	     (*costFunction_).inputSize ());
	  compute_hessian (h, x_, obj_factor, lambda);

	  int idx = 0;
	  for (int i = 0; i < n; ++i)
	    for (int j = 0; j < i + 1; ++j)
	      values[idx++] = h (i, j);
	  assert (idx == nele_hess);
	}

      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::eval_h
    (Index, const Number*, bool,
     Number, Index, const Number*,
     bool, Index, Index*,
     Index*, Number*)
    {
      return false;
    }

#define FILL_RESULT()							\
    res.x = Eigen::Map<const Function::argument_t> (x, n);		\
    res.constraints = Eigen::Map<const Function::vector_t> (g, m);	\
    res.lambda = Eigen::Map<const Function::vector_t> (lambda, m);	\
    fillMultipliers (res.lambda, z_L, z_U, lambda, n, m, solver_);	\
    res.value (0) = obj_value

#define SWITCH_ERROR(NAME, ERROR)			\
    case NAME:                                          \
    {                                                   \
      Result res (n, 1);				\
      FILL_RESULT ();					\
      solver_.result_ = SolverError (ERROR, res);	\
    }                                                   \
    break

#define SWITCH_WARNING(NAME, WARNING)			\
    case NAME:                                          \
    {                                                   \
      ResultWithWarnings res (n, 1);			\
      FILL_RESULT ();					\
      res.warnings.push_back (SolverWarning (WARNING)); \
      solver_.result_ = res;				\
    }                                                   \
    break


#define MAP_IPOPT_ERRORS(MACRO)						\
    MACRO(MAXITER_EXCEEDED, "Max iteration exceeded");                  \
    MACRO(STOP_AT_TINY_STEP,						\
          "Algorithm proceeds with very little progress");		\
    MACRO(LOCAL_INFEASIBILITY,                                          \
          "Algorithm converged to a point of local infeasibility");	\
    MACRO(DIVERGING_ITERATES, "Iterate diverges");			\
    MACRO(RESTORATION_FAILURE, "Restoration phase failed");		\
    MACRO(ERROR_IN_STEP_COMPUTATION,					\
          "Unrecoverable error while Ipopt tried to compute"		\
          " the search direction");					\
    MACRO(INVALID_NUMBER_DETECTED,					\
          "Ipopt received an invalid number");                          \
    MACRO(INTERNAL_ERROR, "Unknown internal error");			\
    MACRO(TOO_FEW_DEGREES_OF_FREEDOM, "Two few degrees of freedom");	\
    MACRO(INVALID_OPTION, "Invalid option");				\
    MACRO(OUT_OF_MEMORY, "Out of memory");				\
    MACRO(CPUTIME_EXCEEDED, "CPU time exceeded")


#define MAP_IPOPT_WARNINGS(MACRO)			\
    MACRO(USER_REQUESTED_STOP, "User-requested stop");	\
    MACRO(STOP_AT_ACCEPTABLE_POINT, "Acceptable point")


    template <typename T>
    void
    Tnlp<T>::finalize_solution
    (SolverReturn status,
     Index n, const Number* x, const Number* z_L,
     const Number* z_U, Index m, const Number* g,
     const Number* lambda, Number obj_value,
     const IpoptData*,
     IpoptCalculatedQuantities*)
    {
      assert ((*costFunction_).inputSize () - n == 0);
      assert (constraintsOutputSize () - m == 0);

      switch (status)
	{
	case FEASIBLE_POINT_FOUND:
	case SUCCESS:
	  {
	    Result res (n, 1);
	    FILL_RESULT ();
	    solver_.result_ = res;
	    break;
	  }

	  MAP_IPOPT_ERRORS(SWITCH_ERROR);
	  MAP_IPOPT_WARNINGS(SWITCH_WARNING);

	case UNASSIGNED:
	  assert (0 && "should never happen");
	}
      assert (solver_.result_.which () != IpoptSolver::SOLVER_NO_SOLUTION);
    }

#undef FILL_RESULT
#undef SWITCH_ERROR
#undef SWITCH_WARNING
#undef MAP_IPOPT_ERRORS
#undef MAP_IPOPT_WARNINGS

    template <typename T>
    bool
    Tnlp<T>::intermediate_callback (AlgorithmMode mode,
                                    Index /*iter*/, Number obj_value,
                                    Number /*inf_pr*/, Number /*inf_du*/,
				    Number /*mu*/, Number /*d_norm*/,
				    Number /*regularization_size*/,
				    Number /*alpha_du*/, Number /*alpha_pr*/,
				    Index /*ls_trials*/,
				    const IpoptData* ip_data,
				    IpoptCalculatedQuantities* ip_cq)
    {
      if (!solver_.callback ())
	return true;
      if (!ip_cq)
	return true;
      Ipopt::OrigIpoptNLP* orignlp = dynamic_cast<OrigIpoptNLP*>
	(GetRawPtr (ip_cq->GetIpoptNLP ()));
      if (!orignlp)
	return true;
      Ipopt::TNLPAdapter* tnlp_adapter = dynamic_cast<TNLPAdapter*>
	(GetRawPtr (orignlp->nlp ()));

      // current scaled optimization parameters
      tnlp_adapter->ResortX (*ip_data->curr ()->x (), solverState_.x ().data ());

      // unscale x
      const Eigen::Map<const typename T::argument_t>
        scaling (solver_.problem ().argumentScaling ().data (),
                 solver_.problem ().function ().inputSize ());
      solverState_.x ().array () /= scaling.array ();

      // unscaled objective value at the current point
      solverState_.cost () = obj_value;

      // unscaled constraint violation at the current point
      solverState_.constraintViolation () = ip_cq->unscaled_curr_nlp_constraint_violation (Ipopt::NORM_MAX);

      // handle extra relevant parameters
      solverState_.parameters()["ipopt.mode"].value =
	(mode == RegularMode)? std::string ("RegularMode") : std::string ("RestorationPhaseMode");
      solverState_.parameters()["ipopt.mode"].description
        = "Indicates the mode in which the algorithm is";

      bool stop_optim = false;
      solverState_.parameters ()["ipopt.stop"].value = stop_optim;
      solverState_.parameters ()["ipopt.stop"].description
        = "Whether to stop the optimization process";

      // call user-defined callback
      solver_.callback () (solver_.problem (), solverState_);

      // ipopt.stop may have been unintentionally removed in the callback
      try
        {
          stop_optim = solverState_.template getParameter<bool> ("ipopt.stop");
        }
      catch (std::out_of_range&)
        {
          stop_optim = false;
        }

      // return true if the solver should continue, false otherwise
      return !stop_optim;
    }

    template <typename T>
    Index
    Tnlp<T>::get_number_of_nonlinear_variables ()
    {
      //FIXME: implement this.
      return -1;
    }

    template <typename T>
    bool
    Tnlp<T>::get_list_of_nonlinear_variables (Index, Index*)
    {
      //FIXME: implement this.
      return false;
    }
  } // end of namespace detail
} // end of namespace roboptim

#endif //! ROBOPTIM_CORE_PLUGING_IPOPT_TNLP_HXX
