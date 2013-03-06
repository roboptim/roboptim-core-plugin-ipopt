#ifndef ROBOPTIM_CORE_PLUGING_IPOPT_TNLP_HXX
# define ROBOPTIM_CORE_PLUGING_IPOPT_TNLP_HXX
# include <boost/mpl/assert.hpp>
# include <boost/mpl/at.hpp>
# include <boost/mpl/size.hpp>
# include <boost/mpl/vector.hpp>

# include <roboptim/core/plugin/ipopt-td.hh>

namespace roboptim
{
  using namespace Ipopt;

  namespace detail
  {
    void
    jacobianFromGradients
    (DerivableFunction::matrix_t& jac,
     const IpoptSolver::problem_t::constraints_t& c,
     const DerivableFunction::vector_t& x);

    template  <typename T>
    Function::size_type
    computeConstraintsOutputSize (const T& pb)
    {
      BOOST_MPL_ASSERT_RELATION
	( (boost::mpl::size<typename T::constraintsList_t>::value), ==, 2);

      // Non-linear function type is supposed to be the second
      // constraint type and the linear function type is the
      // first.
      typedef typename
	boost::mpl::at<typename T::constraintsList_t,
		       boost::mpl::int_<1> >::type
	nonLinearFunction_t;

      typedef typename
	boost::mpl::at<typename T::constraintsList_t,
		       boost::mpl::int_<0> >::type
	linearFunction_t;

      Function::size_type result = 0;
      typedef typename T::constraints_t::const_iterator citer_t;
      for (citer_t it = pb.constraints ().begin ();
	   it != pb.constraints ().end (); ++it)
	{
	  shared_ptr<DifferentiableFunction> g;
	  if (it->which () == IpoptSolver::LINEAR)
	    g = get<shared_ptr<linearFunction_t> > (*it);
	  else
	    g = get<shared_ptr<nonLinearFunction_t> > (*it);
	  result += g->outputSize ();
	}
      return result;
    }

    template <typename T>
    Tnlp<T>::Tnlp (const typename solver_t::problem_t&, solver_t& solver)
      throw ()
      : solver_ (solver),
	cost_ (),
	costGradient_ (),
	constraints_ (),
	jacobian_ ()
    {
      BOOST_MPL_ASSERT_RELATION
	( (boost::mpl::size<
	     typename solver_t::problem_t::constraintsList_t>::value), ==, 2);
    }

    template <typename T>
    Function::size_type
    Tnlp<T>::constraintsOutputSize ()
    {
      return computeConstraintsOutputSize (solver_.problem ());
    }

    template <typename T>
    bool
    Tnlp<T>::get_nlp_info (Index& n, Index& m, Index& nnz_jac_g,
			Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
      throw ()
    {
      n = static_cast<Index> (solver_.problem ().function ().inputSize ());
      m = static_cast<Index> (constraintsOutputSize ());
      nnz_jac_g = n * m; //FIXME: use a dense matrix for now.
      nnz_h_lag = n * (n + 1) / 2; //FIXME: use a dense matrix for now.
      index_style = TNLP::C_STYLE;
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::get_bounds_info (Index n, Number* x_l, Number* x_u,
			   Index m, Number* g_l, Number* g_u)
      throw ()
    {
      assert (solver_.problem ().function ().inputSize () - n == 0);
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
    Tnlp<T>::get_scaling_parameters (Number&,
				  bool& use_x_scaling, Index n,
				  Number* x_scaling,
				  bool& use_g_scaling, Index m,
				  Number* g_scaling)
      throw ()
    {
      std::size_t n_ = static_cast<std::size_t> (n);
      std::size_t m_ = static_cast<std::size_t> (m);

      assert (solver_.problem ().argumentScales ().size () == n_);

      use_x_scaling = true, use_g_scaling = true;
      std::copy (solver_.problem ().argumentScales ().begin (),
		 solver_.problem ().argumentScales ().end (),
		 x_scaling);


      for (std::size_t i = 0; i < m_; ++i)
	     for (std::size_t j = 0;
		  j < solver_.problem ().scalesVector ()[i].size (); ++j)
		    g_scaling[i] = solver_.problem ().scalesVector ()[i][j];
		  return true;
    }

    template <typename T>
    bool
    Tnlp<T>::get_variables_linearity (Index n, LinearityType* var_types) throw ()
    {
      assert (solver_.problem ().function ().inputSize () - n == 0);
      //FIXME: detect from problem.
      for (Index i = 0; i < n; ++i)
	var_types[i] = TNLP::NON_LINEAR;
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::get_function_linearity (Index m, LinearityType* const_types) throw ()
    {
      assert (constraintsOutputSize () - m == 0);

      typedef typename solver_t::problem_t::constraints_t::const_iterator
	citer_t;

      unsigned idx = 0;
      for (citer_t it = solver_.problem ().constraints ().begin ();
	   it != solver_.problem ().constraints ().end (); ++it)

	{
	  LinearityType type =
	    (it->which () == LINEAR) ? TNLP::LINEAR : TNLP::NON_LINEAR;
	  shared_ptr<DifferentiableFunction> g;
	  if (type == LINEAR)
	    g = get<shared_ptr<linearFunction_t> > (*it);
	  else
	    g = get<shared_ptr<nonLinearFunction_t> > (*it);

	  for (unsigned j = 0; j < g->outputSize (); ++j)
	    const_types[idx++] = type;
	}
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::get_starting_point (Index n, bool init_x, Number* x,
			      bool init_z, Number* z_L, Number* z_U,
			      Index m, bool init_lambda,
			      Number*)
	     throw ()
    {
      assert (solver_.problem ().function ().inputSize () - n == 0);
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

      Eigen::Map<Function::result_t> x_ (x, n);
      x_ = *solver_.problem ().startingPoint ();
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::get_warm_start_iterate (IteratesVector&) throw ()
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
    Tnlp<T>::eval_f (Index n, const Number* x, bool new_x, Number& obj_value)
      throw ()
    {
      assert (solver_.problem ().function ().inputSize () - n == 0);

      if (new_x || !cost_)
	{
	  cost_ = Function::vector_t (1);
	  Eigen::Map<const Function::argument_t> x_ (x, n);
	  solver_.problem ().function () (*cost_, x_);
	}
      obj_value = (*cost_)[0];
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::eval_grad_f (Index n, const Number* x, bool new_x, Number* grad_f)
        throw ()
    {
      assert (solver_.problem ().function ().inputSize () - n == 0);

      if (new_x || !costGradient_)
	{
	  if (!costGradient_)
	    costGradient_ = Function::vector_t
	      (solver_.problem ().function ().inputSize ());

	  Eigen::Map<const Function::vector_t> x_ (x, n);
	  solver_.problem ().function ().gradient (*costGradient_, x_, 0);
	}

      Eigen::Map<Function::vector_t> grad_f_ (grad_f, n);
      grad_f_ =  *costGradient_;
      return true;
    }

    template <typename T>
    bool
    Tnlp<T>::eval_g (Index n, const Number* x, bool new_x,
		     Index m, Number* g)
	     throw ()
    {
      using namespace boost;
      Function::size_type n_ = static_cast<Function::size_type> (n);
      std::size_t m_ = static_cast<std::size_t> (m);

      assert (solver_.problem ().function ().inputSize () == n_);
      assert (solver_.problem ().constraints ().size () == m_);

      if (new_x || !constraints_)
	{
	  if (!constraints_)
	    constraints_ = Function::vector_t (constraintsOutputSize ());

#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
	  Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

	  Eigen::Map<const Function::vector_t> x_ (x, n);

	  typedef typename solver_t::problem_t::constraints_t::const_iterator
	    citer_t;

	  Function::size_type idx = 0;
	  for (citer_t it = solver_.problem ().constraints ().begin ();
	       it != solver_.problem ().constraints ().end (); ++it)
	    {
	      shared_ptr<DifferentiableFunction> g;
	      if (it->which () == LINEAR)
		g = get<shared_ptr<linearFunction_t> > (*it);
	      else
		g = get<shared_ptr<nonLinearFunction_t> > (*it);

	      constraints_->segment (idx, g->outputSize ()) = (*g) (x_);
	      idx += g->outputSize ();
	    }
	}

      Eigen::Map<Function::matrix_t> g_ (g, m, 1);
      g_ =  *constraints_;
      return true;
    }


    template <typename T>
    bool
    Tnlp<T>::eval_jac_g(Index n, const Number* x, bool new_x,
			Index m, Index, Index* iRow,
			Index *jCol, Number* values)
	     throw ()
    {
      using namespace boost;
      Function::size_type n_ = static_cast<Function::size_type> (n);
      std::size_t m_ = static_cast<std::size_t> (m);
      assert (solver_.problem ().function ().inputSize () == n_);
      assert (solver_.problem ().constraints ().size () == m_);

      if (!values)
	{
	  //FIXME: always dense for now.
	  int idx = 0;
	  // Eigen matrix are by default in colunmn major
	  // so a table 0 1 2 3 is see as a matrix: 0 2 by Eigen
	  //                                        1 3
	  // so we must fill (iRow,jCol) as 0:(0,0), 1:(1,0), 2:(0,1), 3:(1,1)
	  for (int j = 0; j < n; ++j)
	    for (int i = 0; i < m; ++i)
	      {
		iRow[idx] = i, jCol[idx] = j;
		++idx;
	      }
	}
      else
	{
	  if (new_x || !jacobian_)
	    {
	      if (!jacobian_)
		jacobian_ = Function::matrix_t
		  (constraintsOutputSize (),
		   solver_.problem ().function ().inputSize ());

	      Eigen::Map<const Function::vector_t> x_ (x, n);

	      typedef typename
		solver_t::problem_t::constraints_t::const_iterator
		citer_t;

	      Function::size_type idx = 0;
	      for (citer_t it = solver_.problem ().constraints ().begin ();
		   it != solver_.problem ().constraints ().end (); ++it)
		{
		  shared_ptr<DifferentiableFunction> g;
		  if (it->which () == LINEAR)
		    g = get<shared_ptr<linearFunction_t> > (*it);
		  else
		    g = get<shared_ptr<nonLinearFunction_t> > (*it);

		  jacobian_->block
		    (idx, 0, g->outputSize (), n) = g->jacobian (x_);
		  idx += g->outputSize ();
		}
	    }
	  Eigen::Map<Function::matrix_t> values_ (values, m, n);
	  values_ =  *jacobian_;
	}

      return true;
    }

    /// Compute Ipopt hessian from several hessians.
    template <>
    void
    Tnlp<IpoptSolverTd>::compute_hessian
      (TwiceDifferentiableFunction::hessian_t& h,
       const typename solver_t::vector_t& x,
       Number obj_factor,
       const Number* lambda)
        throw ()
      {
        typedef typename
	  solver_t::problem_t::constraints_t::const_iterator citer_t;

        TwiceDifferentiableFunction::hessian_t fct_h =
          solver_.problem ().function ().hessian (x, 0);
        h = obj_factor * fct_h;

        int i = 0;
        for (citer_t it = solver_.problem ().constraints ().begin ();
             it != solver_.problem ().constraints ().end (); ++it)
	  {
	    shared_ptr<TwiceDifferentiableFunction> g;
	    if (it->which () == LINEAR)
	      g = get<shared_ptr<linearFunction_t> > (*it);
	    else
	      g = get<shared_ptr<nonLinearFunction_t> > (*it);
	    h += lambda[i++] * g->hessian (x, 0);
	  }
      }

    template <>
    bool
    Tnlp<IpoptSolverTd>::eval_h
      (Index n, const Number* x, bool,
       Number obj_factor, Index m, const Number* lambda,
       bool, Index nele_hess, Index* iRow,
       Index* jCol, Number* values)
      throw ()
    {
      Function::size_type n_ = static_cast<Function::size_type> (n);
      std::size_t m_ = static_cast<std::size_t> (m);

      assert (solver_.problem ().function ().inputSize () == n_);
      assert (solver_.problem ().constraints ().size () == m_);

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
	  typename solver_t::vector_t x_ (n);
	  array_to_vector (x_, x);

	  TwiceDifferentiableFunction::hessian_t h
	    (solver_.problem ().function ().inputSize (),
	     solver_.problem ().function ().inputSize ());
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
      throw ()
    {
      return false;
    }


#define SWITCH_ERROR(NAME, ERROR)		\
      case NAME:				\
      solver_.result_ = SolverError (ERROR);	\
      break

#define SWITCH_FATAL(NAME)			\
      case NAME:				\
      assert (0);				\
      break

#define MAP_IPOPT_ERRORS(MACRO)						\
      MACRO(MAXITER_EXCEEDED, "Max iteration exceeded");		\
      MACRO(STOP_AT_TINY_STEP,						\
	    "Algorithm proceeds with very little progress");		\
      MACRO(LOCAL_INFEASIBILITY,					\
	    "Algorithm converged to a point of local infeasibility");	\
      MACRO(DIVERGING_ITERATES, "Iterate diverges");			\
      MACRO(RESTORATION_FAILURE, "Restoration phase failed");		\
      MACRO(ERROR_IN_STEP_COMPUTATION,					\
	    "Unrecoverable error while Ipopt tried to compute"		\
	    " the search direction");					\
      MACRO(INVALID_NUMBER_DETECTED,					\
	    "Ipopt received an invalid number");			\
      MACRO(INTERNAL_ERROR, "Unknown internal error");			\
      MACRO(TOO_FEW_DEGREES_OF_FREEDOM, "Two few degrees of freedom");	\
      MACRO(INVALID_OPTION, "Invalid option");				\
      MACRO(OUT_OF_MEMORY, "Out of memory");				\
      MACRO(CPUTIME_EXCEEDED, "Cpu time exceeded")


#define MAP_IPOPT_FATALS(MACRO)						\
      MACRO(USER_REQUESTED_STOP)

#define FILL_RESULT()					\
      array_to_vector (res.x, x);			\
      res.constraints.resize (m);			\
      array_to_vector (res.constraints, g);		\
      res.lambda.resize (m);				\
      array_to_vector (res.lambda, lambda);		\
      res.value (0) = obj_value


    template <typename T>
    void
    Tnlp<T>::finalize_solution
	     (SolverReturn status,
	      Index n, const Number* x, const Number*,
	      const Number*, Index m, const Number* g,
	      const Number* lambda, Number obj_value,
	      const IpoptData*,
	      IpoptCalculatedQuantities*)
        throw ()
    {
      assert (solver_.problem ().function ().inputSize () - n == 0);
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
	case STOP_AT_ACCEPTABLE_POINT:
	  {
	    ResultWithWarnings res (n, 1);
	    FILL_RESULT ();
	    res.warnings.push_back (SolverWarning ("Acceptable point"));
	    solver_.result_ = res;
	    break;
	  }

	  MAP_IPOPT_ERRORS(SWITCH_ERROR);
	  MAP_IPOPT_FATALS(SWITCH_FATAL);

	case UNASSIGNED:
	  assert (0 && "should never happen");
	}
      assert (solver_.result_.which () != IpoptSolver::SOLVER_NO_SOLUTION);
    }

#undef FILL_RESULT
#undef SWITCH_ERROR
#undef SWITCH_FATAL
#undef MAP_IPOPT_ERRORS
#undef MAP_IPOPT_FATALS

    template <typename T>
    bool
    Tnlp<T>::intermediate_callback (AlgorithmMode mode,
				    Index iter, Number obj_value,
				    Number inf_pr, Number inf_du,
				    Number mu, Number d_norm,
				    Number regularization_size,
				    Number alpha_du, Number alpha_pr,
				    Index ls_trials,
				    const IpoptData* ip_data,
				    IpoptCalculatedQuantities* ip_cq)
      throw ()
    {
      if (solver_.userIntermediateCallback ())
	return (*solver_.userIntermediateCallback ())
	  (mode, iter, obj_value, inf_pr, inf_du, mu, d_norm,
	   regularization_size, alpha_du, alpha_pr, ls_trials,
	   ip_data, ip_cq);
      return true;
    }

    template <typename T>
    Index
    Tnlp<T>::get_number_of_nonlinear_variables () throw ()
    {
      //FIXME: implement this.
      return -1;
    }

    template <typename T>
    bool
    Tnlp<T>::get_list_of_nonlinear_variables (Index, Index*)
	     throw ()
    {
      //FIXME: implement this.
      return false;
    }
  } // end of namespace detail
} // end of namespace roboptim

#endif //! ROBOPTIM_CORE_PLUGING_IPOPT_TNLP_HXX
