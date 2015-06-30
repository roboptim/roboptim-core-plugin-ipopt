// Copyright (C) 2013 by Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_CORE_IPOPT_TNLP_HH
# define ROBOPTIM_CORE_IPOPT_TNLP_HH
# include <boost/mpl/at.hpp>
# include <boost/optional.hpp>

# include <coin/IpSmartPtr.hpp>
# include <coin/IpIpoptApplication.hpp>

# include <roboptim/core/plugin/ipopt/ipopt.hh>
# include <roboptim/core/solver-state.hh>

#ifdef ROBOPTIM_CORE_IPOPT_PLUGIN_CHECK_GRADIENT
# include <boost/format.hpp>
# include <roboptim/core/finite-difference-gradient.hh>
#endif //!ROBOPTIM_CORE_IPOPT_PLUGIN_CHECK_GRADIENT

namespace roboptim
{
  namespace detail
  {
    typedef Ipopt::Index Index;
    typedef Ipopt::Number Number;

    /// \internal
    /// Ipopt non linear problem definition.
    template <typename T>
    class Tnlp : public Ipopt::TNLP
    {
    public:
      typedef T solver_t;
      typedef SolverState<typename solver_t::problem_t> solverState_t;

      Tnlp (const typename solver_t::problem_t& pb, solver_t& solver);

      Function::size_type constraintsOutputSize ();

      virtual bool
      get_nlp_info (Index& n, Index& m, Index& nnz_jac_g,
                    Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style);

      virtual bool
      get_bounds_info (Index n, Number* x_l, Number* x_u,
                       Index m, Number* g_l, Number* g_u);

      virtual bool
      get_scaling_parameters (Number&,
                              bool& use_x_scaling, Index n,
                              Number* x_scaling,
                              bool& use_g_scaling, Index m,
                              Number* g_scaling);

      virtual bool
      get_variables_linearity (Index n, LinearityType* var_types);

      virtual bool
      get_constraints_linearity (Index m, LinearityType* const_types);

      virtual bool
      get_starting_point (Index n, bool init_x, Number* x,
                          bool init_z, Number* z_L, Number* z_U,
                          Index m, bool init_lambda,
                          Number*);

      virtual bool
      get_warm_start_iterate (Ipopt::IteratesVector&);

      virtual bool
      eval_f (Index n, const Number* x, bool new_x, Number& obj_value);

      virtual bool
      eval_grad_f (Index n, const Number* x, bool new_x, Number* grad_f);

      virtual bool
      eval_g (Index n, const Number* x, bool new_x,
              Index m, Number* g);

      virtual bool
      eval_jac_g(Index n, const Number* x, bool new_x,
                 Index m, Index, Index* iRow,
                 Index *jCol, Number* values);

      virtual bool
      eval_h (Index n, const Number* x, bool,
              Number obj_factor, Index m, const Number* lambda,
              bool, Index nele_hess, Index* iRow,
              Index* jCol, Number* values);

      virtual void
      finalize_solution(Ipopt::SolverReturn status,
                        Index n, const Number* x, const Number*,
                        const Number*, Index m, const Number* g,
                        const Number* lambda, Number obj_value,
                        const Ipopt::IpoptData*,
                        Ipopt::IpoptCalculatedQuantities*);

      virtual bool
      intermediate_callback (Ipopt::AlgorithmMode,
                             Index, Number,
                             Number, Number,
                             Number, Number,
                             Number,
                             Number, Number,
                             Index,
                             const Ipopt::IpoptData*,
                             Ipopt::IpoptCalculatedQuantities*);

      virtual Index
      get_number_of_nonlinear_variables ();

      virtual bool
      get_list_of_nonlinear_variables (Index,
                                       Index*);

    protected:
      void compute_hessian (TwiceDifferentiableFunction::hessian_t& h,
			    const typename solver_t::vector_t& x,
			    Number obj_factor,
			    const Number* lambda);

    /// \brief Pointer to function logger (see log4cxx documentation).
    static log4cxx::LoggerPtr logger;

    private:
      /// \brief Non-linear function type.
      ///
      /// The assumption is done this is the second
      /// element of the constraints types vector.
      typedef typename
      boost::mpl::at<typename solver_t::problem_t::constraintsList_t,
		     boost::mpl::int_<1> >::type
      nonLinearFunction_t;

      /// \brief Linear function type.
      ///
      /// The assumption is done this is the first
      /// element of the constraints types vector.
      typedef typename
      boost::mpl::at<typename solver_t::problem_t::constraintsList_t,
		     boost::mpl::int_<0> >::type
      linearFunction_t;

      /// \brief Reference to RobOptim structure which instantiated
      /// this one.
      solver_t& solver_;

      /// \brief Current state of the solver (used by the callback function).
      solverState_t solverState_;

      /// \brief Traits of the cost function of this problem
      typedef typename solver_t::problem_t::function_t::traits_t traits_t;

      /// \brief Cost function buffer.
      boost::optional<typename GenericFunction<traits_t>::result_t> cost_;

      /// \brief Cost gradient buffer.
      boost::optional<typename GenericDifferentiableFunction<traits_t>::gradient_t> costGradient_;

      /// \brief Constraints buffer.
      boost::optional<typename GenericFunction<traits_t>::result_t> constraints_;

      /// \brief Constraints jacobian buffer.
      boost::optional<typename GenericDifferentiableFunction<traits_t>::matrix_t> jacobian_;
    };
  } // end of namespace detail.
} // end of namespace roboptim.

# include "tnlp.hxx"
#endif //! ROBOPTIM_CORE_IPOPT_TNLP_HH
