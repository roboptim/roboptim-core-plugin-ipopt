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

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>

#include <roboptim/core/plugin/ipopt-td.hh>
#include <roboptim/core/plugin/ipopt-sparse.hh>
#include <roboptim/core/debug.hh>

#include "tnlp.hh"

namespace roboptim
{
  using namespace Ipopt;

  namespace detail
  {
    template <>
    bool
    Tnlp<IpoptSolverSparse>::get_nlp_info (Index& n, Index& m, Index& nnz_jac_g,
					   Index& nnz_h_lag,
					   TNLP::IndexStyleEnum& index_style)
      throw ()
    {
      n = static_cast<Index> (solver_.problem ().function ().inputSize ());
      m = static_cast<Index> (constraintsOutputSize ());

      // compute number of non zeros elements in jacobian constraint.
      nnz_jac_g = 0;
      typedef typename solver_t::problem_t::constraints_t::const_iterator
	citer_t;
      for (citer_t it = solver_.problem ().constraints ().begin ();
	   it != solver_.problem ().constraints ().end (); ++it)
	{
	  // FIXME: should make sure we are in the bounds.
	  typename function_t::vector_t x (n);
	  x.setZero ();

	  shared_ptr<typename solver_t::commonConstraintFunction_t> g;
	  if (it->which () == LINEAR)
	    g = get<shared_ptr<linearFunction_t> > (*it);
	  else
	    g = get<shared_ptr<nonLinearFunction_t> > (*it);

	  nnz_jac_g += g->jacobian (x).nonZeros ();
	}


      nnz_h_lag = 0; // unused
      index_style = TNLP::C_STYLE;
      return true;
    }

    template <>
    bool
    Tnlp<IpoptSolverSparse>::eval_jac_g(Index n, const Number* x, bool,
					Index m,
					Index nele_jac,
					Index* iRow, Index *jCol,
					Number* values)
      throw ()
    {
      using namespace boost;
      function_t::size_type n_ = static_cast<function_t::size_type> (n);
      assert (solver_.problem ().function ().inputSize () == n_);
      assert (constraintsOutputSize () == m);

      if (!jacobian_)
	{
	  jacobian_ = function_t::matrix_t
	    (static_cast<function_t::matrix_t::Index> (constraintsOutputSize ()),
	     solver_.problem ().function ().inputSize ());
	  jacobian_->reserve (nele_jac);
	}

      if (!values)
	{
	  LOG4CXX_TRACE
	    (GenericSolver::logger, "Looking for non-zeros elements.");
	  LOG4CXX_TRACE (GenericSolver::logger, "nele_jac = " << nele_jac);

	  // Emptying iRow/jCol arrays.
	  memset (iRow, 0, nele_jac * sizeof (Index));
	  memset (jCol, 0, nele_jac * sizeof (Index));

	  // First evaluate the constraints in zero to build the
	  // constraints jacobian.
	  int idx = 0;
	  typedef typename solver_t::problem_t::constraints_t::const_iterator
	    citer_t;
	  unsigned constraintId = 0;

	  typedef Eigen::Triplet<double> triplet_t;
	  std::vector<triplet_t> coefficients;

	  for (citer_t it = solver_.problem ().constraints ().begin ();
	       it != solver_.problem ().constraints ().end ();
	       ++it, ++constraintId)
	    {
	      LOG4CXX_TRACE
		(GenericSolver::logger,
		 "Compute jacobian of constraint id = " << constraintId
		 << "to count for non-zeros elements");

	      function_t::vector_t x (n);
	      // Look for a place to evaluate the jacobian of the
	      // current constraint.
	      // If we do not have an initial guess...
	      if (!solver_.problem ().startingPoint ())
		for (unsigned i = 0; i < x.size (); ++i)
		  {
		    // if constraint is in an interval, evaluate at middle.
		    if (solver_.problem ().boundsVector ()[constraintId][i].first
			!= Function::infinity ()
			&&
			solver_.problem ().boundsVector ()[constraintId][i].second
			!= Function::infinity ())
		      x[i] =
			(solver_.problem ().boundsVector ()
			 [constraintId][i].second
			 - solver_.problem ().boundsVector ()
			 [constraintId][i].first) / 2.;
		    // otherwise use the non-infinite bound.
		    else if (solver_.problem ().boundsVector ()
			     [constraintId][i].first
			     != Function::infinity ())
		      x[i] = solver_.problem ().boundsVector ()
			[constraintId][i].first;
		    else
		      x[i] = solver_.problem ().boundsVector ()
			[constraintId][i].second;
		  }
	      else // other use initial guess.
		x = *(solver_.problem ().startingPoint ());

	      shared_ptr<typename solver_t::commonConstraintFunction_t> g;
	      if (it->which () == LINEAR)
		g = get<shared_ptr<linearFunction_t> > (*it);
	      else
		g = get<shared_ptr<nonLinearFunction_t> > (*it);

	      typename function_t::jacobian_t jacobian = g->jacobian (x);
	      for (int k = 0; k < jacobian.outerSize (); ++k)
		for (typename function_t::jacobian_t::InnerIterator
		       it (jacobian, k); it; ++it)
		  coefficients.push_back
		    (triplet_t (idx + it.row (), it.col (), it.value ()));
	      idx += g->outputSize ();
	    }

	  jacobian_->setFromTriplets
	    (coefficients.begin (), coefficients.end ());

	  LOG4CXX_TRACE
	    (GenericSolver::logger, "full problem jacobian...\n" << *jacobian_);

	  // Then look for non-zero values.
	  LOG4CXX_TRACE (GenericSolver::logger, "filling iRow and jCol...");
	  idx = 0;

	  for (int k = 0; k < jacobian_->outerSize (); ++k)
	    for (function_t::jacobian_t::InnerIterator it (*jacobian_, k);
		 it; ++it)
	      {
		iRow[idx] = it.row (), jCol[idx] = it.col ();
		LOG4CXX_TRACE
		  (GenericSolver::logger, "row: " << it.row ()
		   << " / col: " << it.col ()
		   << " / index: " << it.index ()
		   << " / value: " << it.value ()
		   << "\nidx: " << idx);
		++idx;
	      }

	  return true;
	}

      Eigen::Map<const function_t::vector_t> x_ (x, n);

      typedef typename
	solver_t::problem_t::constraints_t::const_iterator
	citer_t;

      int idx = 0;
      int constraintId = 0;
      for (citer_t it = solver_.problem ().constraints ().begin ();
	   it != solver_.problem ().constraints ().end (); ++it)
	{
	  shared_ptr<typename solver_t::commonConstraintFunction_t> g;
	  if (it->which () == LINEAR)
	    g = get<shared_ptr<linearFunction_t> > (*it);
	  else
	    g = get<shared_ptr<nonLinearFunction_t> > (*it);

	  jacobian_->middleRows
	    (idx, g->outputSize ()) = g->jacobian (x_);
	  idx += g->outputSize ();

	  IpoptCheckGradient
	    (*g, 0, x_,
	     constraintId++, solver_);
	}

      // Copy jacobian values from interal sparse matrix.
      idx = 0;
      for (int k = 0; k < jacobian_->outerSize (); ++k)
	for (function_t::jacobian_t::InnerIterator it (*jacobian_, k);
	     it; ++it)
	  {
	    assert (idx < nele_jac);
	    values[idx++] = it.value ();
	  }

      return true;
    }
  } // end of namespace detail
} // end of namespace roboptim
