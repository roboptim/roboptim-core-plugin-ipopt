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

#include <cstring>
#include <stdexcept>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpTNLP.hpp>

#include <roboptim/core/plugin/ipopt/ipopt-td.hh>
#include <roboptim/core/plugin/ipopt/ipopt-sparse.hh>

#include <roboptim/core/debug.hh>
#include <roboptim/core/util.hh>

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
    {
      using namespace boost;

      n = static_cast<Index> (costFunction_->inputSize ());
      m = static_cast<Index> (constraintsOutputSize ());

      const function_t::vector_t& x = solver_.startingPoint ();

      // Clear just in case
      constraintJacobians_.clear ();

      // Compute number of nonzeros in constraint Jacobian, and store its value
      // for the first iteration of the solver.
      nnz_jac_g = 0;
      typedef differentiableConstraints_t::const_iterator citer_t;
      for (citer_t it = differentiableConstraints_.begin ();
	   it != differentiableConstraints_.end (); ++it)
	{
          constraintJacobians_.push_back ((*it)->jacobian (x));
          differentiableFunction_t::jacobian_t& jac = constraintJacobians_.back ();
          jac.makeCompressed ();
	  nnz_jac_g += jac.nonZeros ();
	}

      nnz_h_lag = 0; // unused
      index_style = TNLP::C_STYLE;
      return true;
    }

    template <>
    bool
    Tnlp<IpoptSolverSparse>::eval_jac_g(Index n, const Number* x, bool,
					Index ROBOPTIM_DEBUG_ONLY(m),
					Index nele_jac,
					Index* iRow, Index *jCol,
					Number* values)
    {
      using namespace boost;
      ROBOPTIM_DEBUG_ONLY
	(size_type n_ = static_cast<size_type> (n));
      assert (costFunction_->inputSize () == n_);
      assert (constraintsOutputSize () == m);

      if (!jacobianBuf_)
	{
	  jacobianBuf_ = differentiableFunction_t::jacobian_t
	    (static_cast<differentiableFunction_t::matrix_t::Index> (constraintsOutputSize ()),
	     costFunction_->inputSize ());
	  jacobianBuf_->reserve (nele_jac);
	}

      if (!values)
	{
	  // Emptying iRow/jCol arrays.
	  std::memset (iRow, 0, static_cast<std::size_t> (nele_jac) * sizeof (Index));
	  std::memset (jCol, 0, static_cast<std::size_t> (nele_jac) * sizeof (Index));

	  // First evaluate the constraints at the starting point to build the
	  // constraints jacobian.
	  int idx = 0;
	  typedef differentiableConstraints_t::const_iterator citer_t;
	  unsigned constraintId = 0;

	  typedef Eigen::Triplet<double> triplet_t;
	  std::vector<triplet_t> coefficients;

          assert (constraintJacobians_.size () ==
                  differentiableConstraints_.size ());

	  for (citer_t it = differentiableConstraints_.begin ();
	       it != differentiableConstraints_.end ();
	       ++it, ++constraintId)
	    {
              // Using the values already computed in get_nlp_info
	      const differentiableFunction_t::jacobian_t&
                tmp_jac = constraintJacobians_[constraintId];

	      for (int k = 0; k < tmp_jac.outerSize (); ++k)
		for (differentiableFunction_t::jacobian_t::InnerIterator
		       it (tmp_jac, k); it; ++it)
		  {
                    const int row = static_cast<int> (idx + it.row ());
                    const int col = static_cast<int> (it.col ());
		    coefficients.push_back
		      (triplet_t (row, col, it.value ()));
		  }
	      idx += (*it)->outputSize ();
	    }

	  jacobianBuf_->setFromTriplets
	    (coefficients.begin (), coefficients.end ());
	  jacobianBuf_->makeCompressed ();

	  // Then look for non-zero values.
	  idx = 0;

	  for (int k = 0; k < jacobianBuf_->outerSize (); ++k)
	    for (differentiableFunction_t::jacobian_t::InnerIterator it (*jacobianBuf_, k);
		 it; ++it)
	      {
		iRow[idx] = it.row (), jCol[idx] = it.col ();
		++idx;
	      }

	  return true;
	}

      const Eigen::Map<const function_t::argument_t> x_ (x, n);

      typedef differentiableConstraints_t::const_iterator citer_t;

      size_t constraintId = 0;
      for (citer_t it = differentiableConstraints_.begin ();
	   it != differentiableConstraints_.end (); ++it, constraintId++)
	{
	  typename differentiableFunction_t::matrix_t&
	    jac = constraintJacobians_[constraintId];
	  // Set the Jacobian to 0 while keeping its structure (thus do not
	  // call setZero on the full matrix).
	  if (jac.isCompressed ())
	    {
	      Eigen::Map<Eigen::VectorXd> jacValues (jac.valuePtr (), jac.nonZeros ());
	      jacValues.setZero ();
	    }
	  (*it)->jacobian (jac, x_);

	  IpoptCheckGradient
	    (*(*it), 0, x_,
             static_cast<int> (constraintId), solver_);
	}

      // Copy jacobian values from internal sparse matrices.
      int idx = 0;

      if (StorageOrder == Eigen::ColMajor)
	{
	  for (int k = 0; k < solver_.problem ().function ().inputSize (); ++k)
	    for (constraintJacobians_t::const_iterator
		   g  = constraintJacobians_.begin ();
		 g != constraintJacobians_.end (); ++g)
	      {
		for (differentiableFunction_t::jacobian_t::InnerIterator it (*g, k);
		     it; ++it)
		  {
		    assert (idx < nele_jac);
		    values[idx++] = it.value ();
		  }
	      }
	}
      else
	{
	  for (constraintJacobians_t::const_iterator
		 g  = constraintJacobians_.begin ();
	       g != constraintJacobians_.end (); ++g)
	    for (int k = 0; k < g->outerSize(); ++k)
	      for (differentiableFunction_t::jacobian_t::InnerIterator it (*g, k);
		   it; ++it)
		{
		  assert (idx < nele_jac);
		  values[idx++] = it.value ();
		}
	}

      return true;
    }
  } // end of namespace detail
} // end of namespace roboptim
