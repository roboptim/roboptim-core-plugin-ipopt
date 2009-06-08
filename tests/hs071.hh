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


#ifndef OPTIMIZATION_TESTS_HS071_HH
# define OPTIMIZATION_TESTS_HS071_HH
# include <utility>
# include <roboptim/core/twice-derivable-function.hh>

using namespace roboptim;

struct F : public TwiceDerivableFunction
{
  F () : TwiceDerivableFunction (4, 1)
  {
  }

  void
  impl_compute (result_t& result, const argument_t& x) const throw ()
  {
    result.clear ();
    result (0) = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[3];
  }

  void
  impl_gradient (gradient_t& grad, const argument_t& x, size_type) const throw ()
  {
    grad.clear ();
    grad[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
    grad[1] = x[0] * x[3];
    grad[2] = x[0] * x[3] + 1;
    grad[3] = x[0] * (x[0] + x[1] + x[2]);
  }

  void
  impl_hessian (hessian_t& h, const argument_t& x, size_type) const throw ()
  {
    h.clear ();
    h (0, 0) = 2 * x[3];
    h (0, 1) = x[3];
    h (0, 2) = x[3];
    h (0, 3) = 2 * x[0] + x[1] + x[2];

    h (1, 0) = x[3];
    h (1, 1) = 0.;
    h (1, 2) = 0.;
    h (1, 3) = x[0];

    h (2, 0) = x[3];
    h (2, 1) = 0.;
    h (2, 2) = 0.;
    h (2, 3) = x[1];

    h (3, 0) = 2 * x[0] + x[1] + x[2];
    h (3, 1) = x[0];
    h (3, 2) = x[0];
    h (3, 3) = 0.;
  }
};

struct G0 : public TwiceDerivableFunction
{
  G0 ()
    : TwiceDerivableFunction (4, 1)
  {
  }

  void
  impl_compute (result_t& res, const argument_t& x) const throw ()
  {
    res.clear ();
    res (0) = x[0] * x[1] * x[2] * x[3];
  }

  void
  impl_gradient (gradient_t& grad, const argument_t& x, size_type) const throw ()
  {
    grad.clear ();
    grad[0] = x[1] * x[2] * x[3];
    grad[1] = x[0] * x[2] * x[3];
    grad[2] = x[0] * x[1] * x[3];
    grad[3] = x[0] * x[1] * x[2];
  }

  void
  impl_hessian (hessian_t& h, const argument_t& x, size_type) const throw ()
  {
    h.clear ();
    h (0, 0) = 0.;
    h (0, 1) = x[2] * x[3];
    h (0, 2) = x[1] * x[3];
    h (0, 3) = x[1] * x[2];

    h (1, 0) = x[2] * x[3];
    h (1, 1) = 0.;
    h (1, 2) = x[0] * x[3];
    h (1, 3) = x[0] * x[2];

    h (2, 0) = x[1] * x[3];
    h (2, 1) = x[0] * x[3];
    h (2, 2) = 0.;
    h (2, 3) = x[0] * x[1];

    h (3, 0) = x[1] * x[2];
    h (3, 1) = x[0] * x[2];
    h (3, 2) = x[0] * x[1];
    h (3, 3) = 0.;
  }
};

struct G1 : public TwiceDerivableFunction
{
  G1 ()
    : TwiceDerivableFunction (4, 1)
  {
  }

  void
  impl_compute (result_t& res, const argument_t& x) const throw ()
  {
    res.clear ();
    res (0) = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
  }

  void
  impl_gradient (gradient_t& grad, const argument_t& x, size_type) const throw ()
  {
    grad.clear ();
    grad[0] = 2 * x[0];
    grad[1] = 2 * x[1];
    grad[2] = 2 * x[2];
    grad[3] = 2 * x[3];
  }

  void
  impl_hessian (hessian_t& h, const argument_t& x, size_type) const throw ()
  {
    h.clear ();
    h (0, 0) = 2.;
    h (0, 1) = 0.;
    h (0, 2) = 0.;
    h (0, 3) = 0.;

    h (1, 0) = 0.;
    h (1, 1) = 2.;
    h (1, 2) = 0.;
    h (1, 3) = 0.;

    h (2, 0) = 0.;
    h (2, 1) = 0.;
    h (2, 2) = 2.;
    h (2, 3) = 0.;

    h (3, 0) = 0.;
    h (3, 1) = 0.;
    h (3, 2) = 0.;
    h (3, 3) = 2.;
  }
};


template <typename T>
void initialize_problem (T& pb, const G0& g0, const G1& g1)
{
  // Set bound for all variables.
  // 1. < x_i < 5. (x_i in [1.;5.])
  for (Function::size_type i = 0; i < pb.function ().inputSize (); ++i)
    pb.argumentBounds ()[i] = Function::makeInterval (1., 5.);

  // Add constraints.
  pb.addConstraint (&g0, Function::makeLowerInterval (25.));
  pb.addConstraint (&g1, Function::makeInterval (40., 40.));

  // Set the starting point.
  Function::vector_t start (pb.function ().inputSize ());
  start[0] = 1., start[1] = 5., start[2] = 5., start[3] = 1.;
  pb.startingPoint () = start;
}

#endif //! OPTIMIZATION_TESTS_HS071_HH
