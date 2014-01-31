# Copyright 2010, Thomas Moulard, LAAS-CNRS
#
# This file is part of roboptim-core.
# roboptim-core is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# roboptim-core is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Lesser Public License for more details.
# You should have received a copy of the GNU Lesser General Public License
# along with roboptim-core.  If not, see <http://www.gnu.org/licenses/>.

# Prepare test suite
SET(PLUGIN_PATH "${CMAKE_BINARY_DIR}/src")
INCLUDE(shared-tests/tests.cmake)

# Schittkowski problems
SET(SOLVER_NAME "ipopt")
SET(FUNCTION_TYPE ::roboptim::EigenMatrixDense)
SET(PROGRAM_SUFFIX "")
SET(COST_FUNCTION_TYPE ::roboptim::GenericDifferentiableFunction)
SET(CONSTRAINT_TYPE_1 ::roboptim::GenericLinearFunction)
SET(CONSTRAINT_TYPE_2 ::roboptim::GenericDifferentiableFunction)
BUILD_SCHITTKOWSKI_PROBLEMS()

SET(SOLVER_NAME "ipopt-sparse")
SET(FUNCTION_TYPE ::roboptim::EigenMatrixSparse)
SET(PROGRAM_SUFFIX "-sparse")
BUILD_SCHITTKOWSKI_PROBLEMS()

# RobOptim test problems
SET(SOLVER_NAME "ipopt")
SET(FUNCTION_TYPE ::roboptim::EigenMatrixDense)
SET(PROGRAM_SUFFIX "")
SET(COST_FUNCTION_TYPE ::roboptim::GenericDifferentiableFunction)
BUILD_ROBOPTIM_PROBLEMS()

SET(SOLVER_NAME "ipopt-sparse")
SET(FUNCTION_TYPE ::roboptim::EigenMatrixSparse)
SET(PROGRAM_SUFFIX "-sparse")
SET(COST_FUNCTION_TYPE ::roboptim::GenericDifferentiableFunction)
BUILD_ROBOPTIM_PROBLEMS()
