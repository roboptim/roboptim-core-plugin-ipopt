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

# Which linear solver should Ipopt use (unless explicitly set)?
# Default: mumps (open-source linear solver)
SET(IPOPT_DEFAULT_LINEAR_SOLVER "mumps"
    CACHE STRING "Which linear solver should be used by Ipopt by default?")

# Define the directory where plug-ins will be installed.
SET(PLUGINDIR ${CMAKE_INSTALL_LIBDIR}/${ROBOPTIM_CORE_RELPLUGINDIR})
MESSAGE("-- Plug-in directory: ${PLUGINDIR}")

MACRO(IPOPT_PLUGIN NAME)
  ADD_LIBRARY(roboptim-core-plugin-${NAME} MODULE
    ${NAME}.cc tnlp.cc tnlp.hh tnlp.hxx doc.hh ${HEADERS}
    )
  SET_TARGET_PROPERTIES(roboptim-core-plugin-${NAME} PROPERTIES
    COMPILE_DEFINITIONS
      IPOPT_DEFAULT_LINEAR_SOLVER="${IPOPT_DEFAULT_LINEAR_SOLVER}"
    PREFIX ""
    SOVERSION 3 VERSION 3.2.0)
  INSTALL(TARGETS roboptim-core-plugin-${NAME}
    DESTINATION ${PLUGINDIR})
  PKG_CONFIG_USE_DEPENDENCY(roboptim-core-plugin-${NAME} ipopt)
  PKG_CONFIG_USE_COMPILE_DEPENDENCY(roboptim-core-plugin-${NAME} roboptim-core)

  # Make sure all symbols are defined.
  # See:
  # http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=592752
  IF(CMAKE_SYSTEM_NAME MATCHES "Linux"
    AND EXISTS "/etc/debian_version"
    AND IPOPT_PREFIX MATCHES "/usr")
  TARGET_LINK_LIBRARIES(roboptim-core-plugin-${NAME} dmumps_seq)
  ENDIF()
ENDMACRO()

# Ipopt plug-in (with only gradient computation).
IPOPT_PLUGIN(ipopt)
# Ipopt plug-in (with only gradient computation and sparse matrices).
IPOPT_PLUGIN(ipopt-sparse)
# Ipopt twice differentiable (td) plug-in (require hessian computation).
IPOPT_PLUGIN(ipopt-td)
