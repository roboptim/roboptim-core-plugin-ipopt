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

MACRO(IPOPT_PLUGIN NAME)
  ADD_LIBRARY(roboptim-core-plugin-${NAME} MODULE
    ${NAME}.cc tnlp.cc tnlp.hh tnlp.hxx doc.hh ${HEADERS}
    )
  SET_TARGET_PROPERTIES(roboptim-core-plugin-${NAME} PROPERTIES
    PREFIX ""
    SOVERSION 2 VERSION 2.0.0)
  INSTALL(TARGETS roboptim-core-plugin-${NAME}
    DESTINATION ${ROBOPTIM_CORE_PLUGINDIR})
  PKG_CONFIG_USE_DEPENDENCY(roboptim-core-plugin-${NAME} ipopt)
  PKG_CONFIG_USE_DEPENDENCY(roboptim-core-plugin-${NAME} roboptim-core)
ENDMACRO()

# Ipopt plug-in (with only gradient computation).
IPOPT_PLUGIN(ipopt)
# Ipopt plug-in (with only gradient computation and sparse matrices).
IPOPT_PLUGIN(ipopt-sparse)
# Ipopt twice differentiable (td) plug-in (require hessian computation).
IPOPT_PLUGIN(ipopt-td)
