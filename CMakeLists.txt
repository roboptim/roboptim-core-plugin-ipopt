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

# Requires at least CMake 2.6 to configure the package.
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/cpack.cmake)

SET(PROJECT_NAME roboptim-core-plugin-ipopt)
SET(PROJECT_DESCRIPTION "RobOptim core IPOPT plug-in")
SET(PROJECT_URL "http://github.com/laas/roboptim-core-plugin-ipopt")

SET(HEADERS
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/plugin/ipopt/ipopt-common.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/plugin/ipopt/ipopt-parameters-updater.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/plugin/ipopt/ipopt-td.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/plugin/ipopt/ipopt.hh
  ${CMAKE_SOURCE_DIR}/include/roboptim/core/plugin/ipopt/ipopt-sparse.hh
  )

SET(PKG_CONFIG_ADDITIONAL_VARIABLES plugindir ${PKG_CONFIG_ADDITIONAL_VARIABLES})

SETUP_PROJECT()

# Search for dependencies.
SEARCH_FOR_BOOST()

# Libtool dynamic loading
# This project does not use Libtool directly but still uses ltdl for
# plug-in loading.
SET(LIB_LTDL LIB_LTDL-NOTFOUND)
FIND_LIBRARY(LIB_LTDL libltdl${CMAKE_SHARED_LIBRARY_SUFFIX})
IF (NOT LIB_LTDL)
  MESSAGE(FATAL_ERROR
    "Failed to find Libtool ltdl library, check that Libtool ltdl is installed.")
ENDIF()
INCLUDE(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(ltdl.h LTDL_H_FOUND)
IF (NOT LTDL_H_FOUND)
  MESSAGE(FATAL_ERROR
    "Failed to find ltdl.h, check that Libtool ltdl is installed.")
ENDIF()

# Search for Ipopt.
ADD_REQUIRED_DEPENDENCY("roboptim-core >= 0.5")
ADD_REQUIRED_DEPENDENCY("ipopt")

# patch for ipopt: in the .pc file, the header path may be IpoptPath/include/coin
#  for clarity, we want to have IpoptPath/include/ instead
# This patch assumes that the initial header path has the following form :
#  "install-folder/include/coin"
SET(IPOPT_INCLUDE_DIRS_PATCHED ${IPOPT_PREFIX}/include)

#add the path to ipopt to the current include path
INCLUDE_DIRECTORIES(SYSTEM ${IPOPT_INCLUDE_DIRS_PATCHED})

#add this path to the .pc file
PKG_CONFIG_APPEND_CFLAGS(-I${IPOPT_INCLUDE_DIRS_PATCHED})

HEADER_INSTALL("${HEADERS}")

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()
