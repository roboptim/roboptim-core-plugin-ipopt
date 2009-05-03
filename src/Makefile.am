include $(top_srcdir)/build-aux/init.mk

# Add available warnings flags.
AM_CXXFLAGS += $(WARNING_CXXFLAGS)

lib_LTLIBRARIES = roboptim-core-ipopt-plugin.la

roboptim_core_ipopt_plugin_la_SOURCES = doc.hh roboptim-core-ipopt-plugin.cc
roboptim_core_ipopt_plugin_la_LIBADD = -lipopt
roboptim_core_ipopt_plugin_la_LDFLAGS = -module -version-info 1:1:0
