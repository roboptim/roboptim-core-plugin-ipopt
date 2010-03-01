include $(top_srcdir)/build-aux/init.mk
include $(top_srcdir)/ipopt.mk

# Add pthread CFLAGS.
AM_CPPFLAGS += $(PTHREAD_CFLAGS)

# Add RobOptim flags.
AM_CPPFLAGS += $(ROBOPTIMCORE_CFLAGS)

# Add available warnings flags.
AM_CXXFLAGS += $(WARNING_CXXFLAGS)

# Define plug-in directory.
plugindir = $(libdir)/roboptim-core

plugin_LTLIBRARIES = roboptim-core-ipopt-plugin.la

roboptim_core_ipopt_plugin_la_SOURCES = doc.hh roboptim-core-ipopt-plugin.cc


# Link against roboptim-core.
roboptim_core_ipopt_plugin_la_LIBADD = $(ROBOPTIMCORE_LIBS)

# Link against Ipopt.
roboptim_core_ipopt_plugin_la_LIBADD += $(IPOPT_LIBS)

roboptim_core_ipopt_plugin_la_LDFLAGS = -module -version-info 1:1:0

# Make sure one does not compile a plug-in with test coverage support as
# it seems to trigger a segv at exit.
roboptim_core_ipopt_plugin_la_CXXFLAGS = \
	`echo $(CXXFLAGS) | sed -e 's/-ftest-coverage//' -e 's/-fprofile-arcs//'`
