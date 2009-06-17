include $(top_srcdir)/build-aux/init.mk

# Add pthread CFLAGS.
AM_CPPFLAGS += $(PTHREAD_CFLAGS)

# Add RobOptim flags.
AM_CPPFLAGS += $(ROBOPTIMCORE_CFLAGS)

# Add available warnings flags.
AM_CXXFLAGS += $(WARNING_CXXFLAGS)

lib_LTLIBRARIES = roboptim-core-ipopt-plugin.la

roboptim_core_ipopt_plugin_la_SOURCES = doc.hh roboptim-core-ipopt-plugin.cc


# Link against roboptim-core.
roboptim_core_ipopt_plugin_la_LIBADD = $(ROBOPTIMCORE_LIBS)

# Link against Ipopt.
roboptim_core_ipopt_plugin_la_LIBADD += -lipopt

# Fix yet another Ipopt problem...
# Ipopt may use libraries like pthreads and lapack
# but does not link against them even if they are needed.
# The plug-in have to do it to avoid installing libraries with
# undefined symbols.
# The work around is to link against lapack here, hopefully it will
# be solved upstream and this useless dependency will be removed.
roboptim_core_ipopt_plugin_la_LIBADD +=			\
			$(PTHREAD_LIBS) 		\
			-llapack -lblas -lm -ldl	\
			-lgfortranbegin -lgfortran


roboptim_core_ipopt_plugin_la_LDFLAGS = -module -version-info 1:1:0
