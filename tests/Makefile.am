include $(top_srcdir)/build-aux/init.mk
include $(top_srcdir)/build-aux/autotest.mk
include $(top_srcdir)/ipopt.mk

# Add pthread CFLAGS.
AM_CPPFLAGS += $(PTHREAD_CFLAGS)

# Add RobOptim flags.
AM_CPPFLAGS += $(ROBOPTIMCORE_CFLAGS)

# Search for headers in shared tests.
AM_CPPFLAGS += -I$(srcdir)/shared-tests

# Add available warnings flags.
AM_CXXFLAGS += $(WARNING_CXXFLAGS)

# Package the following files.
EXTRA_DIST += package.m4.in

# Clean atconfig (do not distribute it).
CLEANFILES += atconfig

# Source which are shared by all tests.
COMMON_SOURCES = common.hh shared-tests/hs071.hh

# Define tests.

# Simple
check_PROGRAMS = simple
simple_SOURCES = simple.cc $(COMMON_SOURCES)  \
	$(top_builddir)/src/roboptim-core-ipopt-plugin.cc
simple_LDADD = $(ROBOPTIMCORE_LIBS) $(IPOPT_LIBS)

# Plug-in
check_PROGRAMS += plugin
plugin_SOURCES = plugin.cc $(COMMON_SOURCES)
plugin_LDADD = $(ROBOPTIMCORE_LIBS)


.PHONY: generate-reference

# Re-generate reference outputs.
generate-reference:
	@echo "Generating references...";	\
	for i in $(check_PROGRAMS); do		\
	 IN="$(srcdir)/$$i.stdin";		\
	 OUT="$(srcdir)/$$i.stdout";		\
	 if test -f "$$IN"; then 		\
	  "./$$i" < $$IN > $$OUT || true;	\
	 else					\
	  "./$$i" > $$OUT || true;		\
	 fi					\
	done


# ------------ #
# Test suite.  #
# ------------ #

# Distribute reference files.
EXTRA_DIST +=					\
	simple.stdout				\
	plugin.stdout
