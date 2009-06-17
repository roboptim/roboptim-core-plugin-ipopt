#								-*- Automake -*-

# Link against Ipopt.
IPOPT_LIBS = -lipopt

# Fix yet another Ipopt problem...
# Ipopt may use libraries like pthreads and lapack
# but does not link against them even if they are needed.
# The plug-in have to do it to avoid installing libraries with
# undefined symbols.
# The work around is to link against lapack here, hopefully it will
# be solved upstream and this useless dependency will be removed.
IPOPT_LIBS +=				\
	$(PTHREAD_LIBS) 		\
	-llapack -lblas -lm -ldl	\
	-lgfortranbegin -lgfortran
