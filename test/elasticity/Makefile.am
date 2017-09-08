bin_PROGRAMS = test_elasticity

test_elasticity_SOURCES = test_elasticity.cc
test_elasticity_LDADD = $(top_builddir)/lib/libConstitutiveModel.la
test_elasticity_CPPFLAGS = \
  -I$(top_srcdir)/local_include \
  -I$(top_srcdir)/elasticity/include \
  -I$(top_srcdir)/utils/include \
  -I$(top_srcdir)/material/include 

TESTS = test_elasticity.sh

test_elasticity.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../run.sh test_elasticity elasticity stress" >> $@
	chmod +x $@
