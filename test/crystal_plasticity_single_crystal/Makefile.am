bin_PROGRAMS = test_crystal_plasticity_single_crystal

test_crystal_plasticity_single_crystal_SOURCES = test_crystal_plasticity_single_crystal.cc
test_crystal_plasticity_single_crystal_LDADD = $(top_builddir)/lib/libConstitutiveModel.la
test_crystal_plasticity_single_crystal_CPPFLAGS = \
  -I$(top_srcdir)/local_include \
  -I$(top_srcdir)/crystal_plasticity/include \
  -I$(top_srcdir)/elasticity/include \
  -I$(top_srcdir)/utils/include \
  -I$(top_srcdir)/material/include

TESTS = test_crystal_plasticity.sh

test_crystal_plasticity.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../run.sh test_crystal_plasticity_single_crystal crystal_plasticity_single_crystal single_crystal_results" >> $@
	chmod +x $@
