bin_PROGRAMS = test_J2_plasticity

test_J2_plasticity_SOURCES = test_J2_plasticity.cc
test_J2_plasticity_LDADD = $(top_builddir)/lib/libConstitutiveModel.la
test_J2_plasticity_CPPFLAGS = -I$(top_srcdir)/local_include \
                              -I$(top_srcdir)/elasticity/include \
                              -I$(top_srcdir)/utils/include \
                              -I$(top_srcdir)/material/include \
                              -I$(top_srcdir)/J2_plasticity/include

TESTS = test_J2_plasticity.sh

test_J2_plasticity.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../run.sh test_J2_plasticity J2_plasticity stress" >> $@
	chmod +x $@
