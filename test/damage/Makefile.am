bin_PROGRAMS = test_damage_model

test_damage_model_SOURCES = test_damage_model.cc F_of_t.cc
test_damage_model_LDADD = $(top_builddir)/lib/libConstitutiveModel.la
test_damage_model_CPPFLAGS = \
  -I$(top_srcdir)/local_include \
  -I$(top_srcdir)/damage/include \
  -I$(top_srcdir)/elasticity/include \
  -I$(top_srcdir)/crystal_plasticity/include \
  -I$(top_srcdir)/utils/include \
  -I$(top_srcdir)/constitutive_model_handle/include \
  -I$(top_srcdir)/material/include 

TESTS = test_damage_model.sh

test_damage_model.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../run.sh 'test_damage_model params_compression.in' damage stress.sigma" >> $@
	chmod +x $@
