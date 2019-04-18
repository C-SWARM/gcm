bin_PROGRAMS = test_damage_split_model run_intg_alg_sdamage

test_damage_split_model_SOURCES = test_damage_split_model.cc
run_intg_alg_sdamage_SOURCES = run_intg_alg_sdamage.cc F_of_t.cc

LDADD = $(top_builddir)/lib/libConstitutiveModel.la

AM_CPPFLAGS = \
  -I$(top_srcdir)/local_include \
  -I$(top_srcdir)/damage/include \
  -I$(top_srcdir)/elasticity/include \
  -I$(top_srcdir)/crystal_plasticity/include \
  -I$(top_srcdir)/utils/include \
  -I$(top_srcdir)/constitutive_model_handle/include \
  -I$(top_srcdir)/material/include

