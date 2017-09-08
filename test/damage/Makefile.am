bin_PROGRAMS = test_damage_model

test_damage_model_SOURCES = test_damage_model.cc
test_damage_model_LDADD = $(top_builddir)/lib/libConstitutiveModel.la
test_damage_model_CPPFLAGS = \
  -I$(top_srcdir)/local_include \
  -I$(top_srcdir)/damage/include \
  -I$(top_srcdir)/elasticity/include \
  -I$(top_srcdir)/utils/include \
  -I$(top_srcdir)/material/include 
