bin_PROGRAMS = test_crystal_plasticity_with_thermal

test_crystal_plasticity_with_thermal_SOURCES = test_crystal_plasticity_with_thermal.cc
test_crystal_plasticity_with_thermal_LDADD = $(top_builddir)/lib/libConstitutiveModel.la
test_crystal_plasticity_with_thermal_CPPFLAGS = \
  -I$(top_srcdir)/local_include \
  -I$(top_srcdir)/crystal_plasticity/include \
  -I$(top_srcdir)/elasticity/include \
  -I$(top_srcdir)/utils/include \
  -I$(top_srcdir)/material/include 
