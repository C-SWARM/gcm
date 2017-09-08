bin_PROGRAMS = test_slip_system

test_slip_system_SOURCES = test_slip_system.cc
test_slip_system_LDADD = $(top_builddir)/lib/libConstitutiveModel.la
test_slip_system_CPPFLAGS = \
  -I$(top_srcdir)/local_include \
  -I$(top_srcdir)/crystal_plasticity/include \
  -I$(top_srcdir)/elasticity/include \
  -I$(top_srcdir)/utils/include \
  -I$(top_srcdir)/material/include 
