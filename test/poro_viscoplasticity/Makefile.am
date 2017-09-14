bin_PROGRAMS = poro_viscoplasticity

poro_viscoplasticity_SOURCES = poro_viscoplasticity.cc F_of_t.cc
poro_viscoplasticity_LDADD = $(top_builddir)/lib/libConstitutiveModel.la
poro_viscoplasticity_CPPFLAGS = -I$(top_srcdir)/material/include \
                                -I$(top_srcdir)/utils/include \
                                -I$(top_srcdir)/elasticity/include \
                                -I$(top_srcdir)/crystal_plasticity/include \
                                -I$(top_srcdir)/local_include/include \
                                -I$(top_srcdir)/poro_viscoplasticity/include

