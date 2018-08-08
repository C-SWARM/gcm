noinst_PROGRAMS = poro_viscoplasticity \
                  single_test \
                  compute_stress \
                  compute_pc

poro_viscoplasticity_SOURCES = poro_viscoplasticity.cc F_of_t.cc
single_test_SOURCES          = single_test.cc
compute_stress_SOURCES       = compute_stress.cc
compute_pc_SOURCES       = compute_pc.cc
 
LDADD = $(top_builddir)/lib/libConstitutiveModel.la

AM_CPPFLAGS = -I$(top_srcdir)/material/include \
              -I$(top_srcdir)/utils/include \
              -I$(top_srcdir)/elasticity/include \
              -I$(top_srcdir)/crystal_plasticity/include \
              -I$(top_srcdir)/local_include/include \
              -I$(top_srcdir)/constitutive_model_handle/include \
              -I$(top_srcdir)/poro_viscoplasticity/include

