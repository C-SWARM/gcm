bin_PROGRAMS = crystal_plasticity_grains

crystal_plasticity_grains_SOURCES = crystal_plasticity_grains.cc
crystal_plasticity_grains_LDADD = $(top_builddir)/lib/libConstitutiveModel.la
crystal_plasticity_grains_CPPFLAGS = -I$(top_srcdir)/local_include \
                                     -I$(top_srcdir)/crystal_plasticity/include \
                                     -I$(top_srcdir)/elasticity/include \
                                     -I$(top_srcdir)/utils/include \
                                     -I$(top_srcdir)/constitutive_model_handle/include \
                                     -I$(top_srcdir)/material/include 
