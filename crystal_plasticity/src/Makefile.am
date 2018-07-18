noinst_LTLIBRARIES = libcrystal_plasticity.la

libcrystal_plasticity_la_SOURCES = construct_linearization_parts.cc \
                                   flowlaw.cc \
                                   hardening.cc \
                                   solve_system.cc
libcrystal_plasticity_la_CPPFLAGS = -I$(top_srcdir)/local_include \
                                    -I$(top_srcdir)/elasticity/include \
                                    -I$(top_srcdir)/material/include \
                                    -I$(top_srcdir)/utils/include \
                                    -I$(top_srcdir)/constitutive_model_handle/include \
                                    -I$(srcdir)/../include
