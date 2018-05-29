noinst_LTLIBRARIES = libporo_viscoplasticity.la

libporo_viscoplasticity_la_SOURCES = FeFpModels.cc  \
                                     FiniteStrainModels.cc  \
                                     KMS-IJSS2017.cc \
                                     MechanicalModels.cc  \
                                     TimeIntegration_Manager.cc \
                                     pvp_interface.cc \
                                     pvp_interface_2.cc \
                                     poro_visco_plasticity.cc

libporo_viscoplasticity_la_CPPFLAGS = -I$(top_srcdir)/material/include \
                                      -I$(top_srcdir)/local_include \
                                      -I$(top_srcdir)/utils/include \
                                      -I$(srcdir)/../include

