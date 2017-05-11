
add_library(synergia_foundation four_momentum.cc reference_particle.cc)

#add_library(synergia_foundation four_momentum.cc reference_particle.cc
#    distribution.cc multi_diagnostics.cc
#    diagnostics_write_helper.cc)

#target_link_libraries(synergia_foundation ${GSL_CONFIG_LIBS}
#    synergia_parallel_utils synergia_hdf5_utils synergia_serialization
#    lsexpr
#    ${Boost_SERIALIZATION_LIBRARIES})
