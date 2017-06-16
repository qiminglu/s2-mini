add_executable(benchmark benchmark.cc)

target_link_libraries(benchmark 
    synergia_bunch 
    synergia_collective
    synergia_simulation 
    synergia_libFF
    synergia_lattice
    synergia_utils 
    synergia_foundation
    )
