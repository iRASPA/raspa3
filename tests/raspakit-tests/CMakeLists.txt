add_executable(unit_tests_raspakit
               cbmc_interactions.cpp
               energy_decomposition.cpp 
               gradients.cpp
               pbc.cpp
               spacegroup.cpp
               dudlambda.cpp
               ewald.cpp
               insertion_deletion.cpp
               mc.cpp
               pressure.cpp
               static_energy.cpp
               main.cpp)

target_link_libraries(unit_tests_raspakit PRIVATE libs::raspakit pybind11::module ${GTEST_LIBRARIES})

add_test(NAME RaspaKit COMMAND unit_tests_raspakit)