# ---- libraries ----------------------------------------------------------

add_library(fcm_core
        src/core/legendre.cpp
        src/core/gauss.cpp
        src/core/linalg.cpp
        src/fcm3d/shape.cpp
)
target_include_directories(fcm_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(fcm1d_lib
        src/fcm1d/config.cpp
        src/fcm1d/mesh.cpp
        src/fcm1d/partition.cpp
        src/fcm1d/quadrature.cpp
        src/fcm1d/assembly.cpp
        src/fcm1d/solver.cpp
)
target_link_libraries(fcm1d_lib PUBLIC fcm_core)

add_library(fcm3d_lib
        src/fcm3d/shape.cpp
        src/fcm3d/config.cpp
        src/fcm3d/quadrature.cpp
        src/fcm3d/assembly.cpp
        tests/fcm3d/assembly.cpp
        src/fcm3d/mesh.cpp
        tests/fcm3d/mesh.cpp
        src/fcm3d/partition.cpp
        tests/fcm3d/partition.cpp
        src/fcm3d/solver.cpp
        tests/fcm3d/solve.cpp
        apps/fcm3d.cpp
)
target_link_libraries(fcm3d_lib PUBLIC fcm_core)

# ---- executables ----------------------------------------------------------------

add_executable(fcm1d apps/fcm1d.cpp)
target_link_libraries(fcm1d PRIVATE fcm1d_lib)

add_executable(fcm3d apps/fcm3d.cpp)
target_link_libraries(fcm3d PRIVATE fcm3d_lib)

add_executable(fcm3d_sphere apps/fcm3d_sphere.cpp)
target_link_libraries(fcm3d_sphere PRIVATE fcm3d_lib)

# ---- tests ---------------------------------------------------------------

function(fcm_add_test group name lib)
    set(target test_${group}_${name})
    add_executable(${target} tests/${group}/${name}.cpp)
    target_link_libraries(${target} PRIVATE ${lib})
    target_compile_definitions(${target} PRIVATE
            FCM_REFERENCE_DIR="${CMAKE_SOURCE_DIR}/reference")
    add_test(NAME ${target} COMMAND ${target} "${CMAKE_SOURCE_DIR}/reference")
endfunction()

fcm_add_test(core  legendre   fcm_core)
fcm_add_test(core  gauss      fcm_core)

fcm_add_test(fcm1d mesh       fcm1d_lib)
fcm_add_test(fcm1d partition  fcm1d_lib)
fcm_add_test(fcm1d quadrature fcm1d_lib)
fcm_add_test(fcm1d assembly   fcm1d_lib)
fcm_add_test(fcm1d solve      fcm1d_lib)

fcm_add_test(fcm3d shape      fcm3d_lib)
fcm_add_test(fcm3d assembly fcm3d_lib)
fcm_add_test(fcm3d mesh fcm3d_lib)
fcm_add_test(fcm3d partition fcm3d_lib)
fcm_add_test(fcm3d solve fcm3d_lib)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    add_test(NAME validation_1d_cpp
            COMMAND ${Python3_EXECUTABLE} check_validation.py $<TARGET_FILE:fcm1d>
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/python/FCM1D)

    add_test(NAME validation_1d_python
            COMMAND ${Python3_EXECUTABLE} check_validation.py
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/python/FCM1D)
endif()
