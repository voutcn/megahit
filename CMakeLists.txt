cmake_minimum_required(VERSION 2.8)
project(megahit)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 11)

option(USE_POPCNT "Use popcnt hardware instructions" ON)
option(USE_BMI2 "Use bmi2 hardware instructions" ON)

include_directories(src)
include_directories(src/sparsepp)
FILE(GLOB_RECURSE ASMBL_SOURCE "src/assembly/*.cpp")
FILE(GLOB_RECURSE LCASM_SOURCE "src/localasm/*.cpp")
FILE(GLOB_RECURSE IDBA_SOURCE "src/idba/*.cpp")
FILE(GLOB_RECURSE SDBG_SOURCE "src/sdbg/*.cpp")
FILE(GLOB_RECURSE CX1_SOURCE "src/cx1/*.cpp")
FILE(GLOB_RECURSE TOOLKIT_SOURCE "src/tools/*.cpp")

LIST(APPEND OTHER_SOURCE
        src/main.cpp
        src/main_assemble.cpp
        src/main_buildlib.cpp
        src/main_iterate.cpp
        src/main_local_assemble.cpp
        src/main_sdbg_build.cpp
        src/options_description.cpp
        src/sequence_manager.cpp
        )

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

#set(STATIC_LIBS -lz -ldl -lpthread -lgomp -static-libstdc++ -static-libgcc -static)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXXH_INLINE_ALL -D__STDC_FORMAT_MACROS")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprefetch-loop-arrays -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__XFILE__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")


if (USE_POPCNT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
endif (USE_POPCNT)

if (USE_BMI2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -DUSE_BMI2")
endif (USE_BMI2)

set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O1")


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS}")

add_executable(megahit_core ${OTHER_SOURCE} ${ASMBL_SOURCE} ${IDBA_SOURCE} ${SDBG_SOURCE} ${LCASM_SOURCE}
        ${CX1_SOURCE} ${TOOLKIT_SOURCE})
target_link_libraries(megahit_core ${ZLIB_LIBRARIES} ${OpenMP_CXX_LIBRARIES})

add_custom_target(
        megahit
        COMMAND cp ${CMAKE_SOURCE_DIR}/src/megahit .
)

set(TEST_DATA ${CMAKE_SOURCE_DIR}/test_data)

add_custom_target(
        simple_test
        COMMAND rm -rf megahit_out*
        COMMAND ./megahit --12 ${TEST_DATA}/r1.il.fa.gz,${TEST_DATA}/r2.il.fa.bz2,${TEST_DATA}/r3.il.fa -o megahit_out1 -t 2 --keep-tmp-files
        COMMAND rm -rf megahit_out*
        COMMAND ./megahit --12 ${TEST_DATA}/r1.il.fa.gz,${TEST_DATA}/r2.il.fa.bz2,${TEST_DATA}/r3.il.fa -o megahit_out2 -t 2 --keep-tmp-files --kmin-1pass
)

add_dependencies(megahit megahit_core)
add_dependencies(simple_test megahit)

install(TARGETS megahit_core DESTINATION bin)
install(PROGRAMS src/megahit DESTINATION bin)
install(DIRECTORY test_data DESTINATION share/${PROJECT_NAME})
