cmake_minimum_required(VERSION 2.8)
project(megahit)
set(CMAKE_VERBOSE_MAKEFILE ON)

if (CMAKE_VERSION VERSION_LESS "3.1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else ()
    set(CMAKE_CXX_STANDARD 11)
endif ()

option(COVERAGE "Generate coverage report" OFF)
option(STATIC_BUILD "Build static executation" OFF)
option(SANITIZER "Enable sanitizers" OFF)
option(TSAN "Enable thread sanitizers" OFF)

include_directories(src)
FILE(GLOB_RECURSE ASMBL_SOURCE "src/assembly/*.cpp")
FILE(GLOB_RECURSE LCASM_SOURCE "src/localasm/*.cpp")
FILE(GLOB_RECURSE IDBA_SOURCE "src/idba/*.cpp")
FILE(GLOB_RECURSE SDBG_SOURCE "src/sdbg/*.cpp")
FILE(GLOB_RECURSE CX1_SOURCE "src/sorting/*.cpp")
FILE(GLOB_RECURSE SEQ_SOURCE "src/sequence/*.cpp")
FILE(GLOB_RECURSE TOOLKIT_SOURCE "src/tools/*.cpp")

LIST(APPEND OTHER_SOURCE
        src/main.cpp
        src/main_assemble.cpp
        src/main_buildlib.cpp
        src/main_iterate.cpp
        src/main_local_assemble.cpp
        src/main_sdbg_build.cpp
        src/utils/options_description.cpp
        )

if (STATIC_BUILD)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif (STATIC_BUILD)

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXXH_INLINE_ALL -ftemplate-depth=3000")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprefetch-loop-arrays -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__XROOT__='\"${CMAKE_SOURCE_DIR}/src\"'")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__XFILE__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ZLIB_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O1 -D_LIBCPP_DEBUG -D_GLIBCXX_DEBUG")

if (COVERAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
    set(COV_PY "coverage run")
endif (COVERAGE)

if (SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=leak -fsanitize=undefined")
endif (SANITIZER)

if (TSAN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
endif (TSAN)


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS}")

add_executable(megahit_core ${OTHER_SOURCE} ${ASMBL_SOURCE} ${IDBA_SOURCE} ${SDBG_SOURCE} ${LCASM_SOURCE} ${SEQ_SOURCE}
        ${CX1_SOURCE} ${TOOLKIT_SOURCE})
add_executable(megahit_core_popcnt ${OTHER_SOURCE} ${ASMBL_SOURCE} ${IDBA_SOURCE} ${SDBG_SOURCE} ${LCASM_SOURCE} ${SEQ_SOURCE}
        ${CX1_SOURCE} ${TOOLKIT_SOURCE})
add_executable(megahit_core_no_hw_accel ${OTHER_SOURCE} ${ASMBL_SOURCE} ${IDBA_SOURCE} ${SDBG_SOURCE} ${LCASM_SOURCE}
        ${SEQ_SOURCE} ${CX1_SOURCE} ${TOOLKIT_SOURCE})
set_target_properties(megahit_core PROPERTIES COMPILE_FLAGS "-mbmi2 -DUSE_BMI2 -mpopcnt")
set_target_properties(megahit_core_popcnt PROPERTIES COMPILE_FLAGS "-mpopcnt")

if (STATIC_BUILD)
    # TODO dirty
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -static")
    set_target_properties(megahit_core megahit_core_popcnt megahit_core_no_hw_accel PROPERTIES LINK_SEARCH_START_STATIC ON)
    set_target_properties(megahit_core megahit_core_popcnt megahit_core_no_hw_accel PROPERTIES LINK_SEARCH_END_STATIC ON)
endif (STATIC_BUILD)

target_link_libraries(megahit_core ${ZLIB_LIBRARIES})
target_link_libraries(megahit_core_popcnt ${ZLIB_LIBRARIES})
target_link_libraries(megahit_core_no_hw_accel ${ZLIB_LIBRARIES})

add_custom_target(
        megahit ALL
        COMMAND cp ${CMAKE_SOURCE_DIR}/src/megahit .
)

add_custom_target(megahit_toolkit ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink megahit_core_no_hw_accel megahit_toolkit)

set(TEST_DATA ${CMAKE_SOURCE_DIR}/test_data)

add_custom_target(
        simple_test
        COMMAND ./megahit --test -t 2
        COMMAND MEGAHIT_NUM_MERCY_FACTOR=1.5 ./megahit --test -t 4 --mem-flag 0 --no-hw-accel
        COMMAND ./megahit --test -t 2 --kmin-1pass --prune-level 3 --prune-depth 0
        COMMAND rm -rf test-random && python3 ${TEST_DATA}/generate_random_fasta.py > random.fa && ./megahit -r random.fa --k-list 255 --min-count 1 -o test-random
        COMMAND rm -rf test-fastg && ./megahit --test -t 2 --mem-flag 2 --keep-tmp-files -o test-fastg
        COMMAND rm -rf test-empty && ./megahit -r ${TEST_DATA}/empty.fa -o test-empty
        COMMAND rm -rf test-no-contig && ./megahit -r ${TEST_DATA}/r4.fa -o test-no-contig
        COMMAND ./megahit_toolkit contig2fastg 59 test-fastg/intermediate_contigs/k59.contigs.fa > 59.fastg
        COMMAND ./megahit_toolkit readstat < test-fastg/intermediate_contigs/k59.contigs.fa
)

add_dependencies(megahit megahit_core megahit_core_popcnt megahit_core_no_hw_accel megahit_toolkit)
add_dependencies(simple_test megahit)

install(TARGETS megahit_core megahit_core_popcnt megahit_core_no_hw_accel DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/megahit ${CMAKE_CURRENT_BINARY_DIR}/megahit_toolkit DESTINATION bin)
install(DIRECTORY test_data DESTINATION share/${PROJECT_NAME})
