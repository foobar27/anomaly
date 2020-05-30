cmake_minimum_required(VERSION 3.17)

project(anomaly)

find_package(Boost 1.72 REQUIRED COMPONENTS
             date_time filesystem iostreams)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    add_compile_options(-Wall -Wextra -Wwrite-strings -Wenum-compare
                        -Werror=return-type -Wparentheses -Wmissing-include-dirs
                        -Wnon-virtual-dtor
                        -Woverloaded-virtual
                        -pedantic
                        -Wenum-compare
                        -Wlogical-op
                        -Wduplicated-cond
                        -fstack-protector-strong
                        -Wno-missing-field-initializers
                    )

elseif(CMAKE_CXX_COMPILER_ID MATCHES "^.*Clang$")
    add_compile_options(
        -Wall
        -Wextra
        -Wwrite-strings -Wenum-compare
                    -Wno-c++98-compat
                    -Wno-c++98-compat-pedantic
                    -Wno-global-constructors
                    -Wno-exit-time-destructors
                    -Wno-documentation
                    -Wno-padded
                    -Wno-unused-parameter
                    -Wno-weak-vtables
                    -Wno-reserved-id-macro
                    -Wno-attributes -Weverything)
else()
    message(FATAL_ERROR "Unsupported or unknown compiler ${CMAKE_CXX_COMPILER_ID}")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "^.*Clang$")
    add_compile_options(-foptimize-sibling-calls)

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(ANOMALY_SANITIZE "none" CACHE STRING "Sanitizer not all options are available in all compiler versions")
        set_property(CACHE ANOMALY_SANITIZE PROPERTY STRINGS "" none address memory thread undefined leak)
        if(NOT ANOMALY_SANITIZE STREQUAL "none")
            add_compile_options(-fsanitize=${ANOMALY_SANITIZE} -fno-omit-frame-pointer)
            link_libraries(-fsanitize=${ANOMALY_SANITIZE})
        endif()
    endif()
endif()

option(ANOMALY_ENABLE_CLANG_TIDY "" ON)
if(ANOMALY_ENABLE_CLANG_TIDY)
    find_program(
        CLANG_TIDY_EXE
        NAMES "clang-tidy"
        DOC "Path to clang-tidy executable"
        )
    if(CLANG_TIDY_EXE)
        message(STATUS "Enable clang-tidy: ${CLANG_TIDY_EXE}")
        set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
    endif()
endif()

set(ANOMALY_PUBLIC_HEADERS
        src/anomaly/whisper/model.hpp
        src/anomaly/core/timeseries.hpp
        src/anomaly/core/lowess.hpp
        src/anomaly/core/stl.hpp
	src/anomaly/core/dsp.hpp
        src/anomaly/io/tsv.hpp
        )

add_library(anomaly
    ${ANOMALY_PUBLIC_HEADERS}
    src/anomaly/core/timeseries.cpp
    src/anomaly/core/stl.cpp
    src/anomaly/whisper/model.cpp
    )
target_link_libraries(anomaly
  INTERFACE ${CMAKE_DL_LIBS}
  Eigen3::Eigen
  Boost::date_time
  Boost::filesystem
  Boost::iostreams)
target_include_directories(anomaly PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

enable_testing()

add_executable(anomaly_test
    test/test_main.cpp
    test/matchers.hpp
    test/matchers.cpp
    test/lowess_test.cpp
    test/stl_test.cpp
    test/dsp_test.cpp
    test/catch.hpp)
target_link_libraries(anomaly_test anomaly)
target_include_directories(anomaly_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/test)
add_test(anomaly_test anomaly_test -r junit -o junit.xml -s)

add_executable(anomaly_example1
    src/anomaly/main.cpp # TODO move to examples dir
    )
target_link_libraries(anomaly_example1 anomaly)

add_executable(anomaly_experiments
    src/anomaly/experiments.cpp # TODO move to examples dir
    )
