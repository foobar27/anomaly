cmake_minimum_required(VERSION 3.17)

project(anomaly)

find_package(Boost 1.72 REQUIRED COMPONENTS
             date_time filesystem iostreams)

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
                        )
add_library(anomaly_headers INTERFACE)
target_include_directories(anomaly_headers INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>
)

add_executable(anomaly
                        ${ANOMALY_PUBLIC_HEADERS}
                        src/anomaly/whisper/model.cpp
                        src/anomaly/main.cpp
                        )
target_link_libraries(anomaly PUBLIC anomaly_headers
                                       INTERFACE ${CMAKE_DL_LIBS}
                                       Boost::date_time
                                       Boost::filesystem
                                       Boost::iostreams)

set_target_properties(anomaly PROPERTIES PUBLIC_HEADER "${ANOMALY_PUBLIC_HEADERS}")

# enable_testing()

# Test library
#add_library(anomaly_test_lib SHARED test/extension.cpp test/Interfaces.hpp)
#target_link_libraries(anomaly_test_lib PRIVATE extension_system_headers)
# target_link_libraries(anomaly PRIVATE Boost::boost)

# Test program
#add_executable(anomaly_test test/main.cpp test/Interfaces.hpp test/catch.hpp)
#target_link_libraries(anomaly_test anomaly)
#add_test(anomaly_test anomaly_test -r junit -o junit.xml -s)

# Examples
## Example 1
#add_library(anomaly_example1_extension SHARED
#    examples/example1/Extension.cpp
#    examples/example1/Interface.hpp)
#target_link_libraries(anomaly_example1_extension PRIVATE anomaly_headers)

#add_executable(anomaly_example1 examples/example1/main.cpp examples/example1/Interface.hpp)
#target_link_libraries(anomaly_example1 extension_system)
