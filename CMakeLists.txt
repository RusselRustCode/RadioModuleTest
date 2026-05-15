cmake_minimum_required(VERSION 3.14)
project(dsp_project CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/heads/main.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(fft_main fft.cpp)
target_compile_options(fft_main PRIVATE -Wall -Wextra)

add_executable(radio_module radio_module.cpp)
target_compile_options(radio_module PRIVATE -Wall -Wextra)


# Тесты
enable_testing()

add_executable(radio_module_tests radio_module_tests.cpp)
target_include_directories(radio_module_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(radio_module_tests PRIVATE GTest::gtest_main)
target_compile_options(radio_module_tests PRIVATE -Wall -Wextra)

include(GoogleTest)
gtest_discover_tests(radio_module_tests)
