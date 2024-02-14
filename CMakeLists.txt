cmake_minimum_required(VERSION 3.16)

project(
    PMC_HEURISTIC_EXTENDED
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

FILE(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*")
include_directories("${PROJECT_SOURCE_DIR}/include")

add_library(heulib ${SOURCES})

option(WITH_TEST "BUILD THE TEST CODE" ON)

if(WITH_TEST)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    )
    FetchContent_MakeAvailable(googletest)


    enable_testing()
    add_subdirectory(tests)
endif()

