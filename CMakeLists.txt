cmake_minimum_required(VERSION 3.25)

project(atomicsolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_EXTENSIONS False)

add_executable(app 
    src/main.cpp
)

target_include_directories(app PRIVATE include)