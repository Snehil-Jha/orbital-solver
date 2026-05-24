# 1. Require a modern version of CMake (3.25+ is recommended for C++26 support)
cmake_minimum_required(VERSION 3.25)

# 2. Name your project and specify the language
project(atomicsolver LANGUAGES CXX)

# 3. Request a modern C++ standard (20, 23, or 26)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Optional but recommended: Disable compiler-specific extensions (e.g., gnu++23) 
# to ensure your code remains strictly cross-platform standard C++
set(CMAKE_CXX_EXTENSIONS False)

# 4. Tell CMake to build an executable named "my_app" from "main.cpp"
add_executable(app 
    src/main.cpp
)

target_include_directories(app PRIVATE include)