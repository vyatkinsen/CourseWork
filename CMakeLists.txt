cmake_minimum_required(VERSION 3.17)
project(course)

set(CMAKE_CXX_STANDARD 14)
set(PROJECT_HEADERS
        rkf45.h
        quanc8.h
        matrix.h
        zeroin.h
        )
set(PROJECT_SOURCES
        main.cpp
        rkf45.cpp
        quanc8.cpp
        matrix.cpp
        zeroin.cpp
        )

add_executable(course ${PROJECT_SOURCES} ${PROJECT_HEADERS})