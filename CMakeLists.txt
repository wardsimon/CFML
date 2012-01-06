cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# By default the built is in CONSOLE mode.
option(GUI "Build also crysfgl and wcrysfml library" ON)

# The project is declared. It is a Fortran project.
project(crysfml Fortran)

# The subdirectories of the project.
add_subdirectory(Src)
