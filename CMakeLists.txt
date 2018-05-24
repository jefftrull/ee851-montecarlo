cmake_minimum_required( VERSION 3.5.0 )

# Build control file for an ancient semiconductor physics project

set( CMAKE_CXX_STANDARD 17 )

add_executable( monte monte.cpp )
target_compile_options( monte PUBLIC -Wall -Wextra -Werror )
target_link_libraries( monte -lm )
