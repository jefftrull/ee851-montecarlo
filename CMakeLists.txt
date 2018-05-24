cmake_minimum_required( VERSION 3.5.0 )

# Build control file for an ancient semiconductor physics project

add_executable( monte monte.cpp )
target_compile_options( monte PUBLIC -Wall -Wextra -Werror )
target_link_libraries( monte -lm )
