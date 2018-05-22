cmake_minimum_required( VERSION 3.5.0 )

# Build control file for an ancient semiconductor physics project

add_executable( monte monte.c )
target_link_libraries( monte -lm )
