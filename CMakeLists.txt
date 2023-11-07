cmake_minimum_required(VERSION 3.27.7)
project(geo)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi -pedantic -Wall -Wextra -Werror")
add_executable(
	geo geo.c
	main.c
	wgs84.h wgs84.c)
