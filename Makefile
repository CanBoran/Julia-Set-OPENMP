all: julia.cpp
	mpicxx -O3  julia.cpp -o julia
