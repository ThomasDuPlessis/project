CC = g++
CLANG = clang++
FLAGS=-g -std=c++14 -I/include/glpk/include -lglpk -lm -Wextra -pedantic
PROTECT=protect.h protect.cc PASAQ.h PASAQ.cc lin_prog.cc lin_prog.h
MAIN=main.cc

all:
	$(CC) $(FLAGS) $(PROTECT) $(MAIN)
