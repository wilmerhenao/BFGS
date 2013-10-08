
CC       = g++
CFLAGS   = -O2 -fPIC -Wl,-V -Werror -pedantic-errors -Wall -Wextra -Wdouble-promotion -Wunused -Wuninitialized -Wstrict-overflow=5 -Wsuggest-attribute=const -Wshadow -Wconversion -Wsign-conversion -g -I../lib -ldl -lm -I../lib/qd -L$(QD_LIB) -lqd -std=c++11

all:

.cpp.o: 
	$(CC) $(CFLAGS) -c $< -o $@

matrixtest : nummatrix.o
	$(CC) $(CFLAGS) nummatrix.o -o matrixtest -L. -L../lib -lqd -llapack -lblas

clean:
	@rm -f *.o core
	@rm -f *.exe core test
