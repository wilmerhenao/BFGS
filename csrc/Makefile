
CC       = g++
CFLAGS   = -O2 -fPIC -Wl,-V -Werror -pedantic-errors -Wall -Wextra -Wdouble-promotion -Wunused -Wuninitialized -Wstrict-overflow=5 -Wsuggest-attribute=const -Wshadow -Wconversion -Wsign-conversion -g -I../lib -ldl -lm -I../lib/qd -lqd -std=c++11
TESTFUNINCLUDE = -I../testfunctions/

all:

.cpp.o: 
	$(CC) $(CFLAGS) $(TESTFUNINCLUDE) -c $< -o $@

mytest : test.o
	$(CC) $(CFLAGS) $(TESTFUNINCLUDE) test.o -o mytest -L. -L../lib -lqd 

clean:
	@rm -f *.o core
	@rm -f *.exe core test
