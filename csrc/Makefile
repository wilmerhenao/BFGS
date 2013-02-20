
CC       = c++
CFLAGS   = -O2 -fPIC -g -I../lib -ldl -lm -I../lib/qd -lqd
TESTFUNINCLUDE = -I../testfunctions/

all: subsystem

subsystem:
	cd ../testfunctions/ && $(MAKE)

.cpp.o: 
	$(CC) $(CFLAGS) $(TESTFUNINCLUDE) -c $< -o $@

mytest : test.o
	$(CC) $(CFLAGS) $(TESTFUNINCLUDE) test.o -o mytest -L. -L../lib -ltestfun -lqd 

clean:
	@rm -f *.o core
	@rm -f *.exe core test
