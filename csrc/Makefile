
CC       = c++
CFLAGS   = -O2 -fPIC -g -I../lib -ldl -lm -I../lib/qd -lqd
TESTFUNINCLUDE = -I../testfunctions/

#LIBMATRIX = libmatrix_template.o
#COBJS     = bfgs_template.o print_template.o quasinewt_updates_template.o randnums_template.o
#COBJSALL  = $(COBJS) $(LIBMATRIX)

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

#libhansoc.so: $(COBJSALL)
#	$(CC) -shared $(COBJSALL) -o ../lib/libhansoc.so 
