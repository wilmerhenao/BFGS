
CC       = c++
CFLAGS   = -O2 -fPIC -g -I../lib -ldl -lm -I../lib/qd -lqd
TESTFUNINCLUDE = -I../testfunctions/

LIBMATRIX = libmatrix_template.o
COBJS     = bfgs.o linesearch.o print_template.o quasinewt_updates_template.o randnums_template.o
COBJSALL  = $(COBJS) $(LIBMATRIX)

LIBMATRIX_DD = libmatrix_template.o
COBJS_DD     = bfgs_dd.o linesearch_dd.o  print_template.o quasinewt_updates_template.o randnums_template.o
COBJSALL_DD  = $(COBJS_DD) $(LIBMATRIX_DD)

all: subsystem libhansoc.so libhansoc_dd.so

subsystem:
	cd ../testfunctions/ && $(MAKE)

.c.o : 
	$(CC) $(CFLAGS) $(TESTFUNINCLUDE) -c $< -o $@

mytest : libhansoc.so libhansoc_dd.so test.o
	$(CC) $(CFLAGS) $(TESTFUNINCLUDE) test.o -o mytest -L. -L../lib -lhansoc -lhansoc_dd -ltestfun -lqd 

clean:
	@rm -f *.o core
	@rm -f *.exe core test

libhansoc.so: $(COBJSALL)
	$(CC) -shared $(COBJSALL) -o ../lib/libhansoc.so 

libhansoc_dd.so: $(COBJSALL_DD)
	$(CC) -shared $(COBJSALL_DD) -o ../lib/libhansoc_dd.so 
