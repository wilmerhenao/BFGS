CC       = c++
CFLAGS   = -fPIC -g -ldl -I../lib -I../csrc/

TESTFUNOBJS = F10.o T10.o yurirosen.o F10_dd.o T10_dd.o yurirosen_dd.o

all: clean libtestfun.so

.c.o : 
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	@rm -f *.o core

libtestfun.so: $(TESTFUNOBJS)
	$(CC) -shared $(TESTFUNOBJS) -o ../lib/libtestfun.so
