
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%  BFGS tools Copyright (C) 2013  Wilmer Henao
#%%  This program is free software: you can redistribute it and/or modify
#%%  it under the terms of the GNU General Public License as published by
#%%  the Free Software Foundation, either version 3 of the License, or
#%%  (at your option) any later version.
#%%
#%%  This program is distributed in the hope that it will be useful,
#%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
#%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#%%  GNU General Public License for more details.
#%%
#%%  You should have received a copy of the GNU General Public License
#%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC       = g++
CFLAGS   = -O2 -fPIC -Wl,-V -Werror -pedantic-errors -Wall -Wextra -Wdouble-promotion -Wunused -Wuninitialized -Wstrict-overflow=5 -Wsuggest-attribute=const -Wshadow -Wconversion -Wsign-conversion -g -I../lib -ldl -lm -I../lib/qd -lqd -llapack -lblas -std=c++11
TESTFUNINCLUDE = -I../testfunctions/

all:

.cpp.o: 
	$(CC) $(CFLAGS) $(TESTFUNINCLUDE) -c $< -o $@

mytest : test.o
	$(CC) $(CFLAGS) $(TESTFUNINCLUDE) test.o -o mytest -L. -L../lib -lqd -llapack -lblas

clean:
	@rm -f *.o core
	@rm -f *.exe core test
