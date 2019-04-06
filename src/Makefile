CSRCS     = $(wildcard *.cpp) $(wildcard *.c)
CHDRS     = $(wildcard *.h)
COBJS     = $(addsuffix .o, $(basename $(CSRCS)))

CXX       = g++
ECHO      = /bin/echo

CFLAGS    = -g -Wall -std=c++11 -DDEBUG_MODE
CFLAGS    = -O3 -Wall -std=c++11

EXEC      = fm


all: $(COBJS)
	@echo "> building $(EXEC)..."
	@$(CXX) $(CFLAGS) $(COBJS) -o $(EXEC)

%.o : %.cpp
	@$(ECHO) "> compiling: $<"
	@$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	@rm -f $(COBJS)
