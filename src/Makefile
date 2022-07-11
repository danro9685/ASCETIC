SRCS = agony.cpp contract.cpp read.cpp balance.cpp ccomp.cpp
OBJS = $(SRCS:.cpp=.o)


CC = g++
LDFLAGS= -lm -O3 -L/opt/local/lib
CFLAGS= -g -O3 -Wall -Wextra -isystem/opt/local/include/ #-Werror

all : agony 

agony: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)


%.o : %.cpp
	 $(CC) $(CFLAGS) -c $*.cpp

clean:
	-rm *.o agony

zip:
	zip agony.zip *.cpp *.h Makefile

depend:
	makedepend -Y -- $(CFLAGS) -- $(SRCS)
# DO NOT DELETE

agony.o: defines.h graph.h queue.h circulation.h tree.h
contract.o: defines.h graph.h queue.h circulation.h tree.h
read.o: circulation.h defines.h graph.h queue.h tree.h
balance.o: circulation.h defines.h graph.h queue.h tree.h
ccomp.o: circulation.h defines.h graph.h queue.h tree.h
