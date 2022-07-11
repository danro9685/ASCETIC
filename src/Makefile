SRCS = agony.cpp contract.cpp read.cpp balance.cpp ccomp.cpp
OBJS = $(SRCS:.cpp=.o)


CC = g++
LDFLAGS= -lm -O3 #-L/opt/local/lib
CFLAGS= -g -O3 -MMD -Wall -Wextra #-I/opt/local/include/ #-Werror

all : agony 

agony: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)


%.o : %.cpp
	 $(CC) $(CFLAGS) -c $*.cpp

clean:
	-rm *.o agony

zip:
	zip agony.zip *.cpp *.h Makefile


-include $(SRCS:%.cpp=%.d)
