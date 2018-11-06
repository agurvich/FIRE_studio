# Need to make this file twice. Once like this and once with the two
# comments removed.

EXEC  = HsmlAndProject.so
#EXEC  = HsmlAndProject_perbc.so

OBJS   = main.o tree.o peano.o allvars.o

#OPTIONS += -DPERIODIC

INCL   = proto.h allvars.h Makefile

CFLAGS =  $(OPTIONS) -fPIC -Wall -O2 
LNKCMD =  ld -shared

LIBS   =  -lm 

CC     = cc#clang 

$(EXEC): $(OBJS)
	 $(LNKCMD)  $(OBJS) $(LIBS)  -o $(EXEC)

$(OBJS): $(INCL)

clean:
	 rm -f $(OBJS) $(EXEC)
