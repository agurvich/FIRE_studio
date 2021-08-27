
EXEC   = starhsml.so

SRCS   = main.c ngbtree3d.c  \
	 nrutil.c selectb.c


OBJS   = $(SRCS:.c=.o)
INCL   = nrutil.h


CFLAGS =  -shared -m32

LIBS   =   

CC     =  gcc 

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(INCL) 



clean:
	rm $(OBJS) $(EXEC)










