#OPT    += -DDIRECT_RAY_INTEGRATION

EXEC   = attenuate_py.so

SRCS   = attenuate.c


OBJS   = $(SRCS:.c=.o)
INCL   = 


#CFLAGS =  -shared -m32 $(OPT)
#CFLAGS =  -bundle $(OPT)
CFLAGS =  -shared -fPIC -O3 -g  #  -Wall
#SDK    = /Developer/SDKs/MacOSX10.4u.sdk
#CFLAGS  =  -bundle $(OPT) -isysroot $(SDK) -arch ppc -arch ppc64 -arch i386
#LDFLAGS = -isysroot $(SDK) -syslibroot,$(SDK)

LIBS   = -lm

CC     =  gcc 

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(SRCS) 



clean:
	rm $(OBJS) $(EXEC)










