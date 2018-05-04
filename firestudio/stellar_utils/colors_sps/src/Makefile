CC       =  cc            # sets the C-compiler

OPTIMIZE =  #-O2 #-Wall 
OPT1     =  #-xarch=v9

OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2) -fno-common

EXEC   = colors.so


OBJS   = main.o \
	 	 pl.o 

INCL   = pl.h

SDK    = /Developer/SDKs/MacOSX10.4u.sdk
CFLAGS =   $(OPTIONS) -I/sw/include
#CFLAGS  =  -bundle $(OPTIONS) -I/sw/include -isysroot $(SDK) -arch i386 
LDFLAGS = -isysroot $(SDK) -syslibroot,$(SDK)

#LIBS   = -L/sw/include -gsl -gsl_cblas -gsl_interp -lm
LIBS   = -L/sw/lib -lm -lgsl -lgslcblas

$(EXEC): $(OBJS) 
	$(CC) $(CFLAGS) -flat_namespace $(OBJS) $(LIBS)  -o $(EXEC)

#	$(CC) -bundle -flat_namespace -undefined suppress $(OBJS) $(LIBS)  -o $(EXEC)


#	ld  -64  -G $(OBJS) $(LIBS)  -o $(EXEC)
#	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  


$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

