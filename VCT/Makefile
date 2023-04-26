OBJECTS = CPM_dH.o \
	cellmoves.o \
	CM.o \
	init.o \
	mt.o \
	mylib.o \
	cpmfem.o \
	write.o

# include standard C library
LDFLAGS += -lc 
# include standard math library
LDFLAGS += -lm
CFLAGS +=-fPIC


libcpmfem.a: $(OBJECTS)
	ar rcs $@ $^
%.o : %.c, %.h
	gcc -c $(CFLAGS) $< -o -ld -lm $@

clean:
	rm *.o
	rm ./output/*.out
	rm ./imgs/*.png
	rm cpmfem
	rm *.a

