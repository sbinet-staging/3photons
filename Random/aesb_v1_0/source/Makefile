
CC=nvcc

CFLAGS=-arch=sm_20 -maxrregcount 32 -O2

SRC =	*.cu

OBJ =	$(SRC:.cu=.o)

OUT =	../lib/libprand.a

$(OUT): $(OBJ)
	ar rcs $(OUT) $(OBJ)

$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c $(SRC)

clean: 
	rm $(OBJ)
