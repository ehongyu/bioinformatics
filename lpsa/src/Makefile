CC = CC -O3

OBJECT = lpsa.o class_lp.o class_pdb.o class_que.o class_ran.o grid.o envi.o
lpsa: $(OBJECT)
	$(CC) $(OBJECT) -lm -o lpsa
lpsa.o: lpsa.c++ head.h
	$(CC) lpsa.c++ -c
class_lp.o: class_lp.c++ head.h grid.h
	$(CC) class_lp.c++ -c
class_pdb.o: class_pdb.c++ head.h
	$(CC) class_pdb.c++ -c
class_que.o: class_que.c++ head.h
	$(CC) class_que.c++ -c
class_ran.o: class_ran.c++ head.h
	$(CC) class_ran.c++ -c
grid.o: grid.c++ head.h
	$(CC) grid.c++ -c
envi.o: envi.c++ head.h
	$(CC) envi.c++ -c
