# driver: driver.o rkf.o
# 	lcc  driver.o rkf.o -o driver -lm
# makefile
# Tabs *must* be used for the indentations below;
# spaces cause make syntax errors.

CC=g++
#CFLAGS=-fast -xO4 -xdepend -xarch=v8plusa -xprefetch -xvector -xunroll=8 -fsimple=2 -xsafe=mem
LIBS=-lm
GSLLIBS=-lgsl -lgslcblas 

ode:		
		$(CC) -O3 -c -o pkdynamics_ppq.o pkdynamics_ppq.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o pkpd_dha.o pkpd_dha.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o main.o main.cpp $(LIBS) $(GSLLIBS)
		$(CC) $(CFLAGS) -o run_ppq_pk main.o pkpd_dha.o pkdynamics_ppq.o $(LIBS) $(GSLLIBS) 


clean:
		rm -f *.o core run_ppq_pk *~




