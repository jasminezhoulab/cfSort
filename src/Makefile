CC=g++ -w -O3 -static
#CFLAGS=  
getoptDIR=./my_getopt-1.5
INCLUDE= -I/u/local/apps/boost/1_55_0/gcc-4.4.7/include -I/u/local/apps/zlib/1.2.8/gcc-4.4.7/include
OBJECTS1= cfsort.o matrix.o utils.o data_types.o my_getopt.o
LIBS= -L/u/local/apps/boost/1_55_0/gcc-4.4.7/lib -L/u/local/apps/zlib/1.2.8/gcc-4.4.7/lib -L/usr/lib64 -lboost_iostreams -lz -lm
PROG1=cfsort

all: ${PROG1}

cfsort.o: cfsort.cpp utils.h data_types.h matrix.h
	${CC} ${CFLAGS} ${INCLUDE} -c cfsort.cpp -o cfsort.o

%.o: %.cpp %.h
	${CC} ${CFLAGS} ${INCLUDE} -c $< -o $@

${PROG1}: ${OBJECTS1}
	${CC} ${CFLAGS} ${OBJECTS1} ${LIBS} -o ${PROG1}


my_getopt.o: $(getoptDIR)/my_getopt.c $(getoptDIR)/my_getopt.h $(getoptDIR)/getopt.h
	cd $(getoptDIR)
	$(CC) $(CFLAGS) -I$(getoptDIR) -c $(getoptDIR)/my_getopt.c -o my_getopt.o

# --- remove binary and executable files
clean:
	rm -f ${OBJECTS1}
	rm -f ${PROG1}
