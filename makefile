.PHONY = all clean main pLocal pEvol pHybrid

CC = gcc

LINKERFLAG = -lm

SRCS := $(wildcard *.c)
BINS := main pLocal pEvol pHybrid

all: main pLocal pEvol pHybrid
	cp ${BINS} ./testes
	mkdir ./testes/out

main: main.o
	${CC} ${LINKERFLAG} $< -o $@

pLocal: pLocal.o util.o
	${CC} ${LINKERFLAG} $< util.o -o $@

pEvol: pEvol.o util.o
	${CC} ${LINKERFLAG} $< util.o -o $@

pHybrid: pHybrid.o util.o
	${CC} ${LINKERFLAG} $< util.o -o $@

%.o: %.c
	@echo "Creating object.."
	${CC} -c $<

clean:
	@echo "Cleaning up..."
	rm -rvf *.o ${BINS}
	