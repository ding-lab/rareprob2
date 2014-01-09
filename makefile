	
	CC=g++
all: rareprob	
rareprob: probrare.o global.o backward.o baum.o forward.o hmmrand.o \
	hmmutils.o nrutil.o sequence.o viterbi.o
	$(CC) -o $@ probrare.o global.o backward.o baum.o \
	forward.o hmmrand.o hmmutils.o nrutil.o sequence.o viterbi.o 

probrare.o: probrare.cpp hmm.h nrutil.h global.h 
	$(CC) -c probrare.cpp

global.o: global.cpp global.h hmm.h nrutil.h
	$(CC) -c global.cpp

backward.o: backward.c hmm.h
	$(CC) -c backward.c 

baum.o: baum.c hmm.h nrutil.h
	$(CC) -c baum.c 
forward.o: forward.c hmm.h
	$(CC) -c forward.c
hmmrand.o: hmmrand.c
	$(CC) -c hmmrand.c
hmmtuils.o: hmmtuils.c hmm.h nrutil.h
	$(CC) -c hmmutils.c 
nrutil.o: nrutil.c
	$(CC) -c nrutil.c
sequence.o: sequence.c hmm.h nrutil.h
	$(CC) -c sequence.c
	
viterbi.o: viterbi.c hmm.h nrutil.h -lm
	$(CC) -c viterbi.c
clean:
	-rm *.o
