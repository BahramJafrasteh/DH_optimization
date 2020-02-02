
include make.linux

all: DHLocOpt

DHLocOpt: DH_OpAlg.o DHLocOpt.o DH_OpObj.o
	$(LD) ${XLINKERFLAGS} -o DHLocOpt DHLocOpt.o DH_OpAlg.o DH_OpObj.o $(LDFLAGS)
DH_OpAlg.o: DH_OpAlg.cpp DH_OpAlg.h
	$(CC) -c DH_OpAlg.cpp -o DH_OpAlg.o $(CCFLAGS)
DHLocOpt.o: DHLocOpt.cpp DHLocOpt.h DH_OpAlg.cpp DH_OpAlg.h DH_OpObj.cpp DH_OpObj.h CInterface.h
	$(CC) -c DHLocOpt.cpp -o DHLocOpt.o $(CCFLAGS)
DH_OpObj.o: DH_OpObj.cpp DH_OpObj.h lapack.h CInterface.h
	$(CC) -c DH_OpObj.cpp -o DH_OpObj.o $(CCFLAGS)
clean:
	$(RM) *.o
