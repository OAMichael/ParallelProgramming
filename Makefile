MPICC = mpicc
CC = gcc
CXX = g++

all: 6sem_Lesson_1 6sem_Lesson_2 6sem_Lesson_3 6sem_Lesson_4 6sem_Lesson_5 6sem_Lesson_6 6sem_Lesson_7 6sem_Lesson_8 7sem_Lesson_1

6sem_Lesson_1:
	mkdir -p ./build/

	$(MPICC) ./6sem/Week_1/1.1_HelloWorld.c     -o ./build/1.1_HelloWorld.exe
	$(MPICC) ./6sem/Week_1/1.2_ReciprocalSum.c  -o ./build/1.2_ReciprocalSum.exe
	$(MPICC) ./6sem/Week_1/1.3_LoopSend.c       -o ./build/1.3_LoopSend.exe


6sem_Lesson_2:
	mkdir -p ./build/

	$(MPICC) ./6sem/Week_2_STAR/2.STAR_Exp.c -O3 -lm -lgmp -o ./build/2.STAR_Exp.exe


6sem_Lesson_3:
	mkdir -p ./build/

	$(MPICC) ./6sem/Week_3/3.1_SendFunctions.c      -o ./build/3.1_SendFunctions.exe
	$(MPICC) ./6sem/Week_3/3.2_ReciprocalSumComm.c  -o ./build/3.2_ReciprocalSumComm.exe


6sem_Lesson_4:
	mkdir -p ./build/

	$(MPICC) ./6sem/Lab1/Lab1.c -O3 -lm -o ./build/Lab1.exe


6sem_Lesson_5:
	mkdir -p ./build/

	$(MPICC) ./6sem/Week_5_STAR/5.1_Server.c        -o ./build/5.1_Server.exe
	$(MPICC) ./6sem/Week_5_STAR/5.1_Client.c        -o ./build/5.1_Client.exe
	$(MPICC) ./6sem/Week_5_STAR/5.1_ClientServer.c  -o ./build/5.1_ClientServer.exe
	$(MPICC) ./6sem/Week_5_STAR/5.2_OneSideComm.c   -o ./build/5.2_OneSideComm.exe
	$(MPICC) ./6sem/Week_5_STAR/5.3_FilesWithMPI.c  -o ./build/5.3_FilesWithMPI.exe


6sem_Lesson_6:
	mkdir -p ./build/

	$(MPICC) ./6sem/Week_6_STAR/6.1_Sort.c           -o ./build/6.1_Sort.exe
	$(CC)    ./6sem/Week_6_STAR/6.1_GenerateArray.c  -o ./build/6.1_GenerateArray.exe


6sem_Lesson_7:
	mkdir -p ./build/

	$(CC)    ./6sem/Week_7/7.1_HelloWorldThreads.c      -lpthread -o ./build/7.1_HelloWorldThreads.exe
	$(CC)    ./6sem/Week_7/7.2_ReciprocalSumThreads.c   -lpthread -o ./build/7.2_ReciprocalSumThreads.exe
	$(CC)    ./6sem/Week_7/7.3_VariableAccess.c         -lpthread -o ./build/7.3_VariableAccess.exe


6sem_Lesson_8:
	$(CXX)   ./6sem/Week_8/8.1_Integral.cpp         -O3 -lpthread -o ./build/8.1_Integral.exe


7sem_Lesson_1:
	$(CC)    ./7sem/Week_1/1.1_HelloWorldOMP.c 		-fopenmp -o ./build/1.1_HelloWorldOMP.exe
	$(CC)    ./7sem/Week_1/1.2_ReciprocalSumOMP.c 	-fopenmp -o ./build/1.2_ReciprocalSumOMP.exe
	$(CC)    ./7sem/Week_1/1.3_VariableAccessOMP.c 	-fopenmp -o ./build/1.3_VariableAccessOMP.exe

.SILENT clean:
	rm -rf */*.exe
	rm -rf */*.o