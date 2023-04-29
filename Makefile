MPICC = mpicc
CC = gcc

all: FirstWeek SecondWeekStar ThirdWeek Lab_1 FifthWeek SixthWeek SeventhWeek

FirstWeek:
	mkdir -p ./build/

	$(MPICC) ./Week_1/1.1_HelloWorld.c 	    -o ./build/1.1_HelloWorld.exe
	$(MPICC) ./Week_1/1.2_ReciprocalSum.c 	-o ./build/1.2_ReciprocalSum.exe
	$(MPICC) ./Week_1/1.3_LoopSend.c 		-o ./build/1.3_LoopSend.exe


SecondWeekStar:
	mkdir -p ./build/

	$(MPICC) ./Week_2_STAR/2.STAR_Exp.c -O3 -lm -lgmp -o ./build/2.STAR_Exp.exe


ThirdWeek:
	mkdir -p ./build/
	
	$(MPICC) ./Week_3/3.1_SendFunctions.c  	    -o ./build/3.1_SendFunctions.exe
	$(MPICC) ./Week_3/3.2_ReciprocalSumComm.c 	-o ./build/3.2_ReciprocalSumComm.exe


Lab_1:
	mkdir -p ./build/

	$(MPICC) ./Lab1/Lab1.c -O3 -lm -o ./build/Lab1.exe


FifthWeek:
	mkdir -p ./build/

	$(MPICC) ./Week_5_STAR/5.1_Server.c  		-o ./build/5.1_Server.exe
	$(MPICC) ./Week_5_STAR/5.1_Client.c  		-o ./build/5.1_Client.exe
	$(MPICC) ./Week_5_STAR/5.1_ClientServer.c   -o ./build/5.1_ClientServer.exe
	$(MPICC) ./Week_5_STAR/5.2_OneSideComm.c  	-o ./build/5.2_OneSideComm.exe
	$(MPICC) ./Week_5_STAR/5.3_FilesWithMPI.c   -o ./build/5.3_FilesWithMPI.exe


SixthWeek:
	mkdir -p ./build/

	$(MPICC) ./Week_6_STAR/6.1_Sort.c  		    -o ./build/6.1_Sort.exe
	$(CC)    ./Week_6_STAR/6.1_GenerateArray.c  -o ./build/6.1_GenerateArray.exe


SeventhWeek:
	mkdir -p ./build/

	$(CC) 	 ./Week_7/7.1_HelloWorldThreads.c  		-lpthread -o ./build/7.1_HelloWorldThreads.exe
	$(CC) 	 ./Week_7/7.2_ReciprocalSumThreads.c  	-lpthread -o ./build/7.2_ReciprocalSumThreads.exe
	$(CC) 	 ./Week_7/7.3_VariableAccess.c  		-lpthread -o ./build/7.3_VariableAccess.exe

.SILENT clean:
	rm -rf */*.exe
	rm -rf */*.o