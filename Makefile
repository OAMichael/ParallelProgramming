CC = mpicc

all: FirstWeek SecondWeekStar ThirdWeek Lab_1

FirstWeek:
	mkdir -p ./build/

	$(CC) ./Week_1/1.1_HelloWorld.c 	-o ./build/1.1_HelloWorld.exe
	$(CC) ./Week_1/1.2_ReciprocalSum.c 	-o ./build/1.2_ReciprocalSum.exe
	$(CC) ./Week_1/1.3_LoopSend.c 		-o ./build/1.3_LoopSend.exe


SecondWeekStar:
	mkdir -p ./build/

	$(CC) ./Week_2_STAR/2.STAR_Exp.c -O3 -lm -lgmp -o ./build/2.STAR_Exp.exe


ThirdWeek:
	mkdir -p ./build/
	
	$(CC) ./Week_3/3.1_SendFunctions.c  	-o ./build/3.1_SendFunctions.exe
	$(CC) ./Week_3/3.2_ReciprocalSumComm.c 	-o ./build/3.2_ReciprocalSumComm.exe


Lab_1:
	mkdir -p ./build/

	$(CC) ./Lab1/Lab1.c -O3 -lm -o ./build/Lab1.exe

.SILENT clean:
	rm -rf */*.exe
	rm -rf */*.o