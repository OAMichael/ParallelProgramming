CC = mpicc

all: FirstWeek SecondWeekStar

FirstWeek:
	mkdir -p ./build/

	$(CC) ./Week_1/1.1_HelloWorld.c 	-o ./build/1.1_HelloWorld.exe
	$(CC) ./Week_1/1.2_ReciprocalSum.c 	-o ./build/1.2_ReciprocalSum.exe
	$(CC) ./Week_1/1.3_LoopSend.c 		-o ./build/1.3_LoopSend.exe


SecondWeekStar:
	mkdir -p ./build/
	$(CC) ./Week_2_STAR/2.STAR_Exp.c -lgmp -o ./build/2.STAR_Exp.exe

.SILENT clean:
	rm -rf */*.exe
	rm -rf */*.o