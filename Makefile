CC = mpicc

all: FirstWeek

FirstWeek:
	mkdir -p ./build/

	$(CC) ./Week_1/1.1_HelloWorld.c 	-o ./build/1.1_HelloWorld.exe
	$(CC) ./Week_1/1.2_ReciprocalSum.c 	-o ./build/1.2_ReciprocalSum.exe
	$(CC) ./Week_1/1.3_LoopSend.c 		-o ./build/1.3_LoopSend.exe



.SILENT clean:
	rm -rf */*.exe
	rm -rf */*.o