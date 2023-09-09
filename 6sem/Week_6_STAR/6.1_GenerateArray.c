#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char* argv[]) {

    if(argc < 2) {
        printf("Usage: %s [N]\n", argv[0]);
        return 1;
    }

    int arraySize = atoi(argv[1]);

    time_t t;
    srand((unsigned)time(&t));

    printf("%d\n", arraySize);
    for(int i = 0; i < arraySize; i++) {
        printf("%d\n", rand() % arraySize - arraySize / 2);
    }


    return 0;
}