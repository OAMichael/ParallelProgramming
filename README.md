# Parallel Programming

### There are tasks for course of parallel programming. They are sorted by weeks.

## How to build and run
Being in main folder you can type 
```console
make
```
to create `build` folder and build everything into it. To run one of executables you have to go to folder build and run
```console
mpirun -np [num_procs] [executable]
```
For example:
```console
mpirun -np 4 ./1.1_HelloWorld.exe
```


***Note:*** for 1.2_ReciprocalSum.exe you have to specify which number you are summing up to:
```console
mpirun -np 4 ./1.2_ReciprocalSum.exe 1000000
```