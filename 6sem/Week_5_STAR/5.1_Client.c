#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>


int main(int argc, char *argv[]) {
   int size, rank;
 
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   MPI_Comm server;
   char portname[MPI_MAX_PORT_NAME];

#ifndef MPICH
   if(MPI_Lookup_name("my_server", MPI_INFO_NULL, portname) != MPI_SUCCESS) {
      printf("[client] could not find server\n");
      return 1;
   }
   printf("[client] found server with portname:  %s\n", portname);
#else
   strcpy(portname, argv[1]);
#endif

   if(MPI_Comm_connect(portname, MPI_INFO_NULL, 0, MPI_COMM_SELF, &server) != MPI_SUCCESS) {
      printf("[client] could not connect to the server\n");
      return 1;
   }
   

   char message[32] = "Hello!";

   printf("[client] sending message to server...\n");
   if(MPI_Send(message, 32, MPI_CHAR, 0, 0, server) != MPI_SUCCESS)
      printf("[client] could not send data to the server\n");   

   if(MPI_Comm_free(&server) != MPI_SUCCESS)
      printf("[client] could not disconnect from server\n"); 
   else
      printf("[client] disconnected from server\n"); 

   MPI_Finalize();

   return 0;
}