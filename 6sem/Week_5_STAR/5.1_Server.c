#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>


 
int main(int argc, char *argv[]) {
   int size, rank;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   MPI_Comm client;
   MPI_Status status;
   char portname[MPI_MAX_PORT_NAME];

   if(MPI_Open_port(MPI_INFO_NULL, portname) != MPI_SUCCESS) {
      printf("[server] could not open port\n");
      return 1;
   }
   printf("[server] opened server with portname: %s\n", portname);

#ifndef MPICH
   MPI_Info info;
   MPI_Info_create(&info);
   MPI_Info_set(info, "ompi_local_scope", "true");
   MPI_Publish_name("my_server", info, portname);
#endif

   if(MPI_Comm_accept(portname, MPI_INFO_NULL, 0, MPI_COMM_SELF, &client) != MPI_SUCCESS)
      printf("[server] could not accept a client\n");

   printf("[server] client connected\n");

   char message[32];
   MPI_Recv(message, 32, MPI_CHAR, MPI_ANY_SOURCE, 0, client, &status);
   printf("[server] received message = \"%s\"\n", message);

#ifndef MPICH
   sleep(1);
   MPI_Unpublish_name("my_server", MPI_INFO_NULL, portname);
#endif

   MPI_Comm_free(&client);
   MPI_Close_port(portname);
   printf("[server] closed\n");

   MPI_Finalize();

   return 0;
}