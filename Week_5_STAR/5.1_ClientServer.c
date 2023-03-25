#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

 
 
int main(int argc, char *argv[]){
   int size, rank;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // Server part
   if(rank == 0) {

      MPI_Comm client;
      MPI_Status status;
      char portname[MPI_MAX_PORT_NAME];

      MPI_Open_port(MPI_INFO_NULL, portname);
      printf("[server] opened server with portname: %s\n", portname);

      MPI_Info info;
      MPI_Info_create(&info);
      MPI_Info_set(info, "ompi_local_scope", "true");
      MPI_Publish_name("my_server", info, portname);


      if(MPI_Comm_accept(portname, MPI_INFO_NULL, 0, MPI_COMM_SELF, &client) != MPI_SUCCESS)
         printf("[server] could not accept a client\n");

      printf("[server] client connected\n");

      char message[32];
      MPI_Recv(message, 32, MPI_CHAR, MPI_ANY_SOURCE, 0, client, &status);
      printf("[server] received message = \"%s\"\n", message);

      sleep(1);
      MPI_Unpublish_name("my_server", MPI_INFO_NULL, portname);

      MPI_Comm_free(&client);
      MPI_Close_port(portname);
      printf("[server] closed\n");
   }
   // Client part
   if(rank == 1) {

      // Imitate delay as for real server-client
      sleep(1);

      MPI_Comm server;
      char portname[MPI_MAX_PORT_NAME];
      if(MPI_Lookup_name("my_server", MPI_INFO_NULL, portname) != MPI_SUCCESS) {
         printf("[client] could not find server\n");
         return 1;
      }
      printf("[client] found server with portname:  %s\n", portname);

      printf("[client] sending message to server...\n");
      if(MPI_Comm_connect(portname, MPI_INFO_NULL, 0, MPI_COMM_SELF, &server) != MPI_SUCCESS) {
         printf("[client] could not connect to the server\n");
         return 1;
      }
      

      char message[32] = "Hello!";

      if(MPI_Send(message, 32, MPI_CHAR, 0, 0, server) != MPI_SUCCESS)
         printf("[client] could not send data to the server\n");   

      if(MPI_Comm_free(&server) != MPI_SUCCESS)
         printf("[client] could not disconnect from server\n"); 
      else
         printf("[client] disconnected from server\n");  
      
   }


   MPI_Finalize();

   return 0;
}