/* Mandelbrot set program.
 */

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<time.h>
#include<string.h>

#define GO_HOME 0
#define DO_THIS_ROW 1
#define ROW_DONE 2
#define SENDING_ROW 3

int mandelbrot(double cx, double cy, int max_its);

int main(int argc, char * argv[]) {

  int its = 200;
  double startx = -2;
  double starty = -2;
  double endx  = 2;
  double endy = 2;
  int cols = 512;
  int rows = 384;

  int i = 1;
  char *arg = NULL;
  //parse arguments
  while(i < argc) {
    arg = argv[i++];
    if(0 == strcmp(arg, "--startx")) {
      startx=atof(argv[i++]);
      //printf("startx = %f\n", startx);
    } else if(0 == strcmp(arg, "--starty")) {
      starty=atof(argv[i++]);
      //printf("starty = %f\n",starty);
    } else if(0 == strcmp(arg, "--endx")) {
      endx = atof(argv[i++]);
      //printf("endx = %f\n",endx);
    } else if(0 == strcmp(arg, "--endy")) {
      endy = atof(argv[i++]);
      //printf("endy = %f\n",endy);
    } else if(0 == strcmp(arg, "--its")) {
      its = atoi(argv[i++]);
    } else if(0 == strcmp(arg, "--rows")) {
      rows = atoi(argv[i++]);
    } else if(0 == strcmp(arg, "--cols")) {
      cols = atoi(argv[i++]);
    } else {
      printf("Error: Unknown parameter: %s\n", arg);
      return -1;
    }
  }

  MPI_Init(&argc, &argv);
  int myRank;
  int P;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  double startTime;
  double endTime;
  int *gatheredPicture = NULL;
  
  if(myRank == 0)
  {
    gatheredPicture = malloc(sizeof(int)*rows*cols);
  }
  

  int *myRow = malloc(sizeof(int)*rows*cols);
  int col;
  double cx; double cy;
  int color;

  int flag = -1;
  MPI_Status status;
  int source;
  MPI_Barrier(MPI_COMM_WORLD);
  startTime = MPI_Wtime();
  
  //Master
  if (myRank == 0) {
    int sendRow = 0;
    int recvRow = 0;
    int rowsFilled = 0;
    int goHome = GO_HOME;

    
    //First, send rows to all the processors
    for ( i = 1; i < P; i++)
    {
      MPI_Send(&sendRow, 1, MPI_INT, i, DO_THIS_ROW, MPI_COMM_WORLD); sendRow++;
    }

    //Loop will send out row numbers to the other processors to process. 
    while( rowsFilled < rows ) {
      MPI_Iprobe(MPI_ANY_SOURCE, ROW_DONE,  MPI_COMM_WORLD, &flag, &status);
      if(flag == 1) {
        source = status.MPI_SOURCE;
	//Get the row number
	MPI_Recv(&recvRow, 1, MPI_INT, source, ROW_DONE, MPI_COMM_WORLD, &status);
	//Get the row
	MPI_Recv(&gatheredPicture[recvRow*cols], cols, MPI_INT, source, SENDING_ROW, MPI_COMM_WORLD, &status);
	rowsFilled++;

	// send the next row, if there is one.
	if(sendRow < rows) {
	  MPI_Send(&sendRow, 1, MPI_INT, source, DO_THIS_ROW, MPI_COMM_WORLD); sendRow++;
	} else { //Otherwise, tell the worker to go home
	  MPI_Send(&goHome, 1, MPI_INT, source, GO_HOME, MPI_COMM_WORLD);
	}


      }
    }

  }

  //Apprentice
  else {
    int message = -1;
    int currentRow;
    MPI_Recv(&currentRow, 1, MPI_INT, 0, DO_THIS_ROW, MPI_COMM_WORLD, &status);
    while(message != GO_HOME) {
        cy = starty + (endy-starty)*(double)currentRow/(double)rows;
      for (col = 0; col < cols; col++) {
	cx = startx + (endx-startx)*(double)col/(double)cols;
	color = mandelbrot(cx,cy,its);
	if(color == its) { color = 0; }
	myRow[col] = color;
      }
      MPI_Send(&currentRow, 1, MPI_INT, 0, ROW_DONE, MPI_COMM_WORLD);
      MPI_Send(myRow, cols, MPI_INT, 0, SENDING_ROW, MPI_COMM_WORLD);
      MPI_Recv(&message, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      //printf("processors %d received %d\n", myRank, message);
      currentRow = message;
    }
  }



  endTime = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  for(i = 0; i < P; i++) {
    if(myRank == i) {
      printf("Processor %d took %f seconds\n", myRank, (endTime - startTime));
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if(myRank == 0) {
    //printf("Processor %d took %f seconds\n", myRank, (endTime - startTime));
    FILE *outFile = fopen("outFile", "w");
    fwrite(&rows, sizeof(int), 1, outFile);
    fwrite(&cols, sizeof(int), 1, outFile);
    fwrite(&its, sizeof(int), 1, outFile);
    fwrite(gatheredPicture, sizeof(int), rows*cols, outFile);
    fclose(outFile);
    free(gatheredPicture);
  }
  free(myRow);
  MPI_Finalize();
  return 0;
}

int mandelbrot(double cx, double cy, int max_its) {

  double zr = 0;
  double zi = 0;
  int iterations = 0;
  double value = 0;
  double temp;

  while(value <= 4 && iterations < max_its) {
    temp = zr; zr = zr*zr -(zi *zi) + cx;
    zi = 2*temp*zi + cy;
    value = zr*zr+zi*zi;
    iterations++;
  }

  return iterations;
}
