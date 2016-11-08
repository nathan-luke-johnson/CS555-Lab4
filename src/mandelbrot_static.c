/* Mandelbrot set program.
 */

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<time.h>
#include<string.h>
int mandelbrot(double cx, double cy, int max_its);

int main(int argc, char * argv[]) {

  int its = 255;
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
  int *gatheredPicture = NULL;
  double startTime;
  double endTime;
  if(myRank == 0)
  {
    gatheredPicture = malloc(sizeof(int)*rows*cols);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  startTime = MPI_Wtime();
  
  int *picture = malloc(sizeof(int)*rows*cols);
  int row; int col;
  double cx; double cy;
  int color;

  for (row = (rows/P)*myRank; row < (rows/P)*(myRank+1); row++) {
      cy = starty + (endy-starty)*(double)row/(double)rows;
    for(col = 0; col < cols; col++) {
      cx = startx + (endx-startx)*(double)col/(double)cols;
      color = mandelbrot(cx,cy,its);
      if (color == its) { color = 0; }
      picture[cols*row + col] = color;
    }
  }

  
  MPI_Gather(&picture[cols*(rows/P)*myRank], cols*(rows/P), MPI_INT,
      gatheredPicture, cols*(rows/P), MPI_INT,
      0, MPI_COMM_WORLD);
  endTime = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  for(i = 0; i < P; i++) {
    if(myRank == i) {
      printf("Processor %d took %f seconds\n", myRank, (endTime - startTime));
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if(myRank == 0) {
    FILE *outFile = fopen("outFile", "w");
    fwrite(&rows, sizeof(int), 1, outFile);
    fwrite(&cols, sizeof(int), 1, outFile);
    fwrite(&its, sizeof(int), 1, outFile);
    fwrite(gatheredPicture, sizeof(int), rows*cols, outFile);
    fclose(outFile);
    free(gatheredPicture);
  }
  free(picture);
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
