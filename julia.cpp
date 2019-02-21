#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h> 
#include <vector>
#include <time.h>
#include <mpi.h>
#include <cstdio>
#include <cstring>
#include <complex>

using namespace std;

const int domainHeight=1024; 
const int domainWidth=1024; 
const int maxIteration=100; 

const double yMin= -0.135;   // Mandelbrot scene y - range 
const double yMax= -0.115;
const double xMin= -0.79;  // Mandelbrot scene x -range
const double xMax= -0.777;

const double cx = 0.285;
const double cy = 0.013;

int image[domainHeight][domainWidth] = {0};  

void WriteTGA_RGB(const char* filename, unsigned char* data, unsigned int width, unsigned int height)
{
	FILE *f = fopen(filename, "wb");
	if (!f) {
		fprintf(stderr, "Unable to create output TGA image `%s'\n", filename);
		exit(EXIT_FAILURE);
	}

	fputc(0x00, f); /* ID Length, 0 => No ID        */
	fputc(0x00, f); /* Color Map Type, 0 => No color map included   */
	fputc(0x02, f); /* Image Type, 2 => Uncompressed, True-color Image */
	fputc(0x00, f); /* Next five bytes are about the color map entries */
	fputc(0x00, f); /* 2 bytes Index, 2 bytes length, 1 byte size */
	fputc(0x00, f);
	fputc(0x00, f);
	fputc(0x00, f);
	fputc(0x00, f); /* X-origin of Image    */
	fputc(0x00, f);
	fputc(0x00, f); /* Y-origin of Image    */
	fputc(0x00, f);
	fputc(width & 0xff, f); /* Image Width      */
	fputc((width >> 8) & 0xff, f);
	fputc(height & 0xff, f); /* Image Height     */
	fputc((height >> 8) & 0xff, f);
	fputc(0x18, f); /* Pixel Depth, 0x18 => 24 Bits */
	fputc(0x20, f); /* Image Descriptor     */

	for (int y = height - 1; y >= 0; y--) {
		for (size_t x = 0; x < width; x++) {
			const size_t i = (y * width + x) * 3;
			fputc(data[i + 2], f); /* write blue */
			fputc(data[i + 1], f); /* write green */
			fputc(data[i], f); /* write red */
		}
	}
}

int converges(double zx,double zy)
{
	int n=0;

	while( n<maxIteration )
	{
		double new_zx = zx*zx - zy*zy + cx;
		zy = 2*zx*zy + cy;
		zx = new_zx;
		n++;
		if( (zx*zx + zy*zy) > 1.0f ){
			return 255;
		}
	}
	return 0;
}

int main(int argc, char **argv)
{
	int id, nproc;
	MPI_Status status;
	int id_from;

	double resX = 0;  
	double resY = 0; 
	double zx = 0;
	double zy = 0;

	double s_time, e_time; 


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Barrier(MPI_COMM_WORLD); 

  	if(id == 0){ 

    	s_time=MPI_Wtime(); 

    	for(int j=1;j<nproc;++j){
    	  	MPI_Recv(&id_from, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);

	    	for(int i=id_from-1;i<domainHeight;i += (nproc-1))
				MPI_Recv(&image[i][0], domainWidth, MPI_INT, id_from, 2, MPI_COMM_WORLD, &status);
    	}

    	e_time=MPI_Wtime();
    
     	cout << "Time elapsed during calculation: " << e_time-s_time << " secs." << endl;;

    	unsigned char *data = new unsigned char[domainHeight * domainWidth * 3];
		std::memset(data, 0, domainHeight * domainWidth * 3 * sizeof(unsigned char));

     	for(int i=0;i<domainHeight;++i)
    	{
			for(int j=0;j<domainWidth;++j)
			{
				if(image[i][j] == 255)
				{
					data[(i + j * domainWidth) * 3 + 0] = 255;
					data[(i + j * domainWidth) * 3 + 1] = 255;
					data[(i + j * domainWidth) * 3 + 2] = 255;
				}	  			
 		  	}
       	}
     	WriteTGA_RGB("julia.tga", data, domainWidth, domainHeight);

   		e_time=MPI_Wtime(); 
   
   		cout << "Time elapsed total: " << e_time-s_time << " secs \r\n";
  	}
  	else{ 
    	resX = (xMax-xMin) / domainHeight;
    	resY = (yMax-yMin) / domainWidth;

    	zx=zx+resX*(id-1); 

    	for(int i=id-1;i<domainHeight;i += (nproc-1))
      	{
			zy=yMin;  
			for(int j=0;j<domainWidth;j++)
	  		{
	    		image[i][j]=converges(zx,zy);
	    		zy=zy+resY;
	  		}
	
			zx=zx+resX*(nproc-1); 

      	}
    
    	MPI_Send(&id, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    	for(int i=id-1;i<domainHeight;i += (nproc-1))
      	{
			MPI_Send(&image[i][0], domainWidth, MPI_INT, 0, 2, MPI_COMM_WORLD);
      	}
    
  	}
  
	MPI_Finalize();
  	return 0;
}