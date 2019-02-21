// Fractal.cpp : Defines the entry point for the console application.
//

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <complex>
#include <omp.h>
#include <iostream>

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

int main()
{
	const unsigned int domainWidth = 1024;
	const unsigned int domainHeight = 1024;
	unsigned char *data = new unsigned char[domainWidth * domainHeight * 3];
	std::memset(data, 0, domainWidth * domainHeight * 3 * sizeof(unsigned char));

	std::complex<double> c(0.285, 0.013);
	std::complex<double> K(0.353, 0.288);
	std::complex<double> center(-1.23, -1.23);
	double scale = 2.35;

	double startTime, endTime;

    std::cout.precision(12);


	const unsigned int maxIterations = 100;
	startTime = omp_get_wtime();

	for (unsigned int y = 0; y < domainHeight; ++y)
	{
		for (unsigned int x = 0; x < domainWidth; ++x)
		{
			std::complex<double> z(x / (double)domainWidth * scale + center.real(),
				y / (double)domainHeight * scale + center.imag());

			//std::complex<double> z(c);
			for (unsigned int iteration = 0; iteration < maxIterations; ++iteration)
			{
				z = z * z + c;
				if (std::abs(z) > 1.0f)
				{
					data[(x + y * domainWidth) * 3 + 0] = 255;
					data[(x + y * domainWidth) * 3 + 1] = 255;
					data[(x + y * domainWidth) * 3 + 2] = 255;
				}
			}
		}
	}
	endTime = omp_get_wtime();
    std::cout << "Computation time: " << endTime - startTime << std::endl;
    std::cout << "Precision: " << omp_get_wtick() << std::endl;


	WriteTGA_RGB("julia.tga", data, domainWidth, domainHeight);
	delete[] data;
	return 0;
}
