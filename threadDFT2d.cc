// Threaded two-dimensional Discrete FFT transform
// Jake Ashmore
// ECE4122 Project 2


#include <iostream>
#include <string>
#include <math.h>
#include <pthread.h>

#include "Complex.h"
#include "InputImage.h"

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.
Complex* imageData;
int width, height;
int numThreads=16;		//creating 16 threads in addition to main
int N=1024;					//Size of our image
Complex* W;
pthread_barrier_t barrier;
pthread_barrierattr_t attr;
Complex* temp_arr;
int rc;

using namespace std;

void makeW() {
	for(int i=0;i<N/2;i++) {
		W[i]=Complex(cos((2*M_PI*i)/N),-sin((2*M_PI*i)/N));
	}
}

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}
   
void flipthembits(Complex* in_arr, Complex* out_arr) {
	for(int i=0;i<N;i++) {
		int start = i*N;						//looks to the start of each row
		Complex* temp = new Complex[N]; 		//creates a temporary array of size 1024
		for(int j=0;j<N;j++) {
			temp[j]=in_arr[j+start];
		}
		for(int j=start;j< start+N;j++) {
			out_arr[j]=temp[ReverseBits(j)];	//reverses the bits for each row
		}
		delete [] temp;
	}
	for(int j=0;j<N*N;j++) {
		in_arr[j]=out_arr[j]; 					//allows us to keep transforming in place
	}
}

void transpose(Complex* before_t, Complex* after_t) {
	for(int j=0; j<N; j++) {
		for(int i=0;i<N;i++) {
			after_t[i*width+j]=before_t[j*height+i];
		}
	}
	for(int k=0; k<N*N; k++) {
		before_t[k]=after_t[k];			//allows a transform in place
	}
}

void Transform1D(Complex* h, int N)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
	Complex temp;
	for(int n=2; n <= N; n=n*2) {
		for(int i=0; i<N; i=i+n) {
			for(int j=0; j<n/2; j++) {
				int off = n/2;
				temp = h[i+j];
				h[i+j]=h[i+j]+W[j*N/n]*h[i+j+off];
				h[i+j+off]=temp-W[j*N/n]*h[i+j+off];
			}
		}
	}
}

void* Transform2DTHread(void* v)
{ // This is the thread startign point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
	
	unsigned long myId = (unsigned long)v;
	int startRow = myId*(N/numThreads);
	for(int i=0; i<N/numThreads; i++) {
		int startElement = N*(startRow+i);
		Transform1D(imageData + startElement, N);			//thread 0 starts at 0, thread 1 starts at 65536... etc
	}
  // wait for all to complete
	//cout << "waiting" << myId << endl;
	pthread_barrier_wait(&barrier);  
	
  return 0;
}

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
	
  InputImage image(inputFN);  // Create the helper object for reading the image
  imageData = image.GetImageData();
  width = image.GetWidth();
  height = image.GetHeight();
  
  rc = pthread_barrier_init(&barrier, &attr, numThreads+1);			//initialize the barrier function
  if(rc!=0){
	  cout << "error initializing pthread barrier" << endl;
  }
  
  // Create the global pointer to the image array data
  // Create 16 threads
  
  W = new Complex[N/2];
  temp_arr = new Complex[N*N];
  
  flipthembits(imageData, temp_arr);
  makeW();
  
  for(int i=0; i<numThreads;++i) {
	  pthread_t pt;
	  pthread_create(&pt, 0, Transform2DTHread, (void *)i );
	  //cout << "making threads" << endl;
  }
  //cout << "right before I make an image... main is waiting" << endl;
  pthread_barrier_wait(&barrier);
    
  //cout<<"I am making you an image of after1D"<<endl;
  image.SaveImageData("MyAfter1D.txt", imageData, width, height);
  
  int rc = pthread_barrier_init(&barrier, &attr, numThreads+1);			//initialize the barrier function
  if(rc!=0){
	cout << "error initializing pthread barrier" << endl;
  }
  
  transpose(imageData, temp_arr);
  flipthembits(imageData, temp_arr);
  
  for(int i=0; i<numThreads;++i) {
  	  pthread_t pt;
  	  pthread_create(&pt, 0, Transform2DTHread, (void *)i );
  	  //cout << "making threads" << endl;
    }
  
  //cout << "right before I make an image... main is waiting" << endl;
  pthread_barrier_wait(&barrier);
  transpose(imageData, temp_arr);
  
  //cout<<"I am amking you another image after2D"<<endl;
  image.SaveImageData("MyAfter2D.txt", imageData, width, height);
  
  // Wait for all threads complete
  // Write the transformed data
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line

  Transform2D(fn.c_str()); // Perform the transform.
}  
  

  
