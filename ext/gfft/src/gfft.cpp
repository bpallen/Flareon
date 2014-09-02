
#include <GFFT/GFFT.h>

using namespace gfft;
using namespace std;

typedef double Tp;


static int test()
{
	// range of the needed GFFT classes
	const unsigned Min = 1;
	const unsigned Max = 27;

	// initialization of the object factory
	Loki::Factory<AbstractFFT<Tp>,unsigned int> gfft_factory;
	FactoryInit<GFFTList<GFFT,Min,Max>::Result>::apply(gfft_factory);

	unsigned long i;

	// runtime definition of the data length
	int p = 2;

	// create an instance of the GFFT
	AbstractFFT<Tp>* gfft = gfft_factory.CreateObject(p);

	// sample data
	unsigned long n = 1<<p;
	Tp* data = new Tp [2*n];
	for (i=0; i<n; ++i) {
		data[2*i] = 2*i;
		data[2*i+1] = 2*i+1;
	}

	// compute FFT
	gfft->fft(data);

	// print the results
	cout<<"--------------------------------" << endl;
	  for (i=0; i<n; ++i)
		cout << setw(10) << setprecision(5) << data[2*i] << "\t"
			 << data[2*i+1] << "I" << endl;

	delete [] data;

	return 0;
}
