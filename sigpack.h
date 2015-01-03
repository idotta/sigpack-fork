// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Version
//  1.01	Claes Rolén		2015-01-03	First version

#ifndef ARMA_INCLUDES
#include <armadillo>       // Armadillo version 4.320.2 is used
#endif

namespace sp
{
	//////////////////////////////////////////////////////////////////
	// Math functions
	//////////////////////////////////////////////////////////////////
	#define PI   3.14159265358979323846     // ... or use arma::datum::pi
	#define PI_2 6.28318530717958647692

	///////////////////////////////////
	// y = sinc(x)
	//		Sinc function
	double sinc( double x )
	{
		if(x==0.0)
			return 1.0;
		else
			return std::sin(PI*x)/(PI*x);
	}

	///////////////////////////////////
	// y = besseli0()
	//		Modified first kind bessel function order zero
	double besseli0( double x )
	{
		double y=1.0,s=1.0,x2=x*x;
		int n = 1;
		while (s > y*1.0e-9) 
		{
			s *= x2/4.0/(n*n);
			y += s;
			n++;
		}
		return y;
	}

	//////////////////////////////////////////////////////////////////
	// Windows
	//////////////////////////////////////////////////////////////////

	///////////////////////////////////
	// b = cos_win(N,a)
	//		Generic fifth order symmetric cos window
	arma::vec cos_win( int N, arma::vec &a )
	{
		arma::vec h(N);
		for(int i=0;i<N;i++)
		{
			h[i] = a[0] - a[1]*std::cos(1.0*PI_2*i/(N-1)) + a[2]*std::cos(2.0*PI_2*i/(N-1)) \
			 	        - a[3]*std::cos(3.0*PI_2*i/(N-1)) + a[4]*std::cos(4.0*PI_2*i/(N-1));
		}
		return h;
	}

	///////////////////////////////////
	// b = hamming(N)
	arma::vec hamming( int N )
	{
		arma::vec a=arma::zeros<arma::vec>(5);
		a[0] = 0.54;
		a[1] = 0.46;
		return cos_win(N,a);
	}

	///////////////////////////////////
	// b = hann(N)
	arma::vec hann( int N )
	{
		arma::vec a=arma::zeros<arma::vec>(5);
		a[0] = 0.5;
		a[1] = 0.5;
		return cos_win(N,a);
	}

	///////////////////////////////////
	// b = blackman(N)
	arma::vec blackman( int N )
	{
		arma::vec a=arma::zeros<arma::vec>(5);
		a[0] = 0.42; // 7938/18608.0
		a[1] = 0.5;  // 9240/18608.0
		a[2] = 0.08; // 1430/18608.0
		return cos_win(N,a);
	}

	///////////////////////////////////
	// b = blackmanharris(N)
	//      Symmetric BH4 window
	arma::vec blackmanharris( int N )
	{
		arma::vec a=arma::zeros<arma::vec>(5);
		a[0] = 0.35875;
		a[1] = 0.48829;
		a[2] = 0.14128;
		a[3] = 0.01168;
		return cos_win(N,a);
	}

	///////////////////////////////////
	// b = flattopwin(N)
	arma::vec flattopwin( int N )
	{
		arma::vec a=arma::zeros<arma::vec>(5);
		a[0] = 0.21557895;
		a[1] = 0.41663158;
		a[2] = 0.277263158;
		a[3] = 0.083578947;
		a[4] = 0.006947368;
		return cos_win(N,a);
	}

	///////////////////////////////////
	// b = hanning(N)
	arma::vec hanning( int N )
	{
		arma::vec h(N);
		for(int i=0;i<N;i++)
		{
			h[i] = 0.5-0.5*std::cos(PI_2*(i+1)/(N+1));
		}
		return h;
	}

	///////////////////////////////////
	// b = kaiser(N,beta)
	arma::vec kaiser( int N, double beta )
	{
		arma::vec h(N);
		double bb = besseli0(beta);
		for(int i=0;i<N;i++)
		{
			h[i] = besseli0(beta*sqrt(4.0*i*(N-1-i))/(N-1))/bb;
		}
		return h;
	}

	///////////////////////////////////
	// b = triang(N)
	arma::vec triang( int N )
	{
		arma::vec h(N);
		if(N%2)    // Odd
		{
			for(int i=0;i<(N-1)/2;i++)
			{
				h[i]     = 2.0*(i+1)/(N+1);
				h[N-i-1] = h[i];
			}
			h[(N-1)/2] = 1.0;
		}
		else      // Even
		{
			for(int i=0;i<N/2;i++)
			{
				h[i]     = (2.0*i+1)/N;
				h[N-i-1] = h[i];
			}
		}
		return h;
	}

	//////////////////////////////////////////////////////////////////
	// Delay Class
	//////////////////////////////////////////////////////////////////
	template <class T1>
	class Delay
	{
	private:
		int D;                // Delay
		int cur_p;            // Pointer to current sample in buffer
		arma::Col<T1> buf;    // Signal buffer
	public:
		///////////////////////////////////
		// Constructor
		Delay(){}

		///////////////////////////////////
		// Constructor with delay input
		Delay(int _D)
		{
			set_delay(_D);
			clear();
		}

		///////////////////////////////////
		// clear()
		//		Clears internal state
		void clear(void)
		{
			buf.zeros();
			cur_p = 0;
		}

		///////////////////////////////////
		// set_delay(D)
		//		Sets delay
		void set_delay(int _D)
		{
			D = _D+1;
			buf.set_size(D);
		}

		///////////////////////////////////
		// T1 = xxx(T1)
		//		Delay operator, sample input.
		T1 operator()(T1 & in)
		{
			buf[cur_p] = in;                    // Insert new sample		 
			cur_p--;                            // Move inertion point
			if (cur_p < 0) cur_p = D-1;
			return buf[cur_p];
		}

		///////////////////////////////////
		// T1 = delay(T1)
		//		Delay function, vector version.
		arma::Col<T1> delay(arma::Col<T1> & in)
		{
			long int sz = in.size();
			arma::Col<T1> out(sz);
			for(long int n=0;n<sz;n++)
				out[n] = this->operator()(in[n]);
			return out;
		}
	};

	//////////////////////////////////////////////////////////////////
	// FIR Filter Class
	//////////////////////////////////////////////////////////////////
	template <class T1, class T2, class T3>
	class FIR_filt
	{
	private:
		int N;                // Nr of filter taps
		int cur_p;            // Pointer to current sample in buffer
		arma::Col<T1> buf;    // Signal buffer
		arma::Col<T2> b;      // Filter coefficients
	public:
		///////////////////////////////////
		// Constructor
		FIR_filt(){}

		///////////////////////////////////
		// clear()
		//		Clears FIR filter internal state
		void clear(void)
		{
			buf.zeros();
			cur_p = 0;
		}

		///////////////////////////////////
		// set_coeffs(b)
		//		Sets coefficients in FIR filter
		void set_coeffs(arma::Col<T2> &_b)
		{
			N = _b.size();
			buf.set_size(N);
			buf.zeros();
			b = _b;
			cur_p = 0;
		}

		///////////////////////////////////
		// T3 = xxx(T1)
		//		FIR filter operator, sample input.
		T3 operator()(T1 & in)
		{
			T3 out=0;
			int p = 0;
			buf[cur_p] = in;                    // Insert new sample		 
			for(  int n = cur_p; n < N; n++)      
				out += b[p++]*buf[n];           // Calc upper part
			for(int n = 0; n < cur_p; n++)   
				out += b[p++]*buf[n];           // ... and lower

			cur_p--;                            // Move inertion point
			if (cur_p < 0) cur_p = N-1;

			return out;
		}

		///////////////////////////////////
		// T3 = filter(T1)
		//		FIR filter function, vector version.
		arma::Col<T3> filter(arma::Col<T1> & in)
		{
			long int sz = in.size();
			arma::Col<T3> out(sz);
			for(long int n=0;n<sz;n++)
				out[n] = this->operator()(in[n]);
			return out;
		}
	};


	//////////////////////////////////////////////////////////////////
	// IIR Filter Class
	//////////////////////////////////////////////////////////////////
	template <class T1, class T2, class T3>
	class IIR_filt
	{
	private:
		int N;                // Nr of MA filter taps
		int M;                // Nr of AR filter taps
		int b_cur_p;          // Pointer to current sample in MA buffer
		int a_cur_p;          // Pointer to current sample in AR buffer
		arma::Col<T2> b;      // MA Filter coefficients
		arma::Col<T2> a;      // AR Filter coefficients
		arma::Col<T1> b_buf;  // MA Signal buffer
		arma::Col<T1> a_buf;  // AR Signal buffer
	public:
		///////////////////////////////////
		// Constructor
		IIR_filt(){};

		///////////////////////////////////
		// clear()
		//		Clears IIR filter internal state
		void clear(void)
		{
			b_buf.zeros();
			a_buf.zeros();
			b_cur_p = 0;
			a_cur_p = 0;
		}

		///////////////////////////////////
		// set_coeffs(b,a)
		//		Sets coefficients in FIR filter
		void set_coeffs(arma::Col<T2> &_b,arma::Col<T2> &_a)
		{
			N = _b.size();
			M = _a.size();
			b_buf.set_size(N);
			b_buf.zeros();
			a_buf.set_size(M);
			a_buf.zeros();
			b = _b/_a[0];
			a = _a/_a[0];
			b_cur_p = 0;
			a_cur_p = 0;
		}

		///////////////////////////////////
		// T3 = xxx(T1)
		//		IIR filter operator, sample input.
		T3 operator()(T1 & in)
		{
			T3 out=0;
			int p = 0;

			// MA part
			b_buf[b_cur_p] = in;                // Insert new sample		 
			for(int n = b_cur_p; n < N; n++)      
				out += b[p++]*b_buf[n];         // Calc upper part
			for(int n = 0; n < b_cur_p; n++)   
				out += b[p++]*b_buf[n];         // ... and lower

			b_cur_p--;                          // Move insertion point
			if (b_cur_p < 0) b_cur_p = N-1;

			// AR part
			p=1;
			for(int m = a_cur_p+1; m < M; m++)      
				out -= a[p++]*a_buf[m];         // Calc upper part
			for(int m = 0; m < a_cur_p; m++)   
				out -= a[p++]*a_buf[m];         // ... and lower

			a_buf[a_cur_p] = out;		        // Insert output

			a_cur_p--;
			if (a_cur_p < 0) a_cur_p = M-1;     // Move insertion point

			return out;
		}

		///////////////////////////////////
		// T3 = filter(T1)
		//		IIR filter function, vector version.
		arma::Col<T3> filter(arma::Col<T1> & in)
		{
			long int sz = in.size();
			arma::Col<T3> out(sz);
			for(long int n=0;n<sz;n++)
				out[n] = this->operator()(in[n]);
			return out;
		}
	};

    //////////////////////////////////////////////////////////////////
	// Filter design functions
	//////////////////////////////////////////////////////////////////

	///////////////////////////////////
	// vec = fir1(N,f0)
	//		FIR design using windows method.
	//		NB! Returns size N+1
	arma::vec fir1(int N, double f0)
	{
		arma::vec b(N+1), h(N+1);
		h = hamming(N+1);
		double b_sum=0;
		for (int i=0;i<N+1;i++)
		{
			b[i] = h[i]*sinc(f0*(i-N/2.0));
			b_sum += b[i];
		}
		b = b/b_sum;
		return b;
	}

    //////////////////////////////////////////////////////////////////
	// Resampling functions
	//////////////////////////////////////////////////////////////////

	///////////////////////////////////
	// xN = upsample(x,p)
	//	Upsamples x p times
	template <class T1>
	arma::Col<T1> upsample(arma::Col<T1> &x, const int p )
	{
		long int N = x.size();
		arma::Col<T1> y;
		y.set_size(p*N);
		y.zeros();
		for(long int n=0;n<N;n++)
			y[p*n] = x[n];
		return y;
	}

	///////////////////////////////////
	// xN = upfir(x,p,K)
	//	Upsamples x p times using an K tap anti alias filter 
	template <class T1>
	arma::Col<T1> upfir(arma::Col<T1> &x, const int p, const int K )
	{
		FIR_filt<T1,double,T1> AA;
		AA.set_coeffs(fir1(K,1/float(p)));
		return AA.filter(upsample(x,p));
	}

	///////////////////////////////////
	// xN = downsample(x,q)
	//	Downsamples x q times
	template <class T1>
	arma::Col<T1> downsample(arma::Col<T1> &x, const int q )
	{
		arma::Col<T1> y;
		int N = int(floor(1.0*x.size()/q));
		y.set_size(N);
		for(long int n=0;n<N;n++)
			y[n] = x[n*q];
		return y;
	}

	///////////////////////////////////
	// xN = downfir(x,q,K)
	//	Downsamples x q times using an K tap anti alias filter 
	template <class T1>
	arma::Col<T1> downfir(arma::Col<T1> &x, const int q, const int K )
	{
		FIR_filt<T1,double,T1> AA;
		AA.set_coeffs(fir1(K,1/float(q)));
		return downsample(AA.filter(x),q);
	}

    //////////////////////////////////////////////////////////////////
	// Spectrum functions
	//////////////////////////////////////////////////////////////////

	///////////////////////////////////
	// Pxx = psd(x,W)
	//	 Power spectrum density of x using window W
	template <class T1>
	arma::cx_vec psd(arma::Col<T1> &x, arma::vec &W)
	{
		arma::cx_vec Pxx(x.size());
		double wc = sum(W)*sum(W);     // Window correction factor
		Pxx = fft(x % W);              // FFT calc
		Pxx = Pxx % conj(Pxx);         // Calc power spectra and compensate
		return Pxx/wc;
	}

	///////////////////////////////////
	// Pxx = psd(x)
	//	 Power spectrum density of x using Hamming window
	template <class T1>
	arma::cx_vec psd(arma::Col<T1> &x)
	{
		arma::vec W;
		W = hamming(x.size());
		return psd(x,W);
	}


	///////////////////////////////////
	// Pxx = pwelch(x,Nfft,Noverl)
	//	  Spectrum estimation using Welch method
	//    abs(pwelch(x,Nfft,Noverl)) is equivalent 
	//    to Matlab's: pwelch(x,Nfft,Noverl,'twosided','power')
	template <class T1>
	arma::cx_vec pwelch(arma::Col<T1> &x, const int Nfft=512, const int Noverl=256)
	{
		arma::cx_vec Pw;

		//Def params
		int N = x.size();
		int D = Nfft-Noverl;
		int m = 0;
		if(N > Nfft)
		{
			arma::Col<T1> xk(Nfft);
			arma::vec W(Nfft);

			W = hamming(Nfft);
			Pw.set_size(Nfft);
			Pw.zeros();
			// Avg loop
			for(int k=0;k<N-Nfft;k+=D)
			{
				xk = x.rows(k,k+Nfft-1);  // Pick out chunk
				Pw = Pw + psd(xk,W);      // Accumulate spectra
				m++;
			}    
			Pw = Pw/float(m);             // Take average
		}   
		else
		{    
			Pw.set_size(N);
			Pw = psd(x);     
		}
		return Pw;
	}
} // end namespace - sp
