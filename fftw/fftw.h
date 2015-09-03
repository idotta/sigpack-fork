// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_FFTW_H
#define SP_FFTW_H
#include <fftw3.h>

namespace sp
{
    ///
    /// @defgroup fftw FFTW
    /// \brief One dimensional FFT functions using FFTW3 library.
    ///
    /// \note If a single FFT is to be used the Armadillo version is faster. 
    /// FFTW takes longer time at the first calculation but is faster in the following loops 
    /// @{

    ///
    /// \brief FFTW class.
    ///
    /// Implements FFT functions for Armadillo types. For more info see [fftw.org](http://fftw.org/)
    /// 
	class FFTW
	{
	private:
		fftw_plan p;        ///< FFTW plan 
		int N;              ///< FFT length
		int alg;            ///< One of FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE, FFTW_WISDOM_ONLY see [FFTW plans](http://fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags)
		double Ninv;        ///< Inverse of N
		bool replan;        ///< Replan state
	public:
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Constructor.
        /// @param _N FFT length
        /// @param alg FFTW algorithm selection, Default FFTW_MEASURE
        ////////////////////////////////////////////////////////////////////////////////////////////
		FFTW(int _N, int alg = FFTW_MEASURE)
		{
			p = NULL;
			N = _N;
			Ninv = 1.0/N;
			replan = true;
		}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Destructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
		~FFTW()
		{
			fftw_destroy_plan(p);
			fftw_cleanup();
		}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief FFT of complex input.
        /// @param x Complex input data
        /// @param[out] Pxx Vector to hold complex FFT of length N 
        ////////////////////////////////////////////////////////////////////////////////////////////
		void fft_cx( arma::cx_vec &x, arma::cx_vec &Pxx)
		{
			fftw_complex*  in = reinterpret_cast<fftw_complex*>(x.memptr());
			fftw_complex* out = reinterpret_cast<fftw_complex*>(Pxx.memptr());
			if (replan)
			{
				p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, alg);
				replan = false;
			}
			fftw_execute_dft(p, in, out);
		}
        
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief FFT of complex input.
        /// @returns Complex FFT of length N 
        /// @param x Complex input data
        ////////////////////////////////////////////////////////////////////////////////////////////
		arma::cx_vec fft_cx(arma::cx_vec &x)
		{
		   arma::cx_vec Pxx(N);
		   fft_cx(x, Pxx);
		   return Pxx;
		}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Inverse FFT.
        /// @param Pxx Complex FFT
        /// @param[out] x Vector to hold complex data of length N 
        ////////////////////////////////////////////////////////////////////////////////////////////
		void ifft_cx(arma::cx_vec &Pxx, arma::cx_vec &x)
		{
			fftw_complex*  in = reinterpret_cast<fftw_complex*>(Pxx.memptr());
			fftw_complex* out = reinterpret_cast<fftw_complex*>(x.memptr());
			if (replan)
			{
				p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, alg);
				replan = false;
			}
			fftw_execute_dft(p, in, out);
			x *= Ninv;
		}
        
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Inverse FFT.
        /// @returns Complex data vector of length N 
        /// @param Pxx Complex FFT
        ////////////////////////////////////////////////////////////////////////////////////////////
		arma::cx_vec ifft_cx(arma::cx_vec &Pxx)
		{
			arma::cx_vec x(N);
			ifft_cx(Pxx, x);
			return x;
		}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief FFT of real input.
        /// @param x Input data
        /// @param[out] Pxx Vector to hold complex FFT of length N 
        ////////////////////////////////////////////////////////////////////////////////////////////
		void fft( arma::vec &x, arma::cx_vec &Pxx)
		{
			double*        in = x.memptr();
			fftw_complex* out = reinterpret_cast<fftw_complex*>(Pxx.memptr());
			if (replan)
			{
				p = fftw_plan_dft_r2c_1d(N, in, out, alg);
				replan = false;
			}
			fftw_execute_dft_r2c(p, in, out);
			int offset = ceil(N/2.0);
			int n_elem = N - offset;
			for (int i = 0; i < n_elem; ++i) {
			     Pxx(offset + i) = std::conj(Pxx(n_elem - i));
			
			}
		}
        
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief FFT of real input.
        /// @returns Complex FFT of length N 
        /// @param x Real input data
		arma::cx_vec fft(arma::vec &x)
		{
			arma::cx_vec Pxx(N);
			fft(x, Pxx);
			return Pxx;
		}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Inverse FFT.
        /// @param Pxx Complex FFT
        /// @param[out] x Vector to hold real data of length N 
        ////////////////////////////////////////////////////////////////////////////////////////////
		void ifft(arma::cx_vec &Pxx, arma::vec &x)
		{
			fftw_complex* in = reinterpret_cast<fftw_complex*>(Pxx.memptr());
			double*      out = x.memptr();
			if (replan)
			{
				p = fftw_plan_dft_c2r_1d(N, in, out, alg | FFTW_PRESERVE_INPUT);
				replan = false;
			}
			fftw_execute_dft_c2r(p, in, out);
			x *= Ninv;
		}
        
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Inverse FFT.
        /// @returns Real data vector of length N 
        /// @param Pxx Complex FFT
        ////////////////////////////////////////////////////////////////////////////////////////////
		arma::vec ifft(arma::cx_vec &Pxx)
		{
			arma::vec x(N);
			ifft(Pxx, x);
			return x;
		}
	};
    /// @}
   
} // end namespace
#endif