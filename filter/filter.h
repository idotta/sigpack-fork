// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_FILTER_H
#define SP_FILTER_H
namespace sp
{
    ///
    /// @defgroup filter Filter
    /// \brief FIR/MA and IIR/ARMA filter functions.
    /// @{

    ///
    /// \brief FIR/MA filter class.
    ///
    /// Implements FIR/MA filter functions as \f[  y(n) = \sum_{k=0}^{N-1}{b_kx(n-k)}=b_0x(n)+b_1x(n-1)+...+b_{N-1}x(n-(N-1))\f]
    /// where N is the number of taps in the FIR filter. The filter order is N-1.
    /// Adaptive update of filter is possible with LMS or NLMS algorithms
    template <class T1, class T2, class T3>
    class FIR_filt
    {
    private:
        // Ordinary FIR filter
        int N;                ///< Nr of filter taps
        int cur_p;            ///< Pointer to current sample in buffer
        arma::Mat<T1> buf;    ///< Signal buffer
        arma::Mat<T2> b;      ///< Filter coefficients
        // Adaptive LMS FIR filter
        double mu;            ///< Adaptive filter step size
        int L;                ///< Adaptive filter block length
        int blk_ctr;          ///< Adaptive filter block length counter
        T2 c;                 ///< Adaptive filter NLMS regulation const.
        arma::Mat<T1> P;      ///< Adaptive filter Inverse corr matrix
        arma::Mat<T1> K;      ///< Adaptive filter gain vector
        double lmd;           ///< Adaptive filter RLS forgetting factor
        int do_adapt;         ///< Adaptive filter enable flag
    public:
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Constructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        FIR_filt(){}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Destructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        ~FIR_filt(){}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Clears the internal states and pointer.
        ////////////////////////////////////////////////////////////////////////////////////////////
        void clear(void)
        {
            buf.zeros();
            cur_p = 0;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Sets coefficients in FIR filter.
        /// The internal state and pointers are cleared
        /// @param _b Filter coefficients \f$ [b_0 ..b_{N-1}]^T \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void set_coeffs(arma::Mat<T2> &_b)
        {
            N = _b.n_elem;
            buf.set_size(N,1);
            this->clear();
            b.set_size(N,1);
            b = _b;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Sets coefficients in FIR filter (col format)
        /// The internal state and pointers are cleared
        /// @param _b_col Filter coefficients \f$ [b_0 ..b_{N-1}]^T \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void set_coeffs(arma::Col<T2> &_b_col)
        {
          arma::Mat<T2> b_mat = arma::conv_to<arma::Mat<T2> >::from(_b_col);
          set_coeffs(b_mat);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Get coefficients from FIR filter.
        /// @return b Filter coefficients \f$ [b_0 ..b_{N-1}]^T \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        arma::Col<T2> get_coeffs()
        {
           return arma::conv_to<arma::Col<T2> >::from(b);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Updates coefficients in FIR filter without clearing the internal states.
        /// @param _b Filter coefficients \f$ [b_0 ..b_{N-1}] \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void update_coeffs(arma::Mat<T2> &_b)
        {
            b = _b;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter operator.
        /// @return Filtered output
        /// @param in Input sample
        ////////////////////////////////////////////////////////////////////////////////////////////
        T3 operator()(T1 & in)
        {
            T3 out=0;
            int p = 0;
            buf[cur_p] = in;                    // Insert new sample
            for( int n = cur_p; n < N; n++)
                out += b[p++]*buf[n];           // Calc upper part
            for( int n = 0; n < cur_p; n++)
                out += b[p++]*buf[n];           // ... and lower

            cur_p--;                            // Move insertion point
            if (cur_p < 0) cur_p = N-1;

            return out;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter function.
        /// @return Filtered output
        /// @param in Input vector
        ////////////////////////////////////////////////////////////////////////////////////////////
        arma::Mat<T3> filter(arma::Mat<T1> & in)
        {
            long int sz = in.n_elem;
            arma::Mat<T3> out(sz,1);
            for(long int n=0;n<sz;n++)
                out[n] = this->operator()(in[n]);
            return out;
        }
        arma::Col<T3> filter(arma::Col<T1> & in)
        {
           arma::Mat<T1> in_col = arma::conv_to<arma::Mat<T1> >::from(in);
           return arma::conv_to<arma::Col<T3> >::from(filter(in_col));
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief LMS Filter function setup.
        /// @param _N  Number of filter taps
        /// @param _mu Step size
        /// @param _L Block length
        ////////////////////////////////////////////////////////////////////////////////////////////
        void setup_lms(const unsigned int _N, const double _mu, const unsigned int _L=1)
        {
            N  = _N;
            mu = _mu;
            L  = _L;
            buf.set_size(N,1);buf.zeros();
            b.set_size(N,1);b.zeros();
            K.set_size(N,1);K.zeros();
            cur_p = 0;
            blk_ctr = 0;
            do_adapt = 1;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief LMS Filter update function
        ///
        ///  The LMS filter is updated as <br>
        ///  \f$ \mathbf{b(n)} = \mathbf{b(n-1)}+2\mu\mathbf{x(n)}err(n) \f$ <br>
        ///  where <br>\f$ err(n) = d(n)-\mathbf{b(n-1)^Tx(n)} \f$
        /// @param _err  Feedback error
        ////////////////////////////////////////////////////////////////////////////////////////////
        void lms_adapt(T3 _err)
        {
            if(do_adapt)
            {
                // Reshape buf
                arma::Mat<T1> buf_tmp(N,1);
                for(int k=0; k<N; k++)
                {
                    buf_tmp(k) = buf((cur_p+k+1)%N);
                }

                // Accumulate
                K += _err*buf_tmp;

                // Block update
                if(blk_ctr++%L==0)
                {
                      b+=2*mu*K/L;
                      K.zeros();
                }
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief NLMS Filter function setup.
        /// @param _N  Number of filter taps
        /// @param _mu Step size
        /// @param _c Regularization factor
        /// @param _L Block length
        ////////////////////////////////////////////////////////////////////////////////////////////
        void setup_nlms(const unsigned int _N, const double _mu, const T2 _c, const unsigned int _L=1)
        {
            N  = _N;
            mu = _mu;
            L  = _L;
            c  = _c;
            buf.set_size(N,1);buf.zeros();
            b.set_size(N,1);b.zeros();
            K.set_size(N,1);K.zeros();
            cur_p = 0;
            blk_ctr = 0;
            do_adapt = 1;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief NLMS Filter update function
        ///
        ///  The NLMS filter is updated as <br>
        ///  \f$ \mathbf{b(n)} = \mathbf{b(n-1)}+2\mu\frac{\mathbf{x(n)}err(n)}{c+\mathbf{x(n)^Tx(n)}} \f$ <br>
        ///  where <br>\f$ err(n) = d(n)-\mathbf{b(n-1)^Tx(n)} \f$
        /// @param _err  Feedback error
        ////////////////////////////////////////////////////////////////////////////////////////////
        void nlms_adapt(T3 _err)
        {
            if(do_adapt)
            {
                // Reshape buf
                arma::Mat<T1> buf_tmp(N,1);
                for(int k=0; k<N; k++)
                {
                    buf_tmp(k) = buf((cur_p+k+1)%N);
                }

                // Accumulate
                T1 S = c + arma::as_scalar(buf_tmp.t()*buf_tmp);
                K += _err*buf_tmp/S;

                // Block update
                if(blk_ctr++%L==0)
                {
                      b+=2*mu*K/L;
                      K.zeros();
                }
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief RLS Filter function setup.
        /// @param _N  Number of filter taps
        /// @param _lmd Lambda
        /// @param _P0 Inverse corr matrix initializer
        ////////////////////////////////////////////////////////////////////////////////////////////
        void setup_rls(const unsigned int _N, const double _lmd,const double _P0)
        {
            N  = _N;
            lmd  = _lmd;
            L = 1;
            P.eye(N,N);
            P =_P0*P;
            K.set_size(N,1);K.zeros();
            buf.set_size(N,1);buf.zeros();
            b.set_size(N,1);b.zeros();
            cur_p = 0;
            do_adapt = 1;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief RLS Filter update function
        ///
        ///  The RLS filter is updated as <br>
        ///  \f$ \mathbf{b(n)} = \mathbf{b(n-1)}+\mathbf{Kx(n)}\f$ <br>
        ///  where <br>\f$ \mathbf{K} =\frac{\mathbf{Px}}{\lambda+\mathbf{x^TPx}} \f$ <br>
        ///  and <br>\f$ \mathbf{P^+} =\frac{\mathbf{P^-+xP^-x^T }}{\lambda} \f$ <br>
        /// @param _err  Feedback error
        ////////////////////////////////////////////////////////////////////////////////////////////
        void rls_adapt(T3 _err)
        {
            if(do_adapt)
            {
                // Reshape buf
                arma::Mat<T1> buf_tmp(N,1);
                for(int k=0; k<N; k++)
                {
                    buf_tmp(k) = buf((cur_p+k+1)%N);
                }

                // Update P
                T1 S = lmd + arma::as_scalar(buf_tmp.t()*P*buf_tmp);
                K = P*buf_tmp/S;
                P = (P-K*buf_tmp.t()*P)/lmd;

                // Update coeffs
                b = b + K*_err;
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Get step size
        /// @return Step size mu
        ////////////////////////////////////////////////////////////////////////////////////////////
        double get_step_size(void)
        {
            return mu;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Set step size
        /// @param _mu Step size mu
        ////////////////////////////////////////////////////////////////////////////////////////////
        void set_step_size(const double _mu)
        {
            mu = _mu;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Start adapt
        ////////////////////////////////////////////////////////////////////////////////////////////
        void adapt_enable(void)
        {
            do_adapt = 1;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Stop adapt
        ////////////////////////////////////////////////////////////////////////////////////////////
        void adapt_disble(void)
        {
            do_adapt = 0;
        }

    };


    ///
    /// \brief IIR/ARMA filter class.
    ///
    /// Implements IIR/ARMA filter functions as \f[  a_0y(n) = b_0x(n)+b_1x(n-1)+...+b_{N-1}x(n-(N-1))-a_1y(n-1)-...-a_{M-1}y(n-(M-1))\f]
    /// where N is the number of taps in the FIR filter part and M is the number of taps in the IIR filter. The filter order is (N-1,M-1)
    ///
    template <class T1, class T2, class T3>
    class IIR_filt
    {
    private:
        int N;                ///< Nr of MA filter taps
        int M;                ///< Nr of AR filter taps
        int b_cur_p;          ///< Pointer to current sample in MA buffer
        int a_cur_p;          ///< Pointer to current sample in AR buffer
        arma::Col<T2> b;      ///< MA Filter coefficients
        arma::Col<T2> a;      ///< AR Filter coefficients
        arma::Col<T1> b_buf;  ///< MA Signal buffer
        arma::Col<T1> a_buf;  ///< AR Signal buffer
    public:
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Constructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        IIR_filt(){}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Destructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        ~IIR_filt(){}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Clears the internal states and pointers.
        ////////////////////////////////////////////////////////////////////////////////////////////
        void clear(void)
        {
            b_buf.zeros();
            a_buf.zeros();
            b_cur_p = 0;
            a_cur_p = 0;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Sets coefficients in IIR filter.
        /// The internal state and pointers are cleared
        /// @param _b Filter coefficients \f$ [b_0 ..b_N] \f$
        /// @param _a Filter coefficients \f$ [a_0 ..a_M] \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void set_coeffs(arma::Col<T2> &_b,arma::Col<T2> &_a)
        {
            N = _b.size();
            M = _a.size();
            b_buf.set_size(N);
            a_buf.set_size(M);
            this->clear();
            b = _b/_a[0];
            a = _a/_a[0];
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Updates coefficients in filter without clearing the internal states.
        /// @param _b Filter coefficients \f$ [b_0 ..b_N] \f$
        /// @param _a Filter coefficients \f$ [a_0 ..a_M] \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void update_coeffs(arma::Col<T2> &_b,arma::Col<T2> &_a)
        {
            b = _b/_a[0];
            a = _a/_a[0];
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter operator.
        /// @return Filtered output
        /// @param in Input sample
        ////////////////////////////////////////////////////////////////////////////////////////////
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

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter function.
        /// @return Filtered output
        /// @param in Input vector
        ////////////////////////////////////////////////////////////////////////////////////////////
        arma::Col<T3> filter(arma::Col<T1> & in)
        {
            long int sz = in.size();
            arma::Col<T3> out(sz);
            for(long int n=0;n<sz;n++)
                out[n] = this->operator()(in[n]);
            return out;
        }
    };


    ///
    /// Filter design functions
    ///

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief FIR design functions.
    /// FIR design using windows method (hamming window).
    /// NB! Returns size N+1
    /// @return b Filter coefficients \f$ [b_0 ..b_N] \f$
    /// @param N Filter order
    /// @param f0 Filter cutoff frequency in interval [0..1]
    ////////////////////////////////////////////////////////////////////////////////////////////
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

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Fractional delay function.
    /// Fractional delay filter design using windowed sinc method.
    /// Actual delay is N/2+fd samples for even nr of taps and
    /// (N-1)/2+fd for odd nr of taps
    /// Best performance if -1 < fd < 1
    /// @param N Filter length
    /// @param fd Fractional delay
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec fd_filter( const int N, double fd )
    {
        arma::vec h(N);
        arma::vec w = blackmanharris(N);
        if( N % 2 == 1 ) fd = fd-0.5; // Offset for odd nr of taps
        for(int n=0;n<N;n++)
        {
            h(n) = w(n)*sinc(n-N/2.0-fd);
        }
        h = h/sum(h);  // Normalize gain

        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Frequency response function.
    /// Calculates the frequency response
    /// @param b FIR/MA filter coefficients
    /// @param a IIR/AR filter coefficients
    /// @param M Number of evaluation points, Default 512
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::cx_vec freq( const arma::vec b, const arma::vec a, const int M=512)
    {
        arma::cx_vec h(M);
        int Nb = b.size();
        int Na = a.size();
        std::complex<double> b_tmp,a_tmp,i(0,1);
        for(int m=0;m<M;m++)
        {
            b_tmp=std::complex<double>(b(0),0);
            for(int nb=1;nb<Nb;nb++)
                b_tmp+= b(nb)*(cos(nb*PI*m/M)-i*sin(nb*PI*m/M));
            a_tmp=std::complex<double>(a(0),0);
            for(int na=1;na<Na;na++)
                a_tmp+= a(na)*(cos(na*PI*m/M)-i*sin(na*PI*m/M));
            h(m) = b_tmp/a_tmp;
        }
        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Frequency magnitude response function.
    /// Calculates the frequency magnitude response
    /// @param b FIR/MA filter coefficients
    /// @param a IIR/AR filter coefficients
    /// @param M Number of evaluation points, Default 512
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec freqz( const arma::vec b, const arma::vec a, const int M=512)
    {
        arma::cx_vec f = freq(b,a,M);
        return abs(f);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Frequency phase response function.
    /// Calculates the frequency phase response
    /// @param b FIR/MA filter coefficients
    /// @param a IIR/AR filter coefficients
    /// @param M Number of evaluation points, Default 512
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma::vec phasez( const arma::vec b, const arma::vec a, const int M=512)
    {
        arma::cx_vec f = freq(b,a,M);
        return angle(f);
    }
    /// @}

} // end namespace
#endif
