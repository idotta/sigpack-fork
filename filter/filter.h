// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_FILTER_H
#define SP_FILTER_H
namespace sp
{
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
        //      Clears FIR filter internal state
        void clear(void)
        {
            buf.zeros();
            cur_p = 0;
        }

        ///////////////////////////////////
        // set_coeffs(b)
        //      Sets coefficients in FIR filter
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
        //      FIR filter operator, sample input.
        T3 operator()(T1 & in)
        {
            T3 out=0;
            int p = 0;
            buf[cur_p] = in;                    // Insert new sample
            for(  int n = cur_p; n < N; n++)
                out += b[p++]*buf[n];           // Calc upper part
            for(int n = 0; n < cur_p; n++)
                out += b[p++]*buf[n];           // ... and lower

            cur_p--;                            // Move insertion point
            if (cur_p < 0) cur_p = N-1;

            return out;
        }

        ///////////////////////////////////
        // T3 = filter(T1)
        //      FIR filter function, vector version.
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
        //      Clears IIR filter internal state
        void clear(void)
        {
            b_buf.zeros();
            a_buf.zeros();
            b_cur_p = 0;
            a_cur_p = 0;
        }

        ///////////////////////////////////
        // set_coeffs(b,a)
        //      Sets coefficients in FIR filter
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
        //      IIR filter operator, sample input.
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
        //      IIR filter function, vector version.
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
    //      FIR design using windows method.
    //      NB! Returns size N+1
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

    ///////////////////////////////////
    // vec = fd_filter(N,fd)
    //      Fractional delay filter design using windowed sinc method.
    //      Actual delay is N/2+fd samples for even nr of taps and
    //      (N-1)/2+fd for odd nr of taps
    //      Best performance if -1 < fd < 1
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
} // end namespace
#endif