// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SP_SPECTRUM_H
#define SP_SPECTRUM_H
namespace sp
{
    //////////////////////////////////////////////////////////////////
    // Spectrum functions
    //////////////////////////////////////////////////////////////////

    ///////////////////////////////////
    // Pxx = spectrum(x,W)
    //      Spectrum of x using window W
    template <class T1>
    arma::cx_vec spectrum(arma::Col<T1> &x, arma::vec &W)
    {
        arma::cx_vec Pxx(x.size());
        double wc = sum(W);     // Window correction factor
        Pxx = fft(x % W)/wc;    // FFT calc
        return Pxx;
    }

    ///////////////////////////////////
    // Pxx = psd(x,W)
    //      Power spectrum density of x using window W
    template <class T1>
    arma::vec psd(arma::Col<T1> &x, arma::vec &W)
    {
        arma::cx_vec X(x.size());
        arma::vec Pxx(x.size());
        X = spectrum(x,W);          // FFT calc
        Pxx = real(X % conj(X));    // Calc power spectra and compensate
        return Pxx;
    }

    ///////////////////////////////////
    // Pxx = psd(x)
    //      Power spectrum density of x using Hamming window
    template <class T1>
    arma::vec psd(arma::Col<T1> &x)
    {
        arma::vec W;
        W = hamming(x.size());
        return psd(x,W);
    }

    ////////-///////////////////////////
    // Pxx = specgram_cx(x,Nfft,Noverl)
    //      Spectrogram of signal
    //      x      - input signal
    //      Nfft   - FFT size
    //      Noverl - Nr of samples overlap
    template <class T1>
    arma::cx_mat specgram_cx(arma::Col<T1> &x, const int Nfft=512, const int Noverl=256)
    {
        arma::cx_mat Pw;

        //Def params
        int N = x.size();
        int D = Nfft-Noverl;
        int m = 0;
        if(N > Nfft)
        {
            arma::Col<T1> xk(Nfft);
            arma::vec W(Nfft);

            W = hamming(Nfft);
            int U = int(floor((N-Noverl)/double(D)));
            Pw.set_size(Nfft,U);
            Pw.zeros();

            // Avg loop
            for(int k=0;k<N-Nfft;k+=D)
            {
                xk = x.rows(k,k+Nfft-1);       // Pick out chunk
                Pw.col(m++) = spectrum(xk,W);  // Calculate spectrum
            }
        }
        else
        {
			arma::vec W(N);
            W = hamming(N);
            Pw.set_size(N,1);
            Pw = spectrum(x,W);
        }
        return Pw;
    }

    ///////////////////////////////////
    // Pxx = specgram(x,Nfft,Noverl)
    //      Power spectrogram of signal
    //      x      - input signal
    //      Nfft   - FFT size
    //      Noverl - Nr of samples overlap
    template <class T1>
    arma::mat specgram(arma::Col<T1> &x, const int Nfft=512, const int Noverl=256)
    {
        arma::cx_mat Pw;
        arma::mat Sg;
        Pw = specgram_cx(x,Nfft,Noverl);
        Sg = real(Pw % conj(Pw));              // Calculate power spectrum
        return Sg;
    }
    ///////////////////////////////////
    // Pxx = specgram_ph(x,Nfft,Noverl)
    //      Phase spectrogram of signal
    //      x      - input signal
    //      Nfft   - FFT size
    //      Noverl - Nr of samples overlap
    template <class T1>
    arma::mat specgram_ph(arma::Col<T1> &x, const int Nfft=512, const int Noverl=256)
    {
        arma::cx_mat Pw;
        arma::mat Sg;
        Pw = specgram_cx(x,Nfft,Noverl);
        Sg = angle(Pw);                        // Calculate phase spectrum
        return Sg;
    }

    ///////////////////////////////////
    // Pxx = pwelch_ph(x,Nfft,Noverl)
    //      Phase spectrum estimation using Welch method
    template <class T1>
    arma::vec pwelch_ph(arma::Col<T1> &x, const int Nfft=512, const int Noverl=256)
    {
        arma::mat Ph;
        Ph  = specgram_ph(x,Nfft,Noverl);
        return arma::mean(Ph ,1);
    }

    ///////////////////////////////////
    // Pxx = pwelch(x,Nfft,Noverl)
    //      Spectrum estimation using Welch method
    //      abs(pwelch(x,Nfft,Noverl)) is equvivalent
    //      to Matlab's: pwelch(x,Nfft,Noverl,'twosided','power')
    template <class T1>
    arma::vec pwelch(arma::Col<T1> &x, const int Nfft=512, const int Noverl=256)
    {
        arma::mat Pxx;
        Pxx = specgram(x,Nfft,Noverl);
        return arma::mean(Pxx,1);
    }
} // end namepace
#endif
