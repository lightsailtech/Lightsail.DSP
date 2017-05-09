// Code to implement decently performing FFT for complex and real valued                                         
// signals. See www.lomont.org for a derivation of the relevant algorithms                                       
// from first principles. Copyright Chris Lomont 2010-2012.                                                      
// This code and any ports are free for all to use for any reason as long                                        
// as this header is left in place.                                                                              
// Version 1.1, Sept 2011
using System;

namespace Lomont
{
    /// <summary>                                                                                                
    /// Represent a class that performs real or complex valued Fast Fourier                                      
    /// Transforms. Instantiate it and use the FFT or TableFFT methods to                                        
    /// compute complex to complex FFTs. Use FFTReal for real to complex                                         
    /// FFTs which are much faster than standard complex to complex FFTs.                                        
    /// Properties A and B allow selecting various FFT sign and scaling                                          
    /// conventions.                                                                                             
    /// </summary>                                                                                               
    internal class LomontFFT
    {
        /// <summary>                                                                                            
        /// Compute the forward or inverse Fourier Transform of data, with                                       
        /// data containing complex valued data as alternating real and                                          
        /// imaginary parts. The length must be a power of 2. The data is                                        
        /// modified in place.                                                                                   
        /// </summary>                                                                                           
        /// <param name="data">The complex data stored as alternating real                                       
        /// and imaginary parts</param>                                                                          
        /// <param name="forward">true for a forward transform, false for                                        
        /// inverse transform</param>                                                                            
        public void FFT(double[] data, bool forward)
        {
            var n = data.Length;
            // checks n is a power of 2 in 2's complement format                                                 
            if ((n & (n - 1)) != 0)
                throw new ArgumentException(
                    "data length " + n + " in FFT is not a power of 2");
            n /= 2;    // n is the number of samples                                                             

            Reverse(data, n); // bit index data reversal                                                         

            // do transform: so single point transforms, then doubles, etc.                                      
            double sign = forward ? B : -B;
            var mmax = 1;
            while (n > mmax)
            {

                var istep = 2 * mmax;
                var theta = sign * Math.PI / mmax;
                double wr = 1, wi = 0;
                var wpr = Math.Cos(theta);
                var wpi = Math.Sin(theta);
                for (var m = 0; m < istep; m += 2)
                {
                    for (var k = m; k < 2 * n; k += 2 * istep)
                    {
                        var j = k + istep;
                        var tempr = wr * data[j] - wi * data[j + 1];
                        var tempi = wi * data[j] + wr * data[j + 1];
                        data[j] = data[k] - tempr;
                        data[j + 1] = data[k + 1] - tempi;
                        data[k] = data[k] + tempr;
                        data[k + 1] = data[k + 1] + tempi;
                    }
                    var t = wr; // trig recurrence                                                               
                    wr = wr * wpr - wi * wpi;
                    wi = wi * wpr + t * wpi;
                }
                mmax = istep;
            }

            // perform data scaling as needed                                                                    
            Scale(data, n, forward);
        }

        public void NrmlFFT(double[] RlDat, double[] ImDat, bool forward)
        {
            // n is the number of samples                                                             
            var n = RlDat.Length;

            // checks n is a power of 2 in 2's complement format                                                 
            if ((n & (n - 1)) != 0)
                throw new ArgumentException(
                    "data length " + n + " in FFT is not a power of 2");
            if (RlDat.Length != ImDat.Length)
                throw new ArgumentException(
                    "Real data length must match imaginary");

            Reverse(RlDat,ImDat,n); // bit index data reversal                                                         

            // do transform: so single point transforms, then doubles, etc.                                      
            double sign = forward ? B : -B;
            var mmax = 1;
            while (n > mmax)
            {
                var istep = 2 * mmax;
                var theta = sign * Math.PI / mmax;
                double wr = 1, wi = 0;
                var wpr = Math.Cos(theta);
                var wpi = Math.Sin(theta);
                for (var m = 0; m < istep; m += 2)
                {
                    for (var k = m; k < 2 * n; k += 2 * istep)
                    {
                        var j = k + istep;
                        int j_l = j / 2, k_l = k / 2;
                        var tempr = wr * RlDat[j_l] - wi * ImDat[j_l];
                        var tempi = wi * RlDat[j_l] + wr * ImDat[j_l];
                        RlDat[j_l] = RlDat[k_l] - tempr;
                        ImDat[j_l] = ImDat[k_l] - tempi;
                        RlDat[k_l] = RlDat[k_l] + tempr;
                        ImDat[k_l] = ImDat[k_l] + tempi;
                    }
                    var t = wr; // trig recurrence                                                               
                    wr = wr * wpr - wi * wpi;
                    wi = wi * wpr + t * wpi;
                }
                mmax = istep;
            }
            // perform data scaling as needed                                                                    
            Scale(RlDat, n, forward);
            Scale(ImDat, n, forward);
        }

        /// <summary>                                                                                            
        /// Compute the forward or inverse Fourier Transform of data, with data                                  
        /// containing complex valued data as alternating real and imaginary                                     
        /// parts. The length must be a power of 2. This method caches values                                    
        /// and should be slightly faster on than the FFT method for repeated uses.                              
        /// It is also slightly more accurate. Data is transformed in place.                                     
        /// </summary>                                                                                           
        /// <param name="data">The complex data stored as alternating real                                       
        /// and imaginary parts</param>                                                                          
        /// <param name="forward">true for a forward transform, false for                                        
        /// inverse transform</param>                                                                            
        public void TableFFT(double[] data, bool forward)
        {
            var n = data.Length;
            // checks n is a power of 2 in 2's complement format                                                 
            if ((n & (n - 1)) != 0)
                throw new ArgumentException(
                    "data length " + n + " in FFT is not a power of 2"
                    );
            n /= 2;    // n is the number of samples                                                             

            Reverse(data, n); // bit index data reversal                                                         

            // make table if needed                                                                              
            if ((cosTable == null) || (cosTable.Length != n))
                Initialize(n);

            // do transform: so single point transforms, then doubles, etc.                                      
            double sign = forward ? B : -B;
            var mmax = 1;
            var tptr = 0;
            while (n > mmax)
            {
                var istep = 2 * mmax;
                for (var m = 0; m < istep; m += 2)
                {
                    var wr = cosTable[tptr];
                    var wi = sign * sinTable[tptr++];
                    for (var k = m; k < 2 * n; k += 2 * istep)
                    {
                        var j = k + istep;
                        var tempr = wr * data[j] - wi * data[j + 1];
                        var tempi = wi * data[j] + wr * data[j + 1];
                        data[j] = data[k] - tempr;
                        data[j + 1] = data[k + 1] - tempi;
                        data[k] = data[k] + tempr;
                        data[k + 1] = data[k + 1] + tempi;
                    }
                }
                mmax = istep;
            }


            // perform data scaling as needed                                                                    
            Scale(data, n, forward);
        } // end TableFFT

        public void NrmlTableFFT(double[] RlDat, double[] ImDat, bool forward)
        {
            var n = RlDat.Length;
            // checks n is a power of 2 in 2's complement format                                                 
            if ((n & (n - 1)) != 0)
                throw new ArgumentException(
                    "data length " + n + " in FFT is not a power of 2");
            if (RlDat.Length != ImDat.Length)
                throw new ArgumentException(
                    "Real data length must match imaginary");

            // n /= 2;    // n is the number of samples                                                             

            //Reverse(dat, n); // bit index data reversal                                                         
            Reverse(RlDat, ImDat, n);

            // make table if needed                                                                              
            if ((cosTable == null) || (cosTable.Length != n))
                Initialize(n);

            // do transform: so single point transforms, then doubles, etc.                                      
            double sign = forward ? B : -B;
            var mmax = 1;
            var tptr = 0;
            while (n > mmax)
            {
                var istep = 2 * mmax;
                for (var m = 0; m < istep; m += 2)
                {
                    var wr = cosTable[tptr];
                    var wi = sign * sinTable[tptr++];
                    for (var k = m; k < 2 * n; k += 2 * istep)
                    {
                        var j = k + istep;
                        int j_l = j / 2, k_l = k / 2;
                        var tempr = wr * RlDat[j_l] - wi * ImDat[j_l];
                        var tempi = wi * RlDat[j_l] + wr * ImDat[j_l];
                        RlDat[j_l] = RlDat[k_l] - tempr;
                        ImDat[j_l] = ImDat[k_l] - tempi;
                        RlDat[k_l] = RlDat[k_l] + tempr;
                        ImDat[k_l] = ImDat[k_l] + tempi;
                    }
                }
                mmax = istep;
            }
            
            // perform data scaling as needed                                                                    
            Scale(RlDat, n, forward);
            Scale(ImDat, n, forward);
        } // end TableFFT

        /// <summary>                                                                                            
        /// Compute the for-ward or inverse Fourier Trans form of data, with                                       
        /// data containing real-valued data only. The output is complex                                         
        /// valued after the first two entries, stored in alternating real                                       
        /// and imaginary parts. The first two returned entries are the real                                     
        /// parts of the first and last value from the conjugate symmetric                                       
        /// output, which are necessarily real. The length must be a power                                       
        /// of 2.                                                                                                
        /// </summary>                                                                                           
        /// <param name="data">The complex data stored as alternating real                                       
        /// and imaginary parts</param>                                                                          
        /// <param name="forward">true for a forward transform, false for                                        
        /// inverse transform</param>                                                                            
        public void RealFFT(double[] data, bool forward)
        {

            var n = data.Length; // # of real inputs, 1/2 the complex length                                     
            // checks n is a power of 2 in 2's complement format                                                 
            if ((n & (n - 1)) != 0)
                throw new ArgumentException(
                    "data length " + n + " in FFT is not a power of 2"
                    );

            var sign = -1.0; // assume inverse FFT, this controls how algebra below works                        
            if (forward)
            { // do packed FFT. This can be changed to FFT to save memory                                        
                TableFFT(data, true);
                sign = 1.0;
                // scaling - divide by scaling for N/2, then mult by scaling for N                               
                if (A != 1)
                {
                    var scale = Math.Pow(2.0, (A - 1) / 2.0);
                    for (var i = 0; i < data.Length; ++i)
                        data[i] *= scale;
                }
            }

            var theta = B * sign * 2 * Math.PI / n;
            var wpr = Math.Cos(theta);
            var wpi = Math.Sin(theta);
            var wjr = wpr;
            var wji = wpi;

            for (var j = 1; j <= n / 4; ++j)
            {
                var k = n / 2 - j;
                var tkr = data[2 * k];    // real and imaginary parts of t_k  = t_(n/2 - j)                      
                var tki = data[2 * k + 1];
                var tjr = data[2 * j];    // real and imaginary parts of t_j                                     
                var tji = data[2 * j + 1];

                var a = (tjr - tkr) * wji;
                var b = (tji + tki) * wjr;
                var c = (tjr - tkr) * wjr;
                var d = (tji + tki) * wji;
                var e = (tjr + tkr);
                var f = (tji - tki);

                // compute entry y[j]                                                                            
                data[2 * j] = 0.5 * (e + sign * (a + b));
                data[2 * j + 1] = 0.5 * (f + sign * (d - c));

                // compute entry y[k]                                                                            
                data[2 * k] = 0.5 * (e - sign * (b + a));
                data[2 * k + 1] = 0.5 * (sign * (d - c) - f);

                var temp = wjr;
                // todo - allow more accurate version here? make option?                                         
                wjr = wjr * wpr - wji * wpi;
                wji = temp * wpi + wji * wpr;
            }

            if (forward)
            {
                // compute final y0 and y_{N/2}, store in data[0], data[1]                                       
                var temp = data[0];
                data[0] += data[1];
                data[1] = temp - data[1];
            }
            else
            {
                var temp = data[0]; // unpack the y0 and y_{N/2}, then invert FFT                                
                data[0] = 0.5 * (temp + data[1]);
                data[1] = 0.5 * (temp - data[1]);
                // do packed inverse (table based) FFT. This can be changed to regular inverse FFT to save memory
                TableFFT(data, false);
                // scaling - divide by scaling for N, then mult by scaling for N/2                               
                //if (A != -1) // todo - off by factor of 2? this works, but something seems weird               
                {
                    var scale = Math.Pow(2.0, -(A + 1) / 2.0) * 2;
                    for (var i = 0; i < data.Length; ++i)
                        data[i] *= scale;
                }
            }
        } // end RealFFT

        /// <summary>                                                                                            
        /// Determine how scaling works on the forward and inverse transforms.                                   
        /// For size N=2^n transforms, the forward trans form gets divided by                                     
        /// N^((1-a)/2) and the inverse gets divided by N^((1+a)/2). Common                                      
        /// values for (A,B) are                                                                                 
        ///     ( 0, 1)  - default                                                                               
        ///     (-1, 1)  - data processing                                                                       
        ///     ( 1,-1)  - signal processing                                                                     
        /// Usual values for A are 1, 0, or -1                                                                   
        /// </summary>                                                                                           
        public int A { get; set; }

        /// <summary>                                                                                            
        /// Determine how phase works on the forward and inverse transforms.                                     
        /// For size N=2^n transforms, the forward trans form uses an                                             
        /// exp(B*2*pi/N) term and the inverse uses an exp(-B*2*pi/N) term.                                      
        /// Common values for (A,B) are                                                                          
        ///     ( 0, 1)  - default                                                                               
        ///     (-1, 1)  - data processing                                                                       
        ///     ( 1,-1)  - signal processing                                                                     
        /// Abs(B) should be relatively prime to N.                                                              
        /// Setting B=-1 effectively corresponds to conjugating both input and                                   
        /// output data.                                                                                         
        /// Usual values for B are 1 or -1.                                                                      
        /// </summary>                                                                                           
        public int B { get; set; }

        public LomontFFT()
        {
            A = 0;
            B = 1;
        } // end constructor

        #region Internals                                                                                        

        /// <summary>                                                                                            
        /// Scale data using n samples for forward and inverse transforms as needed                              
        /// </summary>                                                                                           
        /// <param name="data"></param>                                                                          
        /// <param name="n"></param>                                                                             
        /// <param name="forward"></param>                                                                       
        void Scale(double[] data, int n, bool forward)
        {
            // forward scaling if needed                                                                         
            if ((forward) && (A != 1))
            {
                var scale = Math.Pow(n, (A - 1) / 2.0);
                for (var i = 0; i < data.Length; ++i)
                    data[i] *= scale;
            }

            // inverse scaling if needed                                                                         
            if ((!forward) && (A != -1))
            {
                var scale = Math.Pow(n, -(A + 1) / 2.0);
                for (var i = 0; i < data.Length; ++i)
                    data[i] *= scale;
            }
        }

        /// <summary>                                                                                            
        /// Call this with the size before using the TableFFT version                                            
        /// Fills in tables for speed. Done automatically in TableFFT                                            
        /// </summary>                                                                                           
        /// <param name="size">The size of the FFT in samples</param>                                            
        void Initialize(int size)
        {
            // NOTE: if you port to non garbage collected languages                                              
            // like C# or Java be sure to free these correctly                                                   
            cosTable = new double[size];
            sinTable = new double[size];

            // forward pass                                                                                      
            var n = size;
            int mmax = 1, pos = 0;
            while (n > mmax)
            {
                var istep = 2 * mmax;
                var theta = Math.PI / mmax;
                double wr = 1, wi = 0;
                var wpi = Math.Sin(theta);
                // compute in a slightly slower yet more accurate manner                                         
                var wpr = Math.Sin(theta / 2);
                wpr = -2 * wpr * wpr;
                for (var m = 0; m < istep; m += 2)
                {
                    cosTable[pos] = wr;
                    sinTable[pos++] = wi;
                    var t = wr;
                    wr = wr * wpr - wi * wpi + wr;
                    wi = wi * wpr + t * wpi + wi;
                }
                mmax = istep;
            }
        }

        /// <summary>                                                                                            
        /// Swap data indices whenever index i has binary                                                        
        /// digits reversed from index j, where data is                                                          
        /// two doubles per index.                                                                               
        /// </summary>                                                                                           
        /// <param name="data"></param>                                                                          
        /// <param name="n"></param>                                                                             
        static void Reverse(double[] data, int n)
        {
            // bit reverse the indices. This is exercise 5 in section                                            
            // 7.2.1.1 of Knuth's TAOCP the idea is a binary counter                                             
            // in k and one with bits reversed in j                                                              
            int j = 0, k = 0; // Knuth R1: initialize                                                            
            var top = n / 2;  // this is Knuth's 2^(n-1)                                                         
            while (true)
            {
                // Knuth R2: swap - swap j+1 and k+2^(n-1), 2 entries each                                       
                var t = data[j + 2];
                data[j + 2] = data[k + n];
                data[k + n] = t;


                t = data[j + 3];
                data[j + 3] = data[k + n + 1];
                data[k + n + 1] = t;

                if (j > k)
                { // swap two more                                     
                                                              
                    // j and k                                                                                   
                    t = data[j];
                    data[j] = data[k];
                    data[k] = t;

                    t = data[j + 1];
                    data[j + 1] = data[k + 1];
                    data[k + 1] = t;
                    
                    // j + top + 1 and k+top + 1                                                                 
                    t = data[j + n + 2];
                    data[j + n + 2] = data[k + n + 2];
                    data[k + n + 2] = t;

                    t = data[j + n + 3];
                    data[j + n + 3] = data[k + n + 3];
                    data[k + n + 3] = t;
                }
                // Knuth R3: advance k                                                                           
                k += 4;
                if (k >= n)
                    break;
                // Knuth R4: advance j                                                                           
                var h = top;
                while (j >= h)
                {
                    j -= h;
                    h /= 2;
                }
                j += h;
            } // bit reverse loop                                                                                
        }

        static void Reverse(double[] RlDat, double[] ImDat, int n)
        {
            // bit reverse the indices. This is exercise 5 in section                                            
            // 7.2.1.1 of Knuth's TAOCP the idea is a binary counter                                             
            // in k and one with bits reversed in j                                                              
            int j = 0, k = 0; // Knuth R1: initialize                                                            
            var top = n / 2;  // this is Knuth's 2^(n-1)                                                         
            while (true)
            {
                int j_l = j / 2, k_l = k / 2;

                // Knuth R2: swap - swap j+1 and k+2^(n-1), 2 entries each                                       
                //var t = data[j + 2];
                var t = RlDat[j_l + 1];
                RlDat[j_l + 1] = RlDat[k_l + top];
                RlDat[k_l + top] = t;


                t = ImDat[j_l + 1];
                ImDat[j_l + 1] = ImDat[k_l + top];
                ImDat[k_l + top] = t;

                if (j > k)
                { // swap two more                                     

                    // j and k                                                                                   
                    t = RlDat[j_l];
                    RlDat[j_l] = RlDat[k_l];
                    RlDat[k_l] = t;

                    t = ImDat[j_l];
                    ImDat[j_l] = ImDat[k_l];
                    ImDat[k_l] = t;

                    // j + top + 1 and k+top + 1                                                                 
                    t = RlDat[j_l + top + 1];
                    RlDat[j_l + top + 1] = RlDat[k_l + top + 1];
                    RlDat[k_l + top + 1] = t;

                    t = ImDat[j_l + top + 1];
                    ImDat[j_l + top + 1] = ImDat[k_l + top + 1];
                    ImDat[k_l + top + 1] = t;
                }
                // Knuth R3: advance k                                                                           
                k += 4;
                if (k >= n)
                    break;
                // Knuth R4: advance j                                                                           
                var h = top;
                while (j >= h)
                {
                    j -= h;
                    h /= 2;
                }
                j += h;
            } // bit reverse loop                                                                                
        }

        /// <summary>                                                                                            
        /// Pre-computed sine/cosine tables for speed                                                            
        /// </summary>                                                                                           
        double[] cosTable;
        double[] sinTable;

        #endregion

    }
}
/*  (C) Copyright Michael Godfrey 2017
 * 
 *      This file is free software, redistributed as part of Lightsail.DSP
 *
 *      See /COPYING for more details
 *  
 *  Chris Lomont is a handsome devil. 
 */
