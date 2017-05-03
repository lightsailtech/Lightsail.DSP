/* (C) Copyright Michael Godfrey 2017
 * 
 * Revised: 4/25/2017
 */
using System;

namespace Lightsail.DSP
{
    public interface ILightsailFFT
    {
        void FFT(        Complex[] data, bool forward = true);
        void FFTReal(    Complex[] data, bool forward = true);
        void FFTTable(   Complex[] data, bool forward = true);
        void FFT2D(     Complex[,] data, bool forward = true);
        void FFT2DReal( Complex[,] data, bool forward = true);
        void FFT2DTable(Complex[,] data, bool forward = true);

        void FFT(       double[]  datRl, double[]  datIm, bool forward = true);
        void FFTReal(   double[]  datRl, double[]  datIm, bool forward = true);
        void FFTTable(  double[]  datRl, double[]  datIm, bool forward = true);
        void FFT2D(     double[,] datRl, double[,] datIm, bool forward = true);
        void FFT2DReal( double[,] datRl, double[,] datIm, bool forward = true);
        void FFT2DTable(double[,] datRl, double[,] datIm, bool forward = true);

        /// <summary>                                                                                            
        ///     Common (A,B) values are: 
        ///        [ 0, 1] - default                                                                               
        ///        [ 1,-1] - signal processing                                                                     
        ///        [-1, 1] - data processing                                                                       
        ///     Common 'A' values are: 1, 0, or -1                                                                   
        /// </summary>                                                                                           
        int A { get; set; }

        /// <summary>                                                                                            
        ///     Common (A,B) values are: 
        ///        [ 0, 1] - default                                                                               
        ///        [ 1,-1] - signal processing                                                                     
        ///        [-1, 1] - data processing                                                                                             
        ///     Math.Abs(B) should be roughly prime to N.                                                              
        ///     B=-1 corresponds to conjugating both input and output data.                                                                                         
        ///     Common 'B' values are: 1 or -1.                                                                      
        /// </summary>                            
        int B { get; set; }
    } // end interface ILightsailFFT

    public class LightsailFFT : ILightsailFFT
    {
        public int A { get { return fft.A; } set { fft.A = value; } }
        public int B { get { return fft.B; } set { fft.B = value; } }

        Lomont.LomontFFT fft;

        public LightsailFFT()
        {
            fft = new Lomont.LomontFFT();
        } // end constructor


        public void FFT(double[] data, bool forward)
        {
            fft.FFT(data, forward);
        } // end FFT

        public void FFTReal(double[] data, bool forward)
        {
            fft.RealFFT(data, forward);
        } // end FFTReal

        public void FFTTable(double[] data, bool forward)
        {
            fft.TableFFT(data, forward);
        } // end FFTTable

        public void FFT2(Complex[,] c, bool forward)
        {

            int idx, jdx;
            int m;

            int nx = c.GetLength(0), 
                ny = c.GetLength(1);

            if (nx == 0 || ny == 0)
                return;

            double[]
                rl_x = new double[nx],
                im_x = new double[nx],
                rl_y = new double[ny],
                im_y = new double[ny];

            if ((nx & (nx - 1)) != 0)
                throw new ArgumentException("x-data length " + nx 
                    + " in FFT2 is not a power of 2");

            if ((ny & (ny - 1)) != 0)
                throw new ArgumentException("y-data length " + ny 
                    + " in FFT2 is not a power of 2");

            // Row Transform
            // ================================================================
            for (jdx = 0; jdx < ny; jdx++)
            {
                for (idx = 0; idx < nx; idx++)
                {
                    rl_x[idx] = c[idx,jdx].RlPart;
                    im_x[idx] = c[idx,jdx].ImPart;
                }

                fft.FFT(FFTUtil.RlImToLomont(rl_x, im_x), forward: forward);

                for (idx = 0; idx < nx; idx++)
                {
                    c[idx,jdx].RlPart = rl_x[idx];
                    c[idx,jdx].ImPart = im_x[idx];
                }
            }

            // Row Transform
            // ================================================================
            for (idx = 0; idx < nx; idx++)
            {
                for (jdx = 0; jdx < ny; jdx++)
                {
                    rl_y[jdx] = c[idx,jdx].RlPart;
                    im_y[jdx] = c[idx,jdx].ImPart;
                }

                //FFT(dir, m, rl_y, im_y);

                for (jdx = 0; jdx < ny; jdx++)
                {
                    c[idx,jdx].RlPart = rl_y[jdx];
                    c[idx,jdx].ImPart = im_y[jdx];
                }
            }
        } // end FFT2

        /**
         * 
         */
        public void IFFT2()
        {

        } // end FFT2

        public void FFT(Complex[] data, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFTReal(Complex[] data, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFTTable(Complex[] data, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFT2D(Complex[,] data, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFT2DReal(Complex[,] data, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFT2DTable(Complex[,] data, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFT(double[] datRl, double[] datIm, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFTReal(double[] datRl, double[] datIm, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFTTable(double[] datRl, double[] datIm, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFT2D(double[,] datRl, double[,] datIm, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFT2DReal(double[,] datRl, double[,] datIm, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFT2DTable(double[,] datRl, double[,] datIm, bool forward = true)
        {
            throw new NotImplementedException();
        }
    } // end class FFT
}

/*  (C) Copyright Michael Godfrey 2017
 * 
 *      This file is free software, part of Lightsail.DSP
 *
 *      See /COPYING for more details
 *  
 *  Peace be with you. Call your parents. You're not as bad as you think you are, when you think you are.
 */
