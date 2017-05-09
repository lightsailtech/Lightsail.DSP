/* (C) Copyright Michael Godfrey 2017
 * 
 * Revised: 4/25/2017
 */
using Lomont;
using System;

namespace Lightsail.DSP
{
    public interface ILightsailFFT
    {
               void FFT(Complex[ ] data, bool forward = true);

          void FFTTable(Complex[ ] data, bool forward = true);

             void FFT2D(Complex[,] data, bool forward = true);

        void FFT2DTable(Complex[,] data, bool forward = true);

              void FFT(double[ ] datRl, double[ ] datIm, bool forward = true);

         void FFTTable(double[ ] datRl, double[ ] datIm, bool forward = true);

            void FFT2D(double[,] datRl, double[,] datIm, bool forward = true);

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

        private LomontFFT fft;

        public LightsailFFT()
        {
            fft = new LomontFFT();
        } // end constructor

        public void FFT(Complex[] data, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFT(double[] datRl, double[] datIm, bool forward = true)
        {
            fft.NrmlFFT(datRl, datIm, forward);
        }

        public void FFTTable(Complex[] data, bool forward = true)
        {
            throw new NotImplementedException();
        }

        public void FFTTable(double[] datRl, double[] datIm, bool forward = true)
        {
            fft.NrmlTableFFT(datRl, datIm, forward);
        } // end FFTTable

        #region 2D FFT Methods

        public void FFT2D(Complex[,] data, bool forward = true)
        {
            int idx, jdx,
                nx = data.GetLength(0),
                ny = data.GetLength(1);

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
                    rl_x[idx] = data[idx, jdx].RlPart;
                    im_x[idx] = data[idx, jdx].ImPart;
                }

                FFT(rl_x, im_x, forward: forward);

                for (idx = 0; idx < nx; idx++)
                {
                    data[idx, jdx].RlPart = rl_x[idx];
                    data[idx, jdx].ImPart = im_x[idx];
                }
            }

            // Row Transform
            // ================================================================
            for (idx = 0; idx < nx; idx++)
            {
                for (jdx = 0; jdx < ny; jdx++)
                {
                    rl_y[jdx] = data[idx, jdx].RlPart;
                    im_y[jdx] = data[idx, jdx].ImPart;
                }

                FFT(rl_y, im_y, forward: forward);

                for (jdx = 0; jdx < ny; jdx++)
                {
                    data[idx, jdx].RlPart = rl_y[jdx];
                    data[idx, jdx].ImPart = im_y[jdx];
                }
            }
        } // end FFT2D

        public void FFT2D(double[,] datRl, double[,] datIm, bool forward = true)
        {
            int idx, jdx,
                nx = datRl.GetLength(0),
                ny = datRl.GetLength(1);

            if (nx == 0 || ny == 0)
                return;

            double[]
                rl_x = new double[nx],
                im_x = new double[nx],
                rl_y = new double[ny],
                im_y = new double[ny];

            if (datRl.GetLength(0) != datIm.GetLength(0)
             || datRl.GetLength(1) != datIm.GetLength(1))
                throw new ArgumentException("Real and Imaginary component"
                    + " arrays must be the same size");
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
                    rl_x[idx] = datRl[idx, jdx];
                    im_x[idx] = datIm[idx, jdx];
                }

                FFT(rl_x, im_x, forward: forward);

                for (idx = 0; idx < nx; idx++)
                {
                    datRl[idx, jdx] = rl_x[idx];
                    datIm[idx, jdx] = im_x[idx];
                }
            }

            // Row Transform
            // ================================================================
            for (idx = 0; idx < nx; idx++)
            {
                for (jdx = 0; jdx < ny; jdx++)
                {
                    rl_y[jdx] = datRl[idx, jdx];
                    im_y[jdx] = datIm[idx, jdx];
                }

                FFT(rl_y, im_y, forward: forward);

                for (jdx = 0; jdx < ny; jdx++)
                {
                    datRl[idx, jdx] = rl_y[jdx];
                    datIm[idx, jdx] = im_y[jdx];
                }
            }
        }

        public void FFT2DTable(Complex[,] data, bool forward = true)
        {
            int idx, jdx,
                nx = data.GetLength(0),
                ny = data.GetLength(1);

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
                    rl_x[idx] = data[idx, jdx].RlPart;
                    im_x[idx] = data[idx, jdx].ImPart;
                }

                FFTTable(rl_x, im_x, forward: forward);

                for (idx = 0; idx < nx; idx++)
                {
                    data[idx, jdx].RlPart = rl_x[idx];
                    data[idx, jdx].ImPart = im_x[idx];
                }
            }

            // Row Transform
            // ================================================================
            for (idx = 0; idx < nx; idx++)
            {
                for (jdx = 0; jdx < ny; jdx++)
                {
                    rl_y[jdx] = data[idx, jdx].RlPart;
                    im_y[jdx] = data[idx, jdx].ImPart;
                }

                FFTTable(rl_y, im_y, forward: forward);

                for (jdx = 0; jdx < ny; jdx++)
                {
                    data[idx, jdx].RlPart = rl_y[jdx];
                    data[idx, jdx].ImPart = im_y[jdx];
                }
            }
        }

        public void FFT2DTable(double[,] datRl, double[,] datIm, bool forward = true)
        {
            int idx, jdx,
                nx = datRl.GetLength(0),
                ny = datRl.GetLength(1);

            if (nx == 0 || ny == 0)
                return;

            double[]
                rl_x = new double[nx],
                im_x = new double[nx],
                rl_y = new double[ny],
                im_y = new double[ny];

            if (datRl.GetLength(0) != datIm.GetLength(0)
             || datRl.GetLength(1) != datIm.GetLength(1))
                throw new ArgumentException("Real and Imaginary component"
                    + " arrays must be the same size");
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
                    rl_x[idx] = datRl[idx, jdx];
                    im_x[idx] = datIm[idx, jdx];
                }

                FFTTable(rl_x, im_x, forward: forward);

                for (idx = 0; idx < nx; idx++)
                {
                    datRl[idx, jdx] = rl_x[idx];
                    datIm[idx, jdx] = im_x[idx];
                }
            }

            // Row Transform
            // ================================================================
            for (idx = 0; idx < nx; idx++)
            {
                for (jdx = 0; jdx < ny; jdx++)
                {
                    rl_y[jdx] = datRl[idx, jdx];
                    im_y[jdx] = datIm[idx, jdx];
                }

                FFTTable(rl_y, im_y, forward: forward);

                for (jdx = 0; jdx < ny; jdx++)
                {
                    datRl[idx, jdx] = rl_y[jdx];
                    datIm[idx, jdx] = im_y[jdx];
                }
            }
        }

        #endregion

    } // end class LightsailFFT
}

/*  (C) Copyright Michael Godfrey 2017
 * 
 *      This file is free software, part of Lightsail.DSP
 *
 *      See /COPYING for more details
 *  
 *  Peace be with you. Call your parents. You're not as bad as you think you are, when you think you are.
 */
