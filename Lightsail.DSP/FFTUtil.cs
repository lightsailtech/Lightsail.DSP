using System;
using System.Collections.Generic;
using System.Text;

namespace Lightsail.DSP
{
    public static class FFTUtil
    {
        // todo: find a way to do this in place
        public static double[] FFTShift(double[] vals)
        {
            var shiftDist = vals.Length / 2;
            double[] ret = new double[vals.Length];
            for (int idx = 0; idx < vals.Length; ++idx)
            {
                ret[shiftDist + idx >= vals.Length ? (shiftDist + idx) - vals.Length : shiftDist + idx] = vals[idx];
            }
            return ret;
        } // end void FFTShift

        // todo: find a way to do this in place
        public static double[] FFTShiftC(double[] vals)
        {
            var shiftDist = vals.Length / 4;
            double[] ret = new double[vals.Length];
            for (int idx = 0; idx < vals.Length / 2; ++idx)
            {
                var newAddr = shiftDist + idx >= vals.Length / 2 ? (shiftDist + idx) - vals.Length / 2 : shiftDist + idx;

                ret[newAddr * 2] = vals[idx * 2];
                ret[newAddr * 2 + 1] = vals[idx * 2 + 1];
            }
            return ret;
        } // end void FFTShift

        // todo: find a way to do this in place
        public static double[] IFFTShiftC(double[] vals)
        {
            var shiftDist = vals.Length / 4;
            double[] ret = new double[vals.Length];
            for (int idx = 0; idx < vals.Length / 2; ++idx)
            {
                var newAddr = shiftDist + idx >= vals.Length / 2 ? (shiftDist + idx) - vals.Length / 2 : shiftDist + idx;

                ret[idx * 2] = vals[newAddr * 2];
                ret[idx * 2 + 1] = vals[newAddr * 2 + 1];
            }
            return ret;
        } // end void IFFTShift

        // todo: find a way to do this in place
        public static double[] IFFTShift(double[] vals)
        {
            var shiftDist = vals.Length / 2;
            double[] ret = new double[vals.Length];
            for (int idx = 0; idx < vals.Length; ++idx)
            {
                ret[idx] = vals[shiftDist + idx >= vals.Length ? (shiftDist + idx) - vals.Length : shiftDist + idx];
            }
            return ret;
        } // end void FFTShift

        /// <summary>
        /// TODO - Needs unit test
        /// </summary>
        /// <param name="dataIn"></param>
        /// <param name="dataOut"></param>
        public static void LomontToComplex(double[] dataIn, out Complex[] dataOut)
        {
            if (dataIn.Length % 2 != 0) //TODO - parameterize argument lengths in exception msg
                throw new ArgumentException("The double[] argument, 'dataIn' must have a length divisible by 2");

            dataOut = new Complex[dataIn.Length / 2];
            for (int idx = 0; idx < dataIn.Length; idx += 2)
                dataOut[idx / 2] = new Complex(dataIn[idx], dataIn[idx + 1]);
        } // end LomontToComplex

        /// <summary>
        /// TODO - Needs unit test
        /// </summary>
        /// <param name="dataIn"></param>
        /// <param name="rlOut"></param>
        /// <param name="imOut"></param>
        public static void LomontToRlIm(double[] dataIn, out double[] rlOut, out double[] imOut)
        {
            if (dataIn.Length % 2 != 0) //TODO - parameterize argument lengths in exception msg
                throw new ArgumentException("The double[] argument, 'dataIn' must have a length divisible by 2.");

            rlOut = new double[dataIn.Length / 2];
            imOut = new double[dataIn.Length / 2];
            for (int idx = 0; idx < dataIn.Length; idx += 2)
            {
                rlOut[idx / 2] = dataIn[idx];
                imOut[idx / 2] = dataIn[idx + 1];
            }
        } // end LomontToRlIm

        /// <summary>
        /// TODO - Needs unit test
        /// </summary>
        /// <param name="dataIn"></param>
        /// <param name="dataOut"></param>
        public static void ComplexToLomont(Complex[] dataIn, out double[] dataOut)
        {
            dataOut = new double[dataIn.Length * 2];
            for(int idx = 0; idx < dataIn.Length; ++ idx)
            {
                dataOut[(idx * 2)    ] = dataIn[idx].RlPart;
                dataOut[(idx * 2) + 1] = dataIn[idx].ImPart;
            }
        } // end ComplexToLomont

        /// <summary>
        /// TODO - Needs unit test
        /// </summary>
        /// <param name="rlIn"></param>
        /// <param name="imIn"></param>
        /// <param name="dataOut"></param>
        public static void RlImToLomont(double[] rlIn, double[] imIn, out double[] dataOut)
        {
            if (rlIn.Length != imIn.Length) //TODO - parameterize argument lengths in exception msg
                throw new ArgumentException("The double[] arguments, 'rlIn' and 'imIn' must have the same length.");

            dataOut = new double[rlIn.Length * 2];
            for(int idx = 0; idx < rlIn.Length; ++ idx)
            {
                dataOut[(idx * 2)    ] = rlIn[idx];
                dataOut[(idx * 2) + 1] = imIn[idx];
            }
        } // end RlImToLomont

        /// <summary>
        /// TODO - Needs unit test
        /// </summary>
        /// <param name="rlIn"></param>
        /// <param name="imIn"></param>
        /// <param name="dataOut"></param>
        public static double[] RlImToLomont(double[] rlIn, double[] imIn)
        {
            double[] dataOut;
            RlImToLomont(rlIn, imIn, out dataOut);
            return dataOut;
        } // end RlImToLomont
    } // end class FFTUtil
}
/*  (C) Copyright Michael Godfrey 2017
 * 
 *      This file is free software, part of Lightsail.DSP
 *
 *      See /COPYING for more details
 *  
 *  Have fun, don't do evil. Try a little harder next time. Tchau, and I love you.
 */
