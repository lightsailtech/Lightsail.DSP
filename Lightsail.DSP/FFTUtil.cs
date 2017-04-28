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
