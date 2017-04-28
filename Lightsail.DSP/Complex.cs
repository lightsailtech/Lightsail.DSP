// Revised: 4/25/2017

namespace Lightsail.DSP
{
    /// <summary>
    /// 
    /// </summary>
    public struct Complex
    {
        public Complex(double rl, double im)
        {
            RlPart = rl; ImPart = im;
        } // end Complex

        public double RlPart { get; set; }
        public double ImPart { get; set; }

        public double Mag { get { return System.Math.Sqrt(RlPart * RlPart + ImPart * ImPart); } }
        public double Phase { get { return System.Math.Atan(ImPart / RlPart); } }

    } // end struct Complex

}

/*  (C) Copyright Michael Godfrey 2017
 * 
 *      This file is free software, part of Lightsail.DSP
 *
 *      See /COPYING for more details
 *  
 *  Imagine your life but lived as another race. There are people too on the other side of the mountain.
 */
   