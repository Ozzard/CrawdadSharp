using System;

namespace CrawdadSharp
{
    internal class crawstats
    {
        public static double median(float[] v)
        {
            if (v.Length == 0)
                throw new Exception("madonna mia! cannot take median of an empty vector");
            if (v.Length == 1)
                return v[0];
            //sort to derive median value
            float[] v_c = new float[v.Length];
            Array.Copy(v, v_c, v.Length);

            Array.Sort(v_c);

            if (v_c.Length % 2 == 1)
                return v_c[(v_c.Length - 1) / 2];
            //even
            int rh_mid = v_c.Length / 2;
            return (v_c[rh_mid] + v_c[rh_mid - 1]) / 2.0;
        }

        ///sum of squares
        public static double ss(float[] v)
        {
            double t = 0.0;
            for (int i = 0; i < v.Length; i++)
                t += v[i] * v[i];
            return t;
        }

        public static double mean(float[] v)
        {
            double t = 0.0;
            for (int i = 0; i < v.Length; i++)
                t += v[i];
            return t / v.Length;
        }

        ///given a vector and its mean , returns the variance
        public static double var_w_mean(float[] v, double m)
        {
            double rt = 0.0;
            for (int i = 0; i < v.Length; i++)
            {
                double dev = (v[i] - m);
                rt += dev * dev;
            }
            return rt / (v.Length - 1);
        }
    }
}