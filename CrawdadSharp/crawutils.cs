using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace CrawdadSharp
{
    internal class crawutils
    {
        public static int CompByStartRTIdx(SlimCrawPeak x, SlimCrawPeak y) => x.start_rt_idx.CompareTo(y.start_rt_idx);
        public static int CompByPeakRTIdx(SlimCrawPeak x, SlimCrawPeak y) => x.peak_rt_idx.CompareTo(y.peak_rt_idx);
        public static int CompByPeakRawHeight(PeaksUsedOrNot x, PeaksUsedOrNot y) => x.peak_raw_height.CompareTo(y.peak_raw_height);

        public static double area_under_curve(float[] vals, int start_idx, int stop_idx)
        {

            double sum = 0.0;
            for (int i = start_idx; i <= stop_idx - 1; i++)
            {
                sum += vals[i]; //rectangle to right of current value
                sum += (vals[i + 1] - vals[i]) / 2.0;
            }
            return sum;
        }

        public static double area_under_curve(float[] vals)
        {
            return area_under_curve(vals, 0, vals.Length - 1);
        }

        ///returns maximum index within bounds [ start, stop ) (not inclusive of stop)
        public static int max_idx_vect_bound(float[] v, int start, int stop)
        {
            int max_idx = start;
            float max = 0;
            for (int i = start; i < stop; i++)
            {
                if (max < v[i])
                {
                    max = v[i];
                    max_idx = i;
                }
            }
            return max_idx;
        }

        ///returns index of minimum values within bounds [ start, stop ) (not inclusive of stop)
        public static int min_idx_vect_bound(float[] v, int start, int stop)
        {
            int min_idx = start;
            float min = 0;
            for (int i = start; i < stop; i++)
            {
                if (min > v[i])
                {
                    min = v[i];
                    min_idx = i;
                }
            }
            return min_idx;
        }

        ///is the vector monotonically increasing?				
        public static bool monotonic_vect(IList<float> v)
        {
            for (int i = 0; i < v.Count - 1; i++)
                if (v[i + 1] < v[i])
                    return false;
            return true;
        }
        public static bool monotonic_vect(IList<int> v)
        {
            for (int i = 0; i < v.Count - 1; i++)
                if (v[i + 1] < v[i])
                    return false;
            return true;
        }

        /** get_lh_idx --
           finds the index(lower_bound) where an element would
           be inserted into a vector -- this should give the left-hand side
           of the lower bound for an element to be inserted into a list  -

           Note -- this interprets lower_bound correctly unlike some other items of the code...

        */
        public static int get_lh_idx(IList<int> v, int lookup)
        {
            Debug.Assert(monotonic_vect(v));
            int lh = lower_bound(v, lookup);
            if (lh == v.Count)
                return v.Count - 1;
            if (lookup == v[lh])
                return lh;
            if (lh == 0)
                return 0;
            return lh - 1;
        }

        public static int get_lh_idx(float[] v, float lookup)
        {
            Debug.Assert(monotonic_vect(v));
            int lh = lower_bound(v, lookup);
            if (lh == v.Length)
                return v.Length - 1;
            if (lookup == v[lh])
                return lh;
            if (lh == 0)
                return 0;
            return lh - 1;
        }

        public static int get_lh_idx(int[] v, int lookup)
        {
            Debug.Assert(monotonic_vect(v));
            int lh = lower_bound(v, lookup);
            if (lh == v.Length)
                return v.Length - 1;
            if (lookup == v[lh])
                return lh;
            if (lh == 0)
                return 0;
            return lh - 1;
        }

        public static int lower_bound<T>(IList<T> list, T value, Comparison<T> comparison)
        {
            // TODO: There could reasonably be a C# 8.0 version of this with Slice
            for (int i = 0; i < list.Count; i++)
                if (comparison(list[i], value) >= 0)
                    return i;
            return list.Count;
        }

        public static int lower_bound(IList<int> list, int value)
        {
            // TODO: There could reasonably be a C# 8.0 version of this with Slice
            for (int i = 0; i < list.Count; i++)
                if (list[i] >= value)
                    return i;
            return list.Count;
        }

        public static int lower_bound(int[] ary, int value)
        {
            // TODO: There could reasonably be a C# 8.0 version of this with Slice
            for (int i = 0; i < ary.Length; i++)
                if (ary[i] >= value)
                    return i;
            return ary.Length;
        }

        public static int lower_bound(float[] ary, float value)
        {
            // TODO: There could reasonably be a C# 8.0 version of this with Slice
            for (int i = 0; i < ary.Length; i++)
                if (ary[i] >= value)
                    return i;
            return ary.Length;
        }
    }
}