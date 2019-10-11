using System;
using System.Collections.Generic;

namespace CrawdadSharp
{
    public abstract class BaseCrawPeakFinder
    {
        ///smoothing object for convoluting with gaussian 2nd derivative
        public bool slim;
        public CrawPeakMethod method;
        public CrawPeakAnnotator annotator;

        GaussSmoother gs_2d;
        GaussSmoother gs_1d;
        GaussSmoother gs_0d;

        ///chromatogram passed into object
        public float[] chrom;
        ///second derivative of the passed-in chromatogram
        public float[] chrom_2d;
        public float[] chrom_1d;
        float[] chrom_0d;
        ///scratch List for the background intensities
        protected float[] chrom_bg_scratch;

        ///peaks found on chromatogram
        protected int mz_idx;
        public List<int> minus_crosses;
        public List<int> plus_crosses;
        readonly List<int> peaks = new List<int>();
        readonly List<int> valleys = new List<int>();

        protected BaseCrawPeakFinder()
        {
            method = new CrawPeakMethod();
            mz_idx = -1;
            slim = false;
            annotator = new CrawPeakAnnotator(this);
        }

        protected BaseCrawPeakFinder(CrawPeakMethod m)
            : this()
        {
            method = m;
        }

        void set_weights()
        {
            set_weights_0d();
            set_weights_1d();
            set_weights_2d();
        }

        ///calls - peaks and edges
        void peaks_and_valleys(float[] c)
        {
            peaks_and_valleys_brendanx(c);
        }

        //stub for brendan's N-pass mean height determinator

        /* virtual methods */
        //making this a pure virtual class...
        protected abstract void peak_voodoo(int lh_valley, int rh_valley, int peak_loc);
        protected abstract int get_num_stored_peaks();
        protected abstract SlimCrawPeak get_peak_ptr(int idx);

        ///returns major, minor version numbers
        (int, int) Version => (0, 2);

        ///is the data >= a threshold over a certain length 
        /*!
        @param c - the vector to search over
        @param idx - start index to search
        @param threshold - threshold value to compare againts
        @param len - starting from idx, search len slots in c to check that
        values are gte than ("greater than or equal to") the threshold
        */
        bool consistent_gte(IList<float> c, int idx, float threshold = 0.0f, int len = 2)
        {
            //find the maximum point to search through the data
            int rh_bound = Math.Min(c.Count, idx + len);
            for (int i = idx; i < rh_bound; i++)
                if (c[i] < threshold)
                    return false;
            return true;
        }


        ///is the data lte a threshold over a certain length 
        /*!
        @param c - the vector to search over
        @param idx - start index to search
        @param threshold - threshold value to compare againts
        @param len - starting from idx, search len slots in c to check that
        values are lte than threshold
        */
        bool consistent_lte(IList<float> c, int idx, float threshold = 0.0f, int len = 2)
        {
            int rh_bound = Math.Min(c.Count, idx + len);
            for (int i = idx; i < rh_bound; i++)
                if (c[i] > threshold)
                    return false;
            return true;
        }

        ///finds points where a data series crosses a threshold value,
        ///populates plus_ and minus_ crosses
        /*!
        @param c - vector to search
        @param threshold - threshold value to search for crossing
        @post - plus_crosses is populated with indices where the second
        ///derivative crosses the threshold line going upwards. Minus crosses
        ///takes the opposite value
        */
        public (List<int> plus_crosses, List<int> minus_crosses) find_cross_points(float[] c, float threshold = 0.0f)
        {
            // Preallocate enough space to keep from doing slow capacity grows
            List<int> plus_crosses_local = new List<int>(c.Length);
            List<int> minus_crosses_local = new List<int>(c.Length);

            for (int i = 0; i < c.Length - (method.switch_len - 1); i++)
            {
                if (c[i] < threshold && consistent_gte(c, i + 1, threshold, method.switch_len))
                    plus_crosses_local.Add(i);
                else if (c[i] > threshold && consistent_lte(c, i + 1, threshold, method.switch_len))
                    minus_crosses_local.Add(i);
            }

            return (plus_crosses_local, minus_crosses_local);
        }

        ///finds points where the second derivative crosses zero - note that
        ///it must stay on the changed-to sign for at least two points
        void set_cross_points(float[] c, float threshold = 0.0f)
        {
            (plus_crosses, minus_crosses) = find_cross_points(c, threshold);
        }

        /*! Finds maxima and minima between crossing points of a signal as described in find_cross_points.
            sets 'peaks' and 'valleys' in the 2D signal in the peaks, valleys
            values

            TODO : this could be sped up
        */
        class PVIndexType
        {
            public PVIndexType()
            {
                index = -1;
                type = '\0';
            }

            public int index;
            public char type;
        }

        void peaks_and_valleys_brendanx(float[] c)
        {

            int iMinus = 0, lenMinus = minus_crosses.Count;
            int iPlus = 0, lenPlus = plus_crosses.Count;

            /* Merge sorted lists of plus and minus crosses */
            /* + and - always have to alternate, since we define in terms of crossing one line, although we have the rule that they can't pass
                 over one gap point apart  */
            List<PVIndexType> indexes = new List<PVIndexType>(lenMinus + lenPlus);
            while (iMinus < lenMinus || iPlus < lenPlus)
            {
                if (iMinus >= lenMinus ||
                    (iPlus < lenPlus && plus_crosses[iPlus] < minus_crosses[iMinus]))
                {
                    // Ignore consecutive plus crosses (fix so this can be an assert)
                    if (indexes.Count > 0 && indexes[indexes.Count - 1].type == 'p')
                    {
                        iPlus++;
                        continue;
                    }
                    indexes.Add(new PVIndexType { type = 'p', index = plus_crosses[iPlus++] });
                }
                else
                {
                    // Ignore consecutive minus crosses (fix so this can be an assert)
                    if (indexes.Count > 0 && indexes[indexes.Count - 1].type == 'm')
                    {
                        iMinus++;
                        continue;
                    }
                    indexes.Add(new PVIndexType { type = 'm', index = minus_crosses[iMinus++] });
                }
            }

            valleys.Clear();
            peaks.Clear();

            //valley at the beginning of the chrom
            int lenIndexes = indexes.Count;
            if (lenIndexes == 0 || indexes[0].type == 'p')
                valleys.Add(0);
            if (lenIndexes > 0)
            {
                for (int i = 0; i < lenIndexes - 1; i++)
                {
                    char trans_type = indexes[i].type;

                    if (trans_type == 'p')
                    {
                        int max_idx = crawutils.max_idx_vect_bound(c, indexes[i].index, indexes[i + 1].index);
                        peaks.Add(max_idx);
                        //find index of maximum value spanning from all_trans[idx] to all_trans[idx+1]
                    }
                    else if (trans_type == 'm')
                    {
                        //find index of minimum value spanning from all_trans[idx] to all_trans[idx+1]
                        int min_idx = crawutils.min_idx_vect_bound(c, indexes[i].index, indexes[i + 1].index);
                        valleys.Add(min_idx);
                    }
                }
                //adding a valley at the end of the chromatogram
                if (valleys.Count == peaks.Count)
                    valleys.Add(c.Length - 1);
            }
        }

        public void get_2d_chrom(float[] in_chrom, float[] out_chrom)
        {
            if (!method.saved_weights)
                set_weights();
            gs_2d.smooth_vect(in_chrom, out_chrom);
        }

        public void get_2d_chrom(float[] in_chrom)
        {
            if (!method.saved_weights)
                set_weights();
            gs_2d.smooth_vect(in_chrom);
        }

        public void get_1d_chrom(float[] in_chrom, float[] out_chrom)
        {
            if (!method.saved_weights)
                set_weights();
            gs_1d.smooth_vect(in_chrom, out_chrom);
        }

        void get_0d_chrom(float[] in_chrom)
        {
            if (!method.saved_weights)
                set_weights();
            gs_0d.smooth_vect(in_chrom);
        }

        void get_0d_chrom(float[] in_chrom, float[] out_chrom)
        {
            if (!method.saved_weights)
                set_weights();
            gs_0d.smooth_vect(in_chrom, out_chrom);
        }

        public void get_1d_chrom(float[] in_chrom)
        {
            if (!method.saved_weights)
                set_weights();
            gs_1d.smooth_vect(in_chrom);
        }

        void set_weights_2d()
        {
            if (null == gs_2d)
                init_smoothers();
            gs_2d.set_gauss_weights(method.Sd, 2);
            method.saved_weights = true;
            gs_2d.invert_weights();
            gs_2d.trim_weights_by_frac_max(0.005f);
        }

        void set_weights_1d()
        {
            if (null == gs_1d)
                init_smoothers();
            gs_1d.set_gauss_weights(method.Sd, 1);
            method.saved_weights = true;
            gs_1d.trim_weights_by_frac_max(0.005f);
        }

        void set_weights_0d()
        {
            if (null == gs_0d)
                init_smoothers();
            gs_0d.set_gauss_weights(method.Sd, 0);
            method.saved_weights = true;
            gs_0d.trim_weights_by_frac_max(0.005f);
        }

        public void call_peaks()
        {

            if (null == gs_2d)
                init_smoothers();

            if (chrom.Length == 0)
                throw new Exception("Come on! Set a chromatogram with set_chrom!");

            if (!method.saved_weights)
                set_weights();

            chrom_2d = (float[])chrom.Clone();
            chrom_1d = (float[])chrom.Clone();
            chrom_0d = (float[])chrom.Clone();

            gs_2d.smooth_vect(chrom_2d);
            gs_1d.smooth_vect(chrom_1d);
            gs_0d.smooth_vect(chrom_0d);

            //find 2nd derivative zero crossing points
            set_cross_points(chrom_2d);
            //annotated peaks and valleys in the 2nd derivative
            peaks_and_valleys(chrom_2d);

            float mean_height_cutoff = 0.0f;
            if (method.mean_cutoff)
                mean_height_cutoff = determine_mean_cutoff();

            for (int i = 0; i < peaks.Count; i++)
            {
                int lh_valley, rh_valley;
                int chrom2d_peak_loc = peaks[i];
                int lh_valley_idx = crawutils.get_lh_idx(valleys, chrom2d_peak_loc);
                if (valleys.Count == 1)
                {
                    if (valleys[0] > peaks[0])
                    {
                        lh_valley = 0;
                        rh_valley = valleys[0];
                    }
                    else
                    {
                        lh_valley = valleys[0];
                        rh_valley = chrom_2d.Length - 1;
                    }
                }

                else if (valleys[lh_valley_idx] > chrom2d_peak_loc)
                {
                    lh_valley = 0;
                    rh_valley = 1;
                }
                else if (lh_valley_idx >= valleys.Count - 1)
                {
                    lh_valley = valleys[lh_valley_idx - 1];
                    rh_valley = chrom_2d.Length - 1;
                }
                else
                {
                    lh_valley = valleys[lh_valley_idx];
                    rh_valley = valleys[lh_valley_idx + 1];
                }

                //now we expand from the peak locations to 

                //note that this alters lh_valley, rh_valley
                delimit_by_minimum_level(ref lh_valley, ref rh_valley, chrom2d_peak_loc);

                //now we have lh,peak,rh
                if (!(method.mean_cutoff) ||
                    (annotator.active_chrom[chrom2d_peak_loc] > mean_height_cutoff))
                {
                    if (passes_min_len(lh_valley, rh_valley))
                        peak_voodoo(lh_valley, rh_valley, chrom2d_peak_loc);
                }
            }
        }

        float determine_mean_cutoff()
        {
            double t = 0.0;
            for (int i = 0; i < peaks.Count; i++)
                t += annotator.active_chrom[peaks[i]];
            return (float)(t / peaks.Count);
        }

        CrawPeak construct_peak(int lh_valley, int rh_valley, int peak_loc)
        {
            return new CrawPeak(lh_valley, rh_valley, peak_loc, chrom, chrom_bg_scratch, mz_idx);
        }

        SlimCrawPeak construct_slim_peak(int lh_valley, int rh_valley, int peak_loc)
        {
            return new SlimCrawPeak(lh_valley, rh_valley, peak_loc, chrom, chrom_bg_scratch, mz_idx);
        }

        bool passes_min_len(int lh_valley, int rh_valley)
        {
            return (rh_valley - lh_valley + 1) >= method.min_len;
        }

        void delimit_by_minimum_level(ref int lh_valley_idx, ref int rh_valley_idx, int peak_loc)
        {
            for (int lh = peak_loc; lh > lh_valley_idx; lh--)
            {
                if (chrom[lh] <= method.minimum_level)
                {
                    lh_valley_idx = lh;
                    break;
                }
            }
            for (int rh = peak_loc; rh <= rh_valley_idx; rh++)
            {
                if (chrom[rh] <= method.minimum_level)
                {
                    rh_valley_idx = rh;
                    break;
                }
            }
        }

        void init_smoothers()
        {
            gs_0d = new GaussSmoother();
            gs_1d = new GaussSmoother();
            gs_2d = new GaussSmoother();
        }

        public void set_chrom(float[] f, int mz_idx)
        {
            this.mz_idx = mz_idx;
            chrom = (float[])f.Clone();
            chrom_bg_scratch = new float[f.Length];
            annotator.active_chrom = chrom;
        }

        public void clear()
        {
            chrom = null;
            chrom_2d = null;
            chrom_1d = null;
            mz_idx = -1;
            minus_crosses?.Clear();
            plus_crosses?.Clear();
            peaks.Clear();
            valleys.Clear();
            clear_sps();
        }

        protected abstract void clear_sps();
    }
}
