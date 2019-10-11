using System;
using System.Linq;

namespace CrawdadSharp
{
    ///A simple structure for peaks
    class CrawPeak
        : SlimCrawPeak
    {
        ///index into the list of mzs to which this peak belongs


        ///indices to peak_rt, start_rt, stop_rt, rt of max intensity


        ///observed intensity values
        public float[] intensities;
        ///background intensity values 
        float[] background_vals;
        ///area below bg area

        ///slope between peak ends

        ///difference between signal and the median above baseline
        float mean_above_baseline;
        float max_above_baseline;
        ///standard dev. of the difference between 
        float stddev_mean_above_baseline;

        //float sharpness;
        //int mean_crossing_cnt;
        int baseline_p_mean_crossing_cnt;
        ///internal method to calculate stats on 

        ///braindead constructor for now
        public CrawPeak()
        {
            start_rt_idx = stop_rt_idx = peak_rt_idx = mz_idx = -1;
            init();
        }


        ///sharpness == area / len
        float get_area_sharpness() => peak_area / len;
        ///height / len
        float get_height_sharpness() => peak_height / len;

        int get_baseline_p_mean_crossing() => baseline_p_mean_crossing_cnt;






        ///calculates slope of the background level as estimated from peak boundaries


        ///calculates nearness of peak location to peak edges


        ///returns peak index as a measure of scans from the leftmost boundary

        ///internal method for calclating peak height

        string as_string_long_header()
        {
            return as_string_header() + "\tpeak_height\tpeak_area\tlen\tpeak_bg_ratio\tmean_above_baseline\tstddev_mean_above_baseline\tmean_crossing\tasymmetry";
        }

        string as_string_long()
        {
            return as_string() + internal_as_string_long();
        }
        
        string internal_as_string_long()
        {
            return $"\t{peak_height:f3.3}\t{peak_area:f3.3}\t{len}\t{get_peak_to_bg():f3.3}\t{mean_above_baseline:f3.3}\t{stddev_mean_above_baseline:f3.3}\t{get_baseline_p_mean_crossing()}\t{assymmetry_stab():f3.3}";
        }

        float assymmetry_stab()
        {
            int peak_len = len;
            int min_distance = Math.Min(stop_rt_idx - peak_rt_idx, peak_rt_idx - start_rt_idx);
            return (min_distance * 1.0f / peak_len) * 2.0f;
        }

        ///extracts data corresponding to the peak's co-ordinates from a float vector to a target vector
        void extract_chrom_regions(float[] chrom, float[] target)
        {
            for (int i = 0; i < len; i++)
                target[i] = chrom[start_rt_idx + i];
        }

        void calc_CV()
        {
        }

        void calc_baseline_stats()
        {
            /* take all intensity values, take median or mean value
               - then use this as a springboard to calculate either CV, or number of crossing times across this point

               1. take baseline , calculate for each point the height above baseline (0 minimum) , then based upon this
                  extend the baseline up by such amount. This can then be used as baseline+(mean_above_baseline),
                  then maybe take number of points above that line


            */
            //TODO these intermediate vectors are probably unnecessary -- let's see what removing them would do..?


            ///difference between intensities and background values
            float[] height_delta = new float[intensities.Length];

            for (int i = 0; i < height_delta.Length; i++)
                height_delta[i] = intensities[i] - background_vals[i];

            float height_median = (float)crawstats.median(height_delta);

            ///background values plus median of height delta
            float[] baseline_plus_median = new float[intensities.Length];
            for (int i = 0; i < height_delta.Length; i++)
                baseline_plus_median[i] = background_vals[i] + height_median;

            //TODO find an efficient means for detecting sign crossing
            for (int i = 0; i < intensities.Length - 1; i++)
            {
                //calculate the number of times signal crosses, or sits on.. the height_delta vector shown before
                if (intensities[i] <= baseline_plus_median[i] && intensities[i + 1] >= baseline_plus_median[i + 1])
                    baseline_p_mean_crossing_cnt += 1;
                else if (intensities[i] >= baseline_plus_median[i] && intensities[i + 1] <= baseline_plus_median[i + 1])
                    baseline_p_mean_crossing_cnt += 1;
            }

            float[] baseline_median_vs_signal = new float[intensities.Length];
            for (int i = 0; i < intensities.Length; i++)
                baseline_median_vs_signal[i] = intensities[i] - baseline_plus_median[i];
            mean_above_baseline = (float)crawstats.mean(baseline_median_vs_signal);
            max_above_baseline = baseline_median_vs_signal.Max();
            stddev_mean_above_baseline = (float)Math.Sqrt(crawstats.var_w_mean(baseline_median_vs_signal, mean_above_baseline));
        }

        ///Constructor taking start,stop,peak,mz indices, and a vector of intensities
        public CrawPeak(int start_idx, int stop_idx, int peak_idx, float[] raw, float[] scratch, int mz_idx)
            : base(start_idx, stop_idx, peak_idx, raw, scratch, mz_idx)
        {
            init();
            intensities = new float[len];
            background_vals = new float[len];
            extract_chrom_regions(raw, intensities);
            calc_baseline_stats();
        }

        protected virtual void init()
        {
            stddev_mean_above_baseline = -1.0f;
            baseline_p_mean_crossing_cnt = 0;
        }
    }
}
