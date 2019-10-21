namespace CrawdadSharp
{
    public class SlimCrawPeak
    {
        public int peak_rt_idx, start_rt_idx, stop_rt_idx, max_rt_idx;
        public int mz_idx;
        public int len;
        public float fwhm;
        public bool fwhm_calculated_ok;
        public float fwfpct;
        public float bg_area;
        public float raw_area; // total area under the curve, including background
        public float peak_area;
        public float bgslope;
        ///cutoff level for extending past the peak

        ///maximum height, calculated above background
        public float peak_height;
        public float raw_height;


        public SlimCrawPeak()
        {
            init();
        }

        public int get_rel_peak_idx() => peak_rt_idx - start_rt_idx;

        /*
        - How much of a case do we need to avoid where the 2D derivative crosses excessively quickly? What
            are the consequences of finding peaks in this fashion?
        - What parameters go into peak finding with this technqiue? Distribution size, which affects smoothing?.


        - It should be a lot easier to find and define baseline using this technique. It should
          work perfectly fine to go and look at frequency distributions of baseline intensity

        */
        protected float get_peak_to_bg()
        {

            if (bg_area <= 0.0f)
                return 99999.0f;
            return peak_area / bg_area;
        }

        public SlimCrawPeak(int start_idx, int stop_idx, int peak_idx, float[] raw, float[] _, int mz_idx)
        {
            peak_rt_idx = peak_idx;
            start_rt_idx = start_idx;
            stop_rt_idx = stop_idx;
            this.mz_idx = mz_idx;
            init();

            // begin erynes BUGBUG DEBUG
            if (raw.Length > 0)
                raw_area = (float)crawutils.area_under_curve(raw, start_idx, stop_idx);
            // end erynes BUGBUG DEBUG
        }

        void init()
        {
            len = stop_rt_idx - start_rt_idx + 1;
            raw_area = bg_area = peak_area = bgslope = peak_height = raw_height = -1.0f;
            fwhm = -1;
            fwhm_calculated_ok = false;
        }

        public float height_norm_slope()
        {
            return (bgslope * len) / peak_height;
        }

        ///internal method for calculating signal and background areas
        void get_sig_bg_areas(float[] raw, float[] bg_scratch)
        {
            float bg_val = raw[start_rt_idx];
            /* TODO -- trapezoidal area calculation in place, rather than copied to a vector and passed to a function */
            float next_used_val, curr_used_val;
            curr_used_val = next_used_val = bg_val;

            //#define FAST_WAY
#if FAST_WAY
		float peak_bg_val;
        int rel_peak_idx = get_rel_peak_idx();
        peak_bg_val = bg_val + rel_peak_idx* bgslope;
		for (int i = 0; i< this->len - 1; i++) {
			total_bg += next_used_val;
			curr_used_val = next_used_val;
			bg_val += bgslope;
			if (raw[start_rt_idx + i + 1] < bg_val) {
				next_used_val = raw[start_rt_idx + i + 1];
			}
			else {
				next_used_val = bg_val;
			}
total_bg += (next_used_val - curr_used_val) / 2.0;
		}
		peak_height = raw[peak_rt_idx] - peak_bg_val;
#else
            bg_val = raw[start_rt_idx];
            for (int i = 0; i < len; i++)
            {
                if (raw[start_rt_idx + i] < bg_val)
                {
                    bg_scratch[i] = raw[start_rt_idx + i];
                }
                else
                {
                    bg_scratch[i] = bg_val;
                }
                bg_val += bgslope;

            }
            bg_area = (float)crawutils.area_under_curve(bg_scratch, 0, len - 1);
            raw_height = raw[peak_rt_idx];
            peak_height = raw[peak_rt_idx] - bg_scratch[get_rel_peak_idx()];
#endif

            peak_area = raw_area - bg_area;

        }

        protected virtual string as_string_header()
        {
            return "mz\tstart_idx\tpeak_idx\tstop_idx\tfwhm";
        }

        public virtual string as_string()
        {
            return $"{mz_idx}\t{start_rt_idx}\t{peak_rt_idx}\t{stop_rt_idx}\t{fwhm:f2}";
        }

        string as_string_long()
        {
            return (as_string() + internal_as_string_long());
        }

        string internal_as_string_long()
        {
            return $"\t{peak_height:3f}\t{peak_area:3f}\t{len}\t{get_peak_to_bg():3f}";
        }

        string as_string_long_header()
        {
            return as_string_header() + "\tpeak_height\tpeak_area\tlen\tpeak_bg_ratio";
        }
    }
}
