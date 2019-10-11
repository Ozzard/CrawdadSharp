namespace CrawdadSharp
{
    class CrawPeakLocated : CrawPeak
    {
        float mz, rt_start, rt_stop, rt_peak;
        CrawPeakLocated(int lh_valley, int rh_valley, int peak_loc, float[] chrom, float[] scratch_chrom, int mz_idx)
             : base(lh_valley, rh_valley, peak_loc, chrom, scratch_chrom, mz_idx)
        {
        }

        protected override void init()
        {
            mz = rt_start = rt_stop = rt_peak = -1.0f;
            base.init();
        }

        void set_rt_mz(float mz, float rt_start, float rt_peak, float rt_stop)
        {
            this.mz = mz;
            this.rt_peak = rt_peak;
            this.rt_start = rt_start;
            this.rt_stop = rt_stop;
        }

        public override string as_string()
        {
            return $"{mz:f4.4}\t{rt_start:f3.3}\t{rt_peak:f3.3}\t{rt_stop:f3.3}\t{start_rt_idx}\t{peak_rt_idx}\t{stop_rt_idx}";
        }

        protected override string as_string_header()
        {
            return "mz\tstart_rt\tpeak_rt\tstop_rt\tstart_idx\tpeak_idx\tstop_idx";
        }
    }
}