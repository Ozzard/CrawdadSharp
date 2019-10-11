namespace CrawdadSharp
{
    public class CrawPeakMethod
    {
        const double FHWM_TO_SD_CONV = 2.35482;

        private float _sd;

        public float minimum_level;
        public bool saved_weights;
        public bool mean_cutoff;
        public int min_len;
        public int switch_len;

        public PeakBoundMethod bound_meth;
        public PeakPeakMethod peak_location_meth;
        public BackgroundEstimationMethod background_estimation_method;
        public ChromSmoothMethod chrom_smooth_method;

        public bool merge_peaks_list_based;
        public bool extend_peak_to_lower_bound;
        public bool extend_to_zero_crossing;

        public bool background_at_lower_boundary;
        public bool background_at_mean_boundary;
        public bool background_at_mid_boundary;

        public bool exclude_extension_overlaps_by_peakrt;
        public bool extend_from_peak_rt;

        public float extend_allowed_asymmetry;
        public float fraction_to_valley;

        public float one_peak_slope_merge_constraint;
        public float mean_slope_merge_constraint;
        ///decrease a peak limit back to fraction of the maximum value
        ///applies during peak extension
        public float ratchet_back_to_frac_maxval;

        public CrawPeakMethod()
        {
            init();
        }

        public float Sd
        {
            get { return _sd; }
            set
            {
                _sd = value;
                saved_weights = false;
            }
        }

        public float Fwhm
        {
            get { return _sd * (float)FHWM_TO_SD_CONV; }
            set
            {
                Sd = value / (float)FHWM_TO_SD_CONV;
            }
        }

        private void init()
        {
            Sd = 4.0f;
            saved_weights = false;
            mean_cutoff = false;
            minimum_level = 0.0f;
            min_len = 3;
            switch_len = 2;
            bound_meth = PeakBoundMethod.GAUSS_2D_BOUND;
            peak_location_meth = PeakPeakMethod.MAXIMUM_PEAK;
            chrom_smooth_method = ChromSmoothMethod.SAVITZKY_GOLAY_SMOOTHER;
            extend_peak_to_lower_bound = false;
            extend_to_zero_crossing = false;
            extend_allowed_asymmetry = 1.0f;
            fraction_to_valley = 0.0f;
            background_estimation_method = BackgroundEstimationMethod.PEAK_BOUNDARY_ESTIMATE;
            one_peak_slope_merge_constraint = 0.0f;
            mean_slope_merge_constraint = 0.0f;
            exclude_extension_overlaps_by_peakrt = false;
            ratchet_back_to_frac_maxval = -1.0f;
            extend_from_peak_rt = false;
            merge_peaks_list_based = false;
        }

        public void set_default_peak_opts()
        {
            init();
        }
    }
}
