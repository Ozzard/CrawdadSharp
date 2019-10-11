namespace CrawdadSharp
{
#if false
    class UnevenChromSmoother
    {
        /// (RT in minutes) if the signal is below 'minimum_baseline' for minimum_spike_len units, then reject the peak
        public float minimum_spike_len;
        /// a threshold for suppressing short spike signals
        public float minimum_baseline;
        /// size of the window in retention time units
        public double rt_filter_half_size;
        public double rt_per_filter_pt;
        /// number of points over which to compute the gaussian weights
        public const int sample_half_width = 5000;
        public const int sample_window_width = sample_half_width * 2 + 1;
        //sample size has to be odd to includ
        public float[] weights;
        public bool fft { get; set; }

        public void init()
        {
            minimum_baseline = 0.0f;
            minimum_spike_len = 0.0f;
        }

        public UnevenChromSmoother(double rt_filter_len)
        {
            rt_filter_half_size = rt_filter_len / 2.0;
            rt_per_filter_pt = rt_filter_len / sample_window_width;
            weights = new float[sample_window_width];
            //re-use the related class
            ChromSmoother weight_factory = new ChromSmoother(sample_window_width);
            weights = weight_factory.Weights;
        }

        public void smooth_vect(float[] raw_vec, float[] RT, ref float[] out_vec);

        public void smooth_vect(float[] smooth_vect, float[] RT);

        void resize(int new_size)
        {
            weights.resize(new_size);
        }

        public void invert_weights()
        {
            for (int i = 0; i < weights.Length; i++)
                weights[i] *= -1;
        }
    }
#endif
}
