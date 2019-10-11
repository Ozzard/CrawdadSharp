using System;
using System.Diagnostics;

namespace CrawdadSharp
{
    class ChromSmoother
    {
        /// a threshold for suppressing short spike signals
        public float baseline;

        ///if the signal is below 'minimum_baseline' for minimum_spike_len units, then reject the peak
        public int spike_len;

        public ChromSmoother()
            : this(0)
        {
        }

        public ChromSmoother(int weight_size)
        {
            baseline = 0.0f;
            spike_len = 1;
            fft = false;
            set_weight_size(weight_size);
        }

        public float[] Weights
        {
            get
            {
                return weights;
            }
            protected set
            {
                Debug.Assert(value.Length == weights.Length);
                Array.Copy(value, weights, value.Length);
            }
        }

        public void invert_weights()
        {
            for (int i = 0; i < weights.Length; i++)
                weights[i] *= -1;
        }

        protected float[] weights;
        protected bool fft;
        protected int half_window_size;
        protected int window_size;

        protected void set_weight_size(int weight_size)
        {
            weights = new float[weight_size];
            half_window_size = weight_size / 2;
            window_size = weight_size;
        }

        private void smooth_vect_discrete(float[] raw_vec, float[] out_vec)
        {

            Debug.Assert(raw_vec.Length == out_vec.Length);
            if (raw_vec.Length <= half_window_size * 2 + 1)
            {
                Array.Copy(raw_vec, out_vec, raw_vec.Length);
            }
            else
            {
                for (int i = 0; i < half_window_size; i++)
                    out_vec[i] = raw_vec[i];
                for (int i = raw_vec.Length - half_window_size; i < (uint)raw_vec.Length; i++)
                    out_vec[i] = raw_vec[i];

                /* we assume the weights are normalized */
                for (int i = half_window_size; i < raw_vec.Length - half_window_size; i++)
                {
                    float t = 0.0f;
                    for (int offset = 0; offset < window_size; offset++)
                    {
                        int raw_idx = i - half_window_size + offset;
                        t += raw_vec[raw_idx] * weights[offset];
                    }
                    out_vec[i] = t;
                }
            }
        }

        ///Filters chromatogram by removing small spikes, where the signal is <= baseline over at least spike_len consecutive units.
        ///Other approaches might be : if in a window, at least M/N points are at 'baseline', then set all N points to 'baseline'.
        ///A bandpass filter may be a better approach. The motivation behind this function is that the matched filtration with
        ///gaussian or savitzky-golay may cause smoothing artifacts.

        private void smooth_vect_fft(float[] raw_vec, float[] out_vec)
        {
            throw new Exception("Forget about trying to call smooth_vect_fft if you don't have FFTPACK");
        }

        public void smooth_vect(float[] raw_vec, float[] out_vec)
        {
            if (fft)
                smooth_vect_fft(raw_vec, out_vec);
            else
                smooth_vect_discrete(raw_vec, out_vec);
        }

        public void smooth_vect(float[] to_smooth)
        {
            float[] tmp_vect = new float[to_smooth.Length];
            smooth_vect(to_smooth, tmp_vect);
            Array.Copy(tmp_vect, to_smooth, to_smooth.Length);
        }
    }
}
