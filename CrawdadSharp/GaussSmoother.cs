using System;
using System.Diagnostics;
using System.Linq;

namespace CrawdadSharp
{
    class GaussSmoother : ChromSmoother
    {

        public GaussSmoother(int filt_size)
            : base(filt_size)
        {
        }

        public GaussSmoother()
            : base()
        {
        }

        ///defines the weights based upon the SD of a gaussian, and the derivative
        public void set_gauss_weights(float sd, int derivative)
        {
            //hw should be width at half-height
            int hw = (int)(4.0 * (sd + 0.5));
            int wlen = hw * 2 + 1;
            /* do the same thing as python where the size of the window is set to 99% of the total area? */
            set_weight_size(wlen);
            float[] weights = new float[wlen]; // Initialised to 0.0
            weights[hw] = 1.0f;
            float sum = weights[hw];
            for (int i = 1; i < hw + 1; i++)
            {
                float t = (float)Math.Exp(-0.5f * i * i / sd);
                weights[hw + i] = t;
                weights[hw - i] = t;
                sum += t * 2;
            }
            for (int i = 0; i < wlen; i++)
                weights[i] /= sum;
            if (derivative > 0)
            {
                if (derivative == 1)
                {
                    weights[hw] = 0.0f;
                    for (int i = 1; i < hw + 1; i++)
                    {
                        float tmp = (i * -1.0f / sd) * weights[hw + i];
                        weights[hw + i] = tmp * -1.0f;
                        weights[hw - i] = tmp;
                    }
                }
                else if (derivative == 2)
                {
                    weights[hw] *= -1.0f / sd;
                    for (int i = 1; i < hw + 1; i++)
                    {
                        float tmp = (i * i / sd - 1.0f) * weights[hw + i] / sd;
                        weights[hw + i] = tmp;
                        weights[hw - i] = tmp;
                    }
                }
                else if (derivative == 3)
                {
                    weights[hw] = 0.0f;
                    float sd2 = sd * sd;
                    for (int i = 1; i < hw + 1; i++)
                    {
                        /* TODO CHECK THIS FORMULA */
                        float tmp = (3.0f - i * i / sd) * i * weights[hw + i] / sd / sd;
                        weights[hw + i] = tmp * -1.0f;
                        weights[hw - i] = tmp;
                    }
                }
                else if (derivative > 3)
                    throw new Exception("gaussian derivative of greater than 3rd order not supported");
            }
            Weights = weights;
        }

        ///trim the weights at the endpoints where the values are [abs(v) <= frac * max(weights)]
        public void trim_weights_by_frac_max(float frac = 0.005f)
        {
            Debug.Assert(frac < 1.0f);
            int first_keep, last_keep;
            first_keep = last_keep = -1;
            float weights_max = weights.Max();
            float thresh = weights_max * frac;
            for (int i = 0; i < weights.Length; i++)
            {
                if (Math.Abs(weights[i]) >= thresh)
                {
                    first_keep = i;
                    break;
                }
            }
            for (int i = weights.Length - 1; i >= first_keep; i--)
            {
                if (Math.Abs(weights[i]) >= thresh)
                {
                    last_keep = i;
                    break;
                }
            }
            Debug.Assert(first_keep > -1 && last_keep >= first_keep);
            float[] new_weights = new float[last_keep - first_keep + 1];
            for (int i = 0; i < new_weights.Length; i++)
                new_weights[i] = weights[first_keep + i];
            set_weight_size(new_weights.Length);
            Weights = new_weights;
        }
    }
}
