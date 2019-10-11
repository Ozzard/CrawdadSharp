using System;
using System.Collections.Generic;

namespace CrawdadSharp
{
    class peak_lists
    {
        public readonly List<SlimCrawPeak> peaks_by_height;
        public readonly List<SlimCrawPeak> peaks_by_rt;

        public peak_lists(List<SlimCrawPeak> h_peaks, List<SlimCrawPeak> r_peaks)
        {
            peaks_by_height = h_peaks;
            peaks_by_rt = r_peaks;
        }

        public void resort_peaks()
        {
            peaks_by_height.Sort((lh, rh) => lh.peak_height.CompareTo(rh.peak_height));
            peaks_by_rt.Sort((lh, rh) => lh.start_rt_idx.CompareTo(rh.start_rt_idx));
        }

        public void clear_peak(SlimCrawPeak p)
        {
            ClearPeaksByHeight(p);
            ClearPeaksByStartRtIdx(p);
        }

        private void ClearPeaksByHeight(SlimCrawPeak p)
        {
            bool found = false;
            for (int i = 0; i < peaks_by_height.Count; i++)
            {
                if (peaks_by_height[i].peak_height > p.peak_height)
                    break;
                if (peaks_by_height[i].peak_height == p.peak_height)
                {
                    peaks_by_height.RemoveAt(i);
                    found = true;
                    break;
                }
            }
            if (!found)
                throw new Exception("peak height must match");
        }

        private void ClearPeaksByStartRtIdx(SlimCrawPeak p)
        {
            bool found = false;
            for (int i = 0; i < peaks_by_rt.Count; i++)
            {
                if (peaks_by_rt[i].start_rt_idx > p.start_rt_idx)
                    break;
                if (peaks_by_rt[i].start_rt_idx == p.start_rt_idx)
                {
                    peaks_by_rt.RemoveAt(i);
                    found = true;
                    break;
                }
            }
            if (!found)
                throw new Exception("peak start_rt_idx must match");
        }
    }
}