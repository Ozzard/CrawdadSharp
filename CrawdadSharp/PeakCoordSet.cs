using System.Collections.Generic;

namespace CrawdadSharp
{
    /*! -- given a chromatogram, finds a series of 'CrawPeaks'
    */

    struct PeakCoordSet
    {
        int[] peak_start_idxs;
        int[] peak_peak_idxs;
        int[] peak_stop_idxs;

        public PeakCoordSet(IList<SlimCrawPeak> peaks)
        {
            peak_start_idxs = new int[peaks.Count];
            peak_peak_idxs = new int[peaks.Count];
            peak_stop_idxs = new int[peaks.Count];
            for (int i = 0; i < peaks.Count; i++)
            {
                peak_start_idxs[i] = peaks[i].start_rt_idx;
                peak_stop_idxs[i] = peaks[i].stop_rt_idx;
                peak_peak_idxs[i] = peaks[i].peak_rt_idx;
            }
        }

        ///returns the leftmost, and rightmost peak index which overlap with this peak
        public (int leftmost, int rightmost) find_overlap_bounds_by_peak_rt_idx(SlimCrawPeak p)
        {
            int lh_peak_idx = crawutils.get_lh_idx(peak_peak_idxs, p.start_rt_idx) + 1;
            int rh_peak_idx;
            for (rh_peak_idx = lh_peak_idx; rh_peak_idx < peak_peak_idxs.Length; rh_peak_idx++)
                if (p.stop_rt_idx < peak_peak_idxs[rh_peak_idx])
                    break;

            // TODO comparison between signed and unsigned integer expressions
            if (rh_peak_idx == peak_peak_idxs.Length)
                rh_peak_idx--;
            return (lh_peak_idx, rh_peak_idx);
        }
    };
}