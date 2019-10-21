using System;
using System.Collections.Generic;
using System.Linq;

namespace CrawdadSharp
{
    public class CrawPeakAnnotator
    {
        const int BGSCRATCH_START_SIZE = 512;

        BaseCrawPeakFinder pf;
        public float[] active_chrom { get; set; }
        float[] bg_scratch;
        float? fixed_background;

        void init()
        {
            bg_scratch = new float[BGSCRATCH_START_SIZE];
            fixed_background = null;
        }

        /*! Calculation of peak, bg areas, bg_subtracted_peak_height for an arbitrary region using the suppeak algorithm as of 03/09 */

        //summary method for methods relating to quantitating peaks.
        void peak_annotate(SlimCrawPeak peak)
        {
            set_peak_slope(peak);
            //calls set_peak_areas onpeak member variables
            set_peak_bg_subtracted_area(peak);
            calc_fwhm(peak);
            calc_fwfpcnt(peak);
        }

        int calc_len(int start_rt_idx, int stop_rt_idx) => stop_rt_idx - start_rt_idx + 1;

        public CrawPeakAnnotator()
        {
            init();
        }

        public CrawPeakAnnotator(BaseCrawPeakFinder pf)
        {
            this.pf = pf;
            init();
        }

        int get_peakloc_in_range(int start_idx, int stop_idx)
        {
            if (pf.method.peak_location_meth == PeakPeakMethod.GAUSS_2D_PEAK)
                return _get_peakloc_by_2d(start_idx, stop_idx);
            if (pf.method.peak_location_meth == PeakPeakMethod.MAXIMUM_PEAK)
                return _get_peakloc_by_local_max(start_idx, stop_idx);
            if (pf.method.peak_location_meth == PeakPeakMethod.CENTROID_PEAK)
                return _get_peakloc_by_centroid(start_idx, stop_idx);
            throw new Exception("no legal peak location has been found");
        }

        int _get_peakloc_by_2d(int start_idx, int stop_idx)
        {
            return crawutils.max_idx_vect_bound(pf.chrom_2d, start_idx, stop_idx);
        }

        int _get_peakloc_by_centroid(int _1, int _2)
        {
            throw new NotImplementedException();
        }

        int _get_peakloc_by_local_max(int start_idx, int stop_idx)
        {
            return crawutils.max_idx_vect_bound(active_chrom, start_idx, stop_idx);
        }

        //ER - USE THIS TO JOIN PEAKS TOGETHER -
        ///adds the right-hand peak to the left-hand peak
        ///a version that returns a new peak would be useful too but that is low priority
        ///based on the determination to merge peaks, join to_glom to changed
        void glom_peak(SlimCrawPeak changed, SlimCrawPeak to_glom)
        {
            if (to_glom.start_rt_idx > changed.start_rt_idx)
                changed.stop_rt_idx = to_glom.stop_rt_idx;
            else
                changed.start_rt_idx = to_glom.start_rt_idx;
            changed.len = changed.stop_rt_idx - changed.start_rt_idx + 1;
            peak_tweak(changed, pf.method);
            set_peak_bg_subtracted_area(changed);
        }


        ///merge peaks based on adjacent peaks with similar slope patterns
        ///merges peaks by:
        ///start with the heighest peak -
        ///find nearest peak
        IList<SlimCrawPeak> merge_peaks_list_based(List<SlimCrawPeak> in_peaks)
        {
            if (!pf.method.merge_peaks_list_based)
                return new List<SlimCrawPeak>(in_peaks);

            List<SlimCrawPeak> peaks_by_height = new List<SlimCrawPeak>(in_peaks);
            //TODO -- switch to using priority queues - not enough time to test at the moment
            in_peaks.Sort(crawutils.CompByStartRTIdx);
            peaks_by_height.Sort((x, y) => x.peak_height.CompareTo(y.peak_height));

            peak_lists plist_struct = new peak_lists(peaks_by_height, in_peaks);

            //delete the rts iterators
            List<SlimCrawPeak> peaks_to_clear = new List<SlimCrawPeak>();
            List<int> rts_its = new List<int>();

            List<SlimCrawPeak> out_peaks = new List<SlimCrawPeak>();
            while (plist_struct.peaks_by_height.Count > 1)
            {
                SlimCrawPeak highest_p_ptr = plist_struct.peaks_by_height[plist_struct.peaks_by_height.Count - 1];
                int lookup = crawutils.lower_bound(plist_struct.peaks_by_rt, highest_p_ptr, crawutils.CompByStartRTIdx);
                if (lookup == plist_struct.peaks_by_rt.Count)
                    throw new Exception("should always find an exact RT peak match");

                SlimCrawPeak candidate = plist_struct.peaks_by_rt[lookup];
                if (!(candidate.start_rt_idx == highest_p_ptr.start_rt_idx))
                    throw new Exception("should have always find an exact match");
                bool comp_lh, comp_rh;
                bool merge_lh, merge_rh;
                comp_lh = comp_rh = true;
                merge_lh = merge_rh = false;

                if (lookup == 0)
                    comp_lh = false;
                if (lookup == (plist_struct.peaks_by_rt.Count - 1))
                    comp_rh = false;
                ///you can merge with up to two adjacent peaks
                rts_its.Add(lookup);
                peaks_to_clear.Add(plist_struct.peaks_by_rt[lookup]);
                if (comp_lh)
                {
                    SlimCrawPeak leftCandidate = plist_struct.peaks_by_rt[lookup - 1];
                    if (susceptible_to_merge(leftCandidate, candidate))
                    {
                        merge_lh = true;
                        rts_its.Add(lookup - 1);
                        peaks_to_clear.Add(leftCandidate);
                    }
                }
                if (comp_rh)
                {
                    SlimCrawPeak rightCandidate = plist_struct.peaks_by_rt[lookup + 1];
                    if (susceptible_to_merge(candidate, rightCandidate))
                    {
                        merge_rh = true;
                        rts_its.Add(lookup + 1);
                        peaks_to_clear.Add(rightCandidate);
                    }
                }
                if (merge_lh)
                    glom_peak(candidate, plist_struct.peaks_by_rt[lookup - 1]);
                if (merge_rh)
                    glom_peak(candidate, plist_struct.peaks_by_rt[lookup + 1]);

                //TODO inefficient..., we need to copy the peak after it is glommed, otherwise it is out of order in these cases...
                plist_struct.resort_peaks();
                //delete this peak from both vectors
                //delete each merged_peak from both vectors
                out_peaks.Add(candidate);
                foreach (SlimCrawPeak peak in peaks_to_clear)
                    plist_struct.clear_peak(peak);

                rts_its.Clear();
                peaks_to_clear.Clear();
            }

            if (plist_struct.peaks_by_height.Count > 0)
                out_peaks.Add(plist_struct.peaks_by_height[0]);
            return out_peaks;
        }

        ///merge peaks if lhs heuristics about slopes apply
        bool susceptible_to_merge(SlimCrawPeak lhs, SlimCrawPeak rhs)
        {
            if (Math.Abs(lhs.stop_rt_idx - rhs.start_rt_idx) > 1)
                return false;

            float max_height = Math.Max(lhs.peak_height, rhs.peak_height);
            float max_len = Math.Max(lhs.len, rhs.len);
            float one_peak_slope_constraint = (pf.method.one_peak_slope_merge_constraint * max_height) / max_len;
            //TODO calcs for mean_peak_slope_constraint in the line below do not make sense
            float mean_peak_slope_constraint = (pf.method.mean_slope_merge_constraint * max_height) / max_len;
            //opposite slopes are not a valid assumption, when merging two shoulders next to each other

            if ((lhs.bgslope * rhs.bgslope) < 0)
            {
                //at least one peak must have a high slope...
                if (one_peak_slope_constraint != 0
                    && ((Math.Abs(lhs.bgslope) >= one_peak_slope_constraint)
                    || (Math.Abs(rhs.bgslope) >= one_peak_slope_constraint)))
                    return true;
                if (mean_peak_slope_constraint != 0 &&
                    ((lhs.bgslope + rhs.bgslope) / 2.0) >= mean_peak_slope_constraint)
                    return true;
            }
            /* HACK BELOW - not configured with options   */

            else if (lhs.raw_height < rhs.raw_height)
            {
                if (lhs.bgslope > 0)
                    if (lhs.height_norm_slope() > 1)
                        return true;
            }
            else if (lhs.raw_height > rhs.raw_height)
            {
                if (rhs.bgslope < 0)
                    if (rhs.height_norm_slope() < 1)
                        return true;
            }
            else
            {
                float peak_to_peak_slope = (rhs.peak_height - lhs.peak_height) / (rhs.peak_rt_idx - lhs.peak_rt_idx);
                if (rhs.height_norm_slope() > 1)
                {
                    if (peak_to_peak_slope * rhs.bgslope > 1)
                        return true;
                }
                else if (lhs.height_norm_slope() > 1)
                {
                    if (peak_to_peak_slope * lhs.bgslope > 1)
                        return true;
                }
            }
            return false;
        }

        ///not implemented yet -- plan is to check the similarity of peaks by a correlation after normalizing for size
        bool susceptible_to_merge_by_shape(SlimCrawPeak _1, SlimCrawPeak _2)
        {
            return false; //ie not done yet
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="start_idx"></param>
        /// <param name="stop_idx"></param>
        /// <returns>The length of the used scratch (stop_idx - start_idx + 1)</returns>
        int set_bg_scratch(int start_idx, int stop_idx)
        {
            int len = stop_idx - start_idx + 1;
            if (len > bg_scratch.Length)
                bg_scratch = new float[len];
            if (pf.method.background_estimation_method == BackgroundEstimationMethod.PEAK_BOUNDARY_ESTIMATE)
            {
                float bgslope = calculate_slope(start_idx, stop_idx);
                double bg_val = active_chrom[start_idx];
                for (int i = 0; i < len; i++)
                {
                    if (active_chrom[start_idx + i] < bg_val)
                        bg_scratch[i] = active_chrom[start_idx + i];
                    else
                        bg_scratch[i] = (float)bg_val;
                    bg_val += bgslope;
                }
            }
            else if (pf.method.background_estimation_method == BackgroundEstimationMethod.LOWER_BOUNDARY)
            {
                float lower_bound = Math.Min(active_chrom[start_idx], active_chrom[stop_idx]);
                for (int i = 0; i < len; i++)
                {
                    if (active_chrom[start_idx + i] < lower_bound)
                        bg_scratch[i] = active_chrom[start_idx + i];
                    else
                        bg_scratch[i] = lower_bound;
                }
            }
            else if (pf.method.background_estimation_method == BackgroundEstimationMethod.MEAN_BOUNDARY)
            {
                float mean_bound = (active_chrom[start_idx] + active_chrom[stop_idx]) / 2.0f;
                for (int i = 0; i < len; i++)
                    bg_scratch[i] = mean_bound;
            }
            else if (pf.method.background_estimation_method == BackgroundEstimationMethod.FIXED_BACKGROUND)
            {
                if (!fixed_background.HasValue)
                    throw new Exception("must set a fixed background");

                for (int i = 0; i < len; i++)
                    bg_scratch[i] = fixed_background.Value;
            }
            else
                throw new Exception("unknown peak background estimate");
            return len;
        }

        /*! get_all_areas - calculate raw,background,and background-subtracted areas under a chromatogram, using
          techniques for 'CrawPeaks' as of 03/09.
          Background (bg) is simply calculated as the area under the pointwise minimum of
            a) a line drawn from the raw values at start,stop or:
            b) the actual data (which may fall below the line in (a)
          Areas are calculated using trapezoids.

          /param start_idx - starting index to chromatogram
          /param stop_idx  - stop index to chromatogram
          /param peak_idx  - 'peak' index in chromatogram
          /param raw_area  - ref float for raw area under curve
          /param bg_area   - ref float for background area
          /param peak_area - ref float for raw - bg
          /param bg_sub_height - bg subtracted height
        */
        void get_all_areas(int start_idx, int stop_idx, int peak_idx, ref float raw_area, ref float bg_area, ref float bg_sub_area, ref float bg_sub_height, ref float raw_height)
        {
            raw_area = (float)get_area(start_idx, stop_idx);
            int len = set_bg_scratch(start_idx, stop_idx);
            bg_sub_height = active_chrom[peak_idx] - bg_scratch[peak_idx - start_idx];
            raw_height = active_chrom[peak_idx];
            bg_area = (float)crawutils.area_under_curve(bg_scratch, 0, len - 1);
            bg_sub_area = raw_area - bg_area;
        }

        public void set_peak_bg_subtracted_area(SlimCrawPeak peak)
        {
            get_all_areas(peak.start_rt_idx, peak.stop_rt_idx, peak.peak_rt_idx,
                ref peak.raw_area, ref peak.bg_area, ref peak.peak_area, ref peak.peak_height, ref peak.raw_height);
        }

        public void calc_fwhm(SlimCrawPeak peak)
        {
            float[] chrom = active_chrom;
            float lh_height = chrom[peak.start_rt_idx];
            float rh_height = chrom[peak.stop_rt_idx];
            float height = peak.raw_height - Math.Min(lh_height, rh_height);
            float half_max = (float)(peak.raw_height - (height / 2.0));
            int lh_pt = -1, rh_pt = -1;
            float lh_hm, rh_hm;
            for (int i = peak.start_rt_idx; i < peak.peak_rt_idx; i++)
            {
                if (chrom[i] <= half_max && chrom[i + 1] >= half_max)
                {
                    lh_pt = i;
                    break;
                }
            }
            for (int i = peak.peak_rt_idx; i < Math.Min(peak.stop_rt_idx, chrom.Length - 2); i++)
            {
                if (chrom[i] >= half_max && chrom[i + 1] <= half_max)
                {
                    rh_pt = i;
                    break;
                }
            }
            peak.fwhm_calculated_ok = (lh_pt > -1 && rh_pt > -1);

            if (lh_pt == -1)
            {
                lh_hm = peak.start_rt_idx;
            }
            else
            {
                float frac_delta = (half_max - chrom[lh_pt]) / (chrom[lh_pt + 1] - chrom[lh_pt]);
                lh_hm = lh_pt + frac_delta;
            }
            if (rh_pt == -1)
            {
                rh_hm = peak.stop_rt_idx;
            }
            else
            {
                float frac_delta = (chrom[rh_pt] - half_max) / (chrom[rh_pt] - chrom[rh_pt + 1]);
                rh_hm = rh_pt + frac_delta;
            }
            peak.fwhm = rh_hm - lh_hm;
        }

        public void calc_fwfpcnt(SlimCrawPeak peak)
        {
            float[] chrom = active_chrom;
            float lh_height = chrom[peak.start_rt_idx];
            float rh_height = chrom[peak.stop_rt_idx];
            float height = peak.raw_height - Math.Min(lh_height, rh_height);
            float fpct_max = (float)(peak.raw_height - (height / 100 * 95));
            int lh_pt = -1, rh_pt = -1;
            float lh_hm, rh_hm;
            for (int i = peak.start_rt_idx; i < peak.peak_rt_idx; i++)
            {
                if (chrom[i] <= fpct_max && chrom[i + 1] >= fpct_max)
                {
                    lh_pt = i;
                    break;
                }
            }
            for (int i = peak.peak_rt_idx; i < Math.Min(peak.stop_rt_idx, chrom.Length - 2); i++)
            {
                if (chrom[i] >= fpct_max && chrom[i + 1] <= fpct_max)
                {
                    rh_pt = i;
                    break;
                }
            }

            if (lh_pt == -1)
            {
                lh_hm = peak.start_rt_idx;
            }
            else
            {
                float frac_delta = (fpct_max - chrom[lh_pt]) / (chrom[lh_pt + 1] - chrom[lh_pt]);
                lh_hm = lh_pt + frac_delta;
            }
            if (rh_pt == -1)
            {
                rh_hm = peak.stop_rt_idx;
            }
            else
            {
                float frac_delta = (chrom[rh_pt] - fpct_max) / (chrom[rh_pt] - chrom[rh_pt + 1]);
                rh_hm = rh_pt + frac_delta;
            }
            peak.fwfpct = rh_hm - lh_hm;
        }
        public void refind_peak_peak(SlimCrawPeak peak)
        {
            peak.peak_rt_idx = get_peakloc_in_range(peak.start_rt_idx, peak.stop_rt_idx);
        }

        public void reannotate_peak(SlimCrawPeak peak, int start_idx, int stop_idx)
        {
            peak.start_rt_idx = start_idx;
            peak.stop_rt_idx = stop_idx;
            peak.peak_rt_idx = get_peakloc_in_range(start_idx, stop_idx);
            set_peak_bg_subtracted_area(peak);
            //update this to incorporate other functions
        }

        public void peak_tweak(SlimCrawPeak peak, CrawPeakMethod m)
        {
            if (m.extend_to_zero_crossing)
                extend_to_zero_crossing(peak, m.fraction_to_valley);
            else if (m.extend_peak_to_lower_bound)
                extend_to_lower_boundary(peak, m.extend_allowed_asymmetry);
        }

        void extend_to_zero_crossing(SlimCrawPeak peak, float perc_towards_valley)
        {
            ////advance until next chrom2d zero crossing from rh_boundary

            float c2d_intensity_at_rh_bound = pf.chrom_2d[peak.stop_rt_idx];
            float c2d_intensity_at_lh_bound = pf.chrom_2d[peak.start_rt_idx];
            float rh_cutoff = c2d_intensity_at_rh_bound * perc_towards_valley;
            float lh_cutoff = c2d_intensity_at_lh_bound * perc_towards_valley;

            int plus_cross_to_right_idx = crawutils.lower_bound(pf.plus_crosses, peak.stop_rt_idx);
            if (plus_cross_to_right_idx == pf.plus_crosses.Count - 1)
            {
                //TODO -- extend out to last point?
            }
            else
            {
                plus_cross_to_right_idx++;
                for (int i = peak.stop_rt_idx; i <= pf.plus_crosses[plus_cross_to_right_idx]; i++)
                    if (pf.chrom_2d[i] < rh_cutoff)
                        peak.stop_rt_idx = i - 1;
            }

            int minus_cross_to_left_idx = crawutils.lower_bound(pf.minus_crosses, peak.start_rt_idx);
            if (minus_cross_to_left_idx == 0)
            {
                //TODO -- extend out to last point?
            }
            else
            {
                for (int i = peak.start_rt_idx; i >= pf.minus_crosses[plus_cross_to_right_idx]; i--)
                    if (pf.chrom_2d[i] < lh_cutoff)
                        peak.stop_rt_idx = i + 1;
            }
        }

        void extend_to_lower_boundary(SlimCrawPeak peak, float allowed_asymmetry = 1.0f)
        {
            //1. find lower point
            //2. find farther end-point
            //3. extend higher point until either: lower point is reached,

            //get intensity at either edge...
            float intensity_at_lh_bound = pf.chrom[peak.start_rt_idx];
            float intensity_at_rh_bound = pf.chrom[peak.stop_rt_idx];

            int max_distance = Math.Max(peak.peak_rt_idx - peak.start_rt_idx, peak.stop_rt_idx - peak.peak_rt_idx);
            int min_distance = Math.Min(peak.peak_rt_idx - peak.start_rt_idx, peak.stop_rt_idx - peak.peak_rt_idx);
            int max_allowed_distance = (int)Math.Round(max_distance * allowed_asymmetry);
            int max_delta_from_edge = max_allowed_distance - min_distance;

            if (intensity_at_lh_bound < intensity_at_rh_bound)
            {
                //extend bound rightwards until either max distance is reached, or
                int delta = 0;
                for (; delta < max_delta_from_edge; delta++)
                {
                    int pos = peak.stop_rt_idx + delta;
                    if (pf.chrom[pos] <= intensity_at_lh_bound)
                        break;
                }
                peak.stop_rt_idx += delta; //changed peak boundary, what would I need to adjust now...?
            }

            if (intensity_at_rh_bound < intensity_at_lh_bound)
            {
                //extend bound leftwards until either max distance is reached, or
                int delta = 0;
                for (; delta < max_delta_from_edge; delta++)
                {
                    int pos = peak.start_rt_idx - delta;
                    if (pf.chrom[pos] <= intensity_at_rh_bound)
                        break;
                }
                peak.start_rt_idx -= delta; //changed peak boundary, what would I need to adjust now...?
            }
        }

        void ratchet_back_to_frac_maxval(SlimCrawPeak peak, float frac = 0.01f, float bg_level = 0.0f)
        {
            float[] chrom = active_chrom;
            if (frac > 1)
                throw new Exception("define the fraction ranging from 0 to 1");

            float min_level = (peak.raw_height - bg_level) * frac;
            for (int i = peak.start_rt_idx; i < peak.peak_rt_idx; i++)
            {
                if (chrom[i] >= min_level)
                {
                    peak.start_rt_idx = i;
                    break;
                }
            }
            for (int i = peak.stop_rt_idx; i > peak.peak_rt_idx; i--)
            {
                if (chrom[i] >= min_level)
                {
                    peak.stop_rt_idx = i;
                    break;
                }
            }
        }

        void extend_to_1d_zero(SlimCrawPeak peak, bool start_at_peak)
        {
            int delta_to_lh, delta_to_rh;
            int lh_start, rh_start;
            if (start_at_peak)
                lh_start = rh_start = peak.peak_rt_idx;
            else
            {
                lh_start = peak.start_rt_idx;
                rh_start = peak.stop_rt_idx;
            }

            //find crossing point in the 1st derivative
            (List<int> c1d_plus_crosses, List<int> c1d_minus_crosses) = pf.find_cross_points(pf.chrom_1d);

            /*ER BUGBUG
              PLACE WHERE PEAKS ARE EXTENDED USING 1st derivative.
              crossing points in chrom_1d are found above.
            so you would need to limit the 'walk' out to the 1st derivative crossing point.
            The other option is to investigate what I was doing in the function below, local minimum.
            (perhaps walk out to the point where either the 1st deriv hits zero or the intensity falls below some factor of the peak
             intensity, like 0.05?.

            -- it may be best later on, after peak_tweak or whatever the heck it is is called, to trim peaks back, rather than rejecting them... hmmm..
        */
            //extend left hand
            if (!(lh_start < c1d_plus_crosses[0]))
            {
                int lh_it = crawutils.lower_bound(c1d_plus_crosses, lh_start);
                lh_it--;
                int leftmost_c1d_plus_cross = c1d_plus_crosses[lh_it];
                //since the cross is defined as going left-to-right, nudge to the other side
                leftmost_c1d_plus_cross++;
                delta_to_lh = lh_start - leftmost_c1d_plus_cross;
            }
            else
            {
                delta_to_lh = 0;
            }

            //extend right hand
            if (!(rh_start > c1d_plus_crosses.Last()))
            {
                int rh_it = crawutils.lower_bound(c1d_plus_crosses, rh_start);
                if (rh_it == c1d_plus_crosses.Count - 1)
                {

                }
                else
                {
                    //rh_it++;
                }
                int rightmost_c1d_plus_cross = c1d_plus_crosses[rh_it];
                delta_to_rh = rightmost_c1d_plus_cross - rh_start;
            }
            else
            {
                delta_to_rh = 0;
            }

            //for now, extend boundaries by delta_to_lh, delta_to_rh - note that we can also use those to enforce some peak asymmetry

            if (start_at_peak)
            {
                peak.start_rt_idx = peak.peak_rt_idx - delta_to_lh;
                peak.stop_rt_idx = peak.peak_rt_idx + delta_to_rh;
            }
            else
            {
                peak.start_rt_idx -= delta_to_lh;
                peak.stop_rt_idx += delta_to_rh;
            }
        }

        void extend_to_1d_zero_local_minimum(SlimCrawPeak peak, bool start_at_peak)
        {
            int delta_to_lh, delta_to_rh;
            int lh_start, rh_start;
            if (start_at_peak)
                lh_start = rh_start = peak.peak_rt_idx;
            else
            {
                lh_start = peak.start_rt_idx;
                rh_start = peak.stop_rt_idx;
            }

            //find crossing point in the 1st derivative
            (List<int> c1d_plus_crosses, List<int> c1d_minus_crosses) = pf.find_cross_points(pf.chrom_1d);

            //extend left hand

            int peak_intensity = (int)active_chrom[peak.peak_rt_idx];
            int last_minimum_intensity = peak_intensity;

            if (!(lh_start < c1d_plus_crosses[0]))
            {
                int lh_it = crawutils.lower_bound(c1d_plus_crosses, lh_start);
                if (lh_it == c1d_plus_crosses.Count)
                    throw new Exception("makes no sense");

                //walk leftwards
                do
                {
                    if (lh_it == 0)
                        break;
                    int this_valley_intensity = (int)active_chrom[c1d_plus_crosses[lh_it]];
                    if (this_valley_intensity > last_minimum_intensity)
                    {
                        //the last valley was the local minimum
                        lh_it++;
                        break;
                    }
                    last_minimum_intensity = this_valley_intensity;
                    lh_it--;
                } while (lh_it != 0);

                int leftmost_c1d_plus_cross = c1d_plus_crosses[lh_it];
                //since the cross is defined as going left-to-right, nudge to the other side
                leftmost_c1d_plus_cross++;
                delta_to_lh = lh_start - leftmost_c1d_plus_cross;
            }
            else
            {
                delta_to_lh = 0;
            }

            //extend right hand
            last_minimum_intensity = peak_intensity;

            if (!(rh_start > c1d_minus_crosses.Last()))
            {
                int rh_it = crawutils.lower_bound(c1d_minus_crosses, rh_start);

                do
                {
                    if (rh_it == c1d_minus_crosses.Count)
                    {
                        rh_it--;
                        break;
                    }
                    int this_valley_intensity = (int)active_chrom[c1d_minus_crosses[rh_it]];
                    if (this_valley_intensity > last_minimum_intensity)
                    {
                        //the last valley was the local minimum
                        rh_it--;
                        break;
                    }
                    last_minimum_intensity = this_valley_intensity;
                    rh_it++;
                } while (rh_it != c1d_plus_crosses.Count);

                int rightmost_c1d_minus_cross = c1d_minus_crosses[rh_it];
                delta_to_rh = rightmost_c1d_minus_cross - rh_start;
            }
            else
            {
                delta_to_rh = 0;
            }

            //for now, extend boundaries by delta_to_lh, delta_to_rh - note that we can also use those to enforce some peak asymmetry

            if (start_at_peak)
            {
                peak.start_rt_idx = peak.peak_rt_idx - delta_to_lh;
                peak.stop_rt_idx = peak.peak_rt_idx + delta_to_rh;
            }
            else
            {
                peak.start_rt_idx -= delta_to_lh;
                peak.stop_rt_idx += delta_to_rh;
            }
        }

        ///this takes the set of peaks and extends them, based on potential shoulder peaks which have odd alignments
        void extend_peak_set(IList<SlimCrawPeak> in_peaks, IList<SlimCrawPeak> out_peaks, bool start_at_peak = false)
        {
            if (in_peaks.Count == 0)
                return;

            List<SlimCrawPeak> peak_ptrs_by_rt = new List<SlimCrawPeak>(in_peaks);

            peak_ptrs_by_rt.Sort(crawutils.CompByPeakRTIdx);
            List<PeaksUsedOrNot> peaks_by_rt = new List<PeaksUsedOrNot>(peak_ptrs_by_rt.Select(p => new PeaksUsedOrNot(p)));
            List<PeaksUsedOrNot> peaks_by_I = new List<PeaksUsedOrNot>(peaks_by_rt);
            peaks_by_I.Sort(crawutils.CompByPeakRawHeight);
            PeakCoordSet peak_coords = new PeakCoordSet(peak_ptrs_by_rt);
            int intensity_index = peaks_by_I.Count - 1;
            do
            {
                intensity_index--;
                PeaksUsedOrNot here = peaks_by_I[intensity_index];
                if (!here.used)
                {
                    //extend peak
                    extend_to_1d_zero(here.p, pf.method.extend_from_peak_rt);
                    if (pf.method.ratchet_back_to_frac_maxval > -1)
                        ratchet_back_to_frac_maxval(here.p, pf.method.ratchet_back_to_frac_maxval);
                    //exclude overlappers -- find peaks from my start to my end
                    //TODO -- deal with overlappers... or convert to linear walk algorithm
                    //safe, since push_back will copy the dereferenced SlimCrawPeakPtr
                    if (pf.method.exclude_extension_overlaps_by_peakrt)
                        mark_overlaps_by_peak_rt(peaks_by_rt, peak_coords, here);

                    out_peaks.Add(here.p);
                    here.used = true;

                }
            } while (intensity_index > 0);
        }

        //assumes peaks are sorted by peakrt
        void mark_overlaps_by_peak_rt(IList<PeaksUsedOrNot> peaks, PeakCoordSet peak_set, PeaksUsedOrNot p)
        {
            var overlap_bounds = peak_set.find_overlap_bounds_by_peak_rt_idx(p.p);
            for (int i = overlap_bounds.leftmost; i <= overlap_bounds.rightmost; i++)
                peaks[i].used = true;
        }

        double get_raw_area(SlimCrawPeak peak) => get_area(peak.start_rt_idx, peak.stop_rt_idx);

        double get_raw_area(int start_rt_idx, int stop_rt_idx) => get_area(start_rt_idx, stop_rt_idx);

        public void set_peak_slope(SlimCrawPeak peak)
        {
            float slope = calculate_slope(peak.start_rt_idx, peak.stop_rt_idx);
            peak.bgslope = slope;
        }

        float calculate_slope(int start_rt_idx, int stop_rt_idx)
        {
            int len = calc_len(start_rt_idx, stop_rt_idx);
            float delta = active_chrom[stop_rt_idx] - active_chrom[start_rt_idx];
            float d = len - 1;
            return delta / d;

        }

        double get_area(int start_rt_idx, int stop_rt_idx)
        {
            return crawutils.area_under_curve(active_chrom, start_rt_idx, stop_rt_idx);
        }
    }
}