using System.Collections.Generic;

namespace CrawdadSharp
{
    public class StackCrawPeakFinder : BaseCrawPeakFinder
    {
        public List<SlimCrawPeak> sps = new List<SlimCrawPeak>();

        protected override SlimCrawPeak get_peak_ptr(int idx) => sps[idx];

        public StackCrawPeakFinder()
            : base()
        {
        }

        public StackCrawPeakFinder(CrawPeakMethod m)
            : base(m)
        {
        }

        protected override void peak_voodoo(int lh_valley, int rh_valley, int peak_loc)
        {
            SlimCrawPeak peak = new SlimCrawPeak(lh_valley, rh_valley, peak_loc, chrom, chrom_bg_scratch, mz_idx);

            annotator.refind_peak_peak(peak);
            annotator.peak_tweak(peak, method);

            //peak.peak_height = 0.0f;
            //note this needs to altered if we do not use the standard of defining the peak area as
            //stretching from one end to the other
            annotator.set_peak_slope(peak);
            annotator.set_peak_bg_subtracted_area(peak);
            annotator.calc_fwhm(peak);
            sps.Add(peak);
        }

        protected override int get_num_stored_peaks() => sps.Count;

        protected override void clear_sps()
        {
            sps.Clear();
        }
    }
}
