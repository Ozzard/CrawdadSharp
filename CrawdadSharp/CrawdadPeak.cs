using System;
using System.Collections.Generic;

namespace CrawdadSharp
{
    public class CrawdadPeak
    {
        public CrawdadPeak(SlimCrawPeak crawPeak)
        {
            TimeIndex = crawPeak.peak_rt_idx;
            StartIndex = crawPeak.start_rt_idx;
            EndIndex = crawPeak.stop_rt_idx;
            Height = crawPeak.peak_height;
            // Note: Crawdad had a bug that caused it to return negative areas.
            //       Using Math::Max() below left in place to protect against that,
            //       even though the bug is believed to be fixed.
            Area = Math.Max(0.0f, crawPeak.peak_area);
            BackgroundArea = Math.Max(0.0f, crawPeak.bg_area);
            Fwhm = crawPeak.fwhm;
            FwhmDegenerate = !crawPeak.fwhm_calculated_ok;
            Fwfpct = crawPeak.fwfpct;
        }

        public int TimeIndex { get; }
        public int StartIndex { get; set; }
        public int EndIndex { get; set; }
        public int Length => EndIndex - StartIndex + 1;
        public int Center => (int)(((StartIndex + EndIndex) / 2.0) + 0.5);
        public float Area { get; }
        public float BackgroundArea { get; }
        public float Height { get; }
        public float Fwhm { get; }
        public float Fwfpct { get; }
        public bool FwhmDegenerate { get; }

        public bool IsIdentified(IEnumerable<int> idIndices)
        {
            foreach (int idIndex in idIndices)
                if (StartIndex <= idIndex && idIndex <= EndIndex)
                    return true;
            return false;
        }

        public override string ToString() => $"a = {Area}, bg = {BackgroundArea}, s = {StartIndex}, e = {EndIndex}, r = {TimeIndex}";
    };
};
