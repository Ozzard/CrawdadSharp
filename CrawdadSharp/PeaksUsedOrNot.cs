namespace CrawdadSharp
{
    class PeaksUsedOrNot
    {
        public readonly SlimCrawPeak p;
        public bool used;

        public PeaksUsedOrNot(SlimCrawPeak in_p)
        {
            p = in_p;
            used = false;
        }

        public float peak_rt => (null == p) ? -1.0f : p.peak_rt_idx;
        public float peak_height => (null == p) ? -1.0f : p.peak_height;
        public float peak_raw_height => (null == p) ? -1.0f : p.raw_height;
    }
}