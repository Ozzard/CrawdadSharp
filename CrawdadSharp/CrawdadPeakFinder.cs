using System;
using System.Collections.Generic;

namespace CrawdadSharp
{
    public class CrawdadPeakFinder
    {
        public CrawdadPeakFinder()
        {
            _pPeakFinder = new StackCrawPeakFinder
            {
                slim = true
            };
            _pPeakFinder.method.peak_location_meth = PeakPeakMethod.MAXIMUM_PEAK;
            _pPeakFinder.method.background_estimation_method = BackgroundEstimationMethod.LOWER_BOUNDARY;
        }

        public float FullWidthHalfMax
        {
            get { return _pPeakFinder.method.Fwhm; }
            set
            {
                _pPeakFinder.method.Fwhm = value * 3;
                _pPeakFinder.method.min_len = (int)(value / 4.0 + 0.5);
            }
        }

        public float StdDev
        {
            get { return _pPeakFinder.method.Sd; }
            set { _pPeakFinder.method.Sd = value; }
        }

        public List<CrawdadPeak> CalcPeaks() { return CalcPeaks(-1); }
        public List<CrawdadPeak> CalcPeaks(int max) { return CalcPeaks(max, new int[0]); }



        /// Padding data added before and after real data to ensure peaks near the edges get detected.
        private int _widthDataWings;

        private StackCrawPeakFinder _pPeakFinder;
        public void SetChromatogram(IList<double> times, IList<double> intensities)
        {
            // TODO: Check times to make sure they are evenly spaced

            // Marshall intensities to vector for Crawdad
            int len = intensities.Count;
            float[] intensitiesCrawdad = new float[len];
            double baselineIntensity = double.MaxValue;
            double maxIntensity = 0;
            int maxIntensityIndex = -1;
            for (int i = 0; i < len; i++)
            {
                float intensity = (float)intensities[i];
                intensitiesCrawdad[i] = intensity;

                // Keep track of where the maximum intensity is
                if (intensity > maxIntensity)
                {
                    maxIntensity = intensity;
                    maxIntensityIndex = i;
                }
                if (intensity < baselineIntensity)
                    baselineIntensity = intensity;
            }
            if (baselineIntensity == double.MaxValue)
                baselineIntensity = 0;

            SetChromatogram(intensitiesCrawdad, maxIntensityIndex, baselineIntensity);
        }

        public void SetChromatogram(IList<float> times, IList<float> intensities)
        {
            // TODO: Check times to make sure they are evenly spaced

            // Marshall intensities to vector for Crawdad
            int len = intensities.Count;
            float[] intensitiesCrawdad = new float[len];
            double baselineIntensity = double.MaxValue;
            double maxIntensity = 0;
            int maxIntensityIndex = -1;
            for (int i = 0; i < len; i++)
            {
                float intensity = intensities[i];
                intensitiesCrawdad[i] = intensity;

                // Keep track of where the maximum intensity is
                if (intensity > maxIntensity)
                {
                    maxIntensity = intensity;
                    maxIntensityIndex = i;
                }
                if (intensity < baselineIntensity)
                    baselineIntensity = intensity;
            }
            if (baselineIntensity == double.MaxValue)
                baselineIntensity = 0;

            SetChromatogram(intensitiesCrawdad, maxIntensityIndex, baselineIntensity);
        }

        void SetChromatogram(float[] intensities, int maxIntensityIndex, double baselineIntensity)
        {
            // Find the peak width of the maximum intensity point at half its height.
            int fwhm = 6;
            if (maxIntensityIndex != -1)
            {
                double halfHeight = (intensities[maxIntensityIndex] - baselineIntensity) / 2 + baselineIntensity;
                int iStart = 0;
                for (int i = maxIntensityIndex - 1; i >= 0; i--)
                {
                    if (intensities[i] < halfHeight)
                    {
                        iStart = i;
                        break;
                    }
                }
                int len = intensities.Length;
                int iEnd = len - 1;
                for (int i = maxIntensityIndex + 1; i < len; i++)
                {
                    if (intensities[i] < halfHeight)
                    {
                        iEnd = i;
                        break;
                    }
                }
                fwhm = Math.Max(fwhm, iEnd - iStart);
            }
            FullWidthHalfMax = fwhm;

            int fwfpcnt = 0;
            if (maxIntensityIndex != -1)
            {
                double fpcntHeight = (intensities[maxIntensityIndex] - baselineIntensity) / 20 + baselineIntensity;
                int iBeginning = 0;
                for (int i = maxIntensityIndex - 1; i >= 0; i--)
                {
                    if (intensities[i] < fpcntHeight)
                    {
                        iBeginning = i;
                        break;
                    }
                }
                int len = intensities.Length;
                int iTail = len - 1;
                for (int i = maxIntensityIndex + 1; i < len; i++)
                {
                    if (intensities[i] < fpcntHeight)
                    {
                        iTail = i;
                        break;
                    }
                }
                fwfpcnt = Math.Max(fwfpcnt, iTail - iBeginning);
            }

                _widthDataWings = (int)(FullWidthHalfMax * 2);

            if (_widthDataWings > 0)
            {
                float[] intensitiesWithWings = new float[intensities.Length + 2 * _widthDataWings];
                for (int i = 0; i < _widthDataWings; i++)
                {
                    intensitiesWithWings[i] = (float)baselineIntensity;
                    intensitiesWithWings[i + _widthDataWings + intensities.Length] = (float)baselineIntensity;
                }
                Array.Copy(intensities, 0, intensitiesWithWings, _widthDataWings, intensities.Length);
                intensities = intensitiesWithWings;
            }

            _pPeakFinder.clear();
            _pPeakFinder.set_chrom(intensities, 0);
        }

        public IList<float> Intensities2d
        {
            get
            {
                // Make sure the 2d chromatogram is populated
                if (_pPeakFinder.chrom_2d.Length != _pPeakFinder.chrom.Length)
                {
                    _pPeakFinder.chrom_2d = new float[_pPeakFinder.chrom.Length];
                    _pPeakFinder.get_2d_chrom(_pPeakFinder.chrom, _pPeakFinder.chrom_2d);
                }

                // Marshal 2nd derivative peaks back to managed list
                float[] intensities2d = new float[_pPeakFinder.chrom.Length - _widthDataWings * 2];
                Array.Copy(_pPeakFinder.chrom_2d, _widthDataWings, intensities2d, 0, intensities2d.Length);
                return intensities2d;
            }
        }

        public IList<float> Intensities1d
        {
            get
            {
                // Make sure the 2d chromatogram is populated
                if (_pPeakFinder.chrom_1d.Length != _pPeakFinder.chrom.Length)
                {
                    _pPeakFinder.chrom_1d = new float[_pPeakFinder.chrom.Length];
                    _pPeakFinder.get_1d_chrom(_pPeakFinder.chrom, _pPeakFinder.chrom_1d);
                }

                // Marshal 1st derivative peaks back to managed list
                float[] intensities1d = new float[_pPeakFinder.chrom.Length - _widthDataWings * 2];
                Array.Copy(_pPeakFinder.chrom_1d, _widthDataWings, intensities1d, 0, intensities1d.Length);
                return intensities1d;
            }
        }

        CrawdadPeak GetPeak(int startIndex, int endIndex)
        {
            startIndex += _widthDataWings;
            endIndex += _widthDataWings;

            SlimCrawPeak peak = new SlimCrawPeak();
            _pPeakFinder.annotator.reannotate_peak(peak, startIndex, endIndex);
            _pPeakFinder.annotator.calc_fwhm(peak);

            peak.start_rt_idx -= _widthDataWings;
            peak.stop_rt_idx -= _widthDataWings;
            peak.peak_rt_idx -= _widthDataWings;

            return new CrawdadPeak(peak);
        }

        List<CrawdadPeak> CalcPeaks(int max, int[] idIndices)
        {
            // Find peaks
            _pPeakFinder.call_peaks();

            // Marshall found peaks to managed list
            List<CrawdadPeak> result = new List<CrawdadPeak>((int)_pPeakFinder.sps.Count);
            double totalArea = 0;
            int stop_rt = _pPeakFinder.chrom.Length - _widthDataWings - 1;
            int adjust_stop_rt = stop_rt - _widthDataWings;
            foreach (SlimCrawPeak itPeak in _pPeakFinder.sps)
            {
                if (itPeak.start_rt_idx < stop_rt && itPeak.stop_rt_idx > _widthDataWings)
                {
                    double rheight = itPeak.peak_height / itPeak.raw_height;
                    double rarea = itPeak.peak_area / itPeak.raw_area;

                    if (rheight > 0.02 && rarea > 0.02)
                    {
                        itPeak.start_rt_idx = Math.Max(_widthDataWings, itPeak.start_rt_idx);
                        itPeak.start_rt_idx -= _widthDataWings;
                        itPeak.peak_rt_idx = Math.Max(_widthDataWings, Math.Min(stop_rt, itPeak.peak_rt_idx));
                        itPeak.peak_rt_idx -= _widthDataWings;
                        itPeak.stop_rt_idx = Math.Max(_widthDataWings, Math.Min(stop_rt, itPeak.stop_rt_idx));
                        itPeak.stop_rt_idx -= _widthDataWings;

                        result.Add(new CrawdadPeak(itPeak));

                        totalArea += itPeak.peak_area;
                    }
                }
            }

            // If max is not -1, then return the max most intense peaks, plus any
            // peaks that have been identified with MS/MS peptide search results
            if (max != -1)
            {
                // Shorten the list before performing the slow sort by intensity.
                // The sort shows up as bottleneck in a profiler.
                int lenResult = result.Count;
                float intensityCutoff = 0;
                FindIntensityCutoff(result, 0, (float)(totalArea / lenResult) * 2, max, 1, ref intensityCutoff, ref lenResult);

                List<KeyValuePair<CrawdadPeak, bool>> resultFiltered =
                    new List<KeyValuePair<CrawdadPeak, bool>>(lenResult);
                for (int i = 0, lenOrig = result.Count; i < lenOrig; i++)
                {
                    CrawdadPeak peak = result[i];
                    bool isIdentified = peak.IsIdentified(idIndices);
                    if (isIdentified || peak.Area >= intensityCutoff || intensityCutoff == 0)
                        resultFiltered.Add(new KeyValuePair<CrawdadPeak, bool>(peak, isIdentified));
                }

                resultFiltered.Sort(new Comparison<KeyValuePair<CrawdadPeak, bool>>(OrderIdAreaDesc));
                if (max < resultFiltered.Count)

                    resultFiltered.RemoveRange(max, resultFiltered.Count - max);

                result = new List<CrawdadPeak>(resultFiltered.Count);
                foreach (KeyValuePair<CrawdadPeak, bool> peakId in resultFiltered)
                    result.Add(peakId.Key);
            }

            return result;
        }

        void FindIntensityCutoff(List<CrawdadPeak> listPeaks, float left, float right, int minPeaks, int calls, ref float cutoff, ref int len)
        {
            if (calls < 3)
            {
                float mid = (left + right) / 2;
                int count = FilterPeaks(listPeaks, mid);
                if (count < minPeaks)
                    FindIntensityCutoff(listPeaks, left, mid, minPeaks, calls + 1, ref cutoff, ref len);
                else
                {
                    cutoff = mid;
                    len = count;
                    if (count > minPeaks * 1.5)
                        FindIntensityCutoff(listPeaks, mid, right, minPeaks, calls + 1, ref cutoff, ref len);
                }
            }
        }

        int FilterPeaks(List<CrawdadPeak> listPeaks, float intensityCutoff)
        {
            int nonNoise = 0;
            foreach (CrawdadPeak peak in listPeaks)
            {
                if (peak.Area >= intensityCutoff)
                    nonNoise++;
            }
            return nonNoise;
        }

        int OrderIdAreaDesc(KeyValuePair<CrawdadPeak, bool> peakId1,
            KeyValuePair<CrawdadPeak, bool> peakId2)
        {
            // Identified peaks come first
            bool id1 = peakId1.Value, id2 = peakId2.Value;
            if (id1 != id2)
                return id1 ? -1 : 1;

            // Then area descending
            float a1 = peakId1.Key.Area, a2 = peakId2.Key.Area;
            return (a1 > a2 ? -1 : (a1 < a2 ? 1 : 0));
        }
    }
}
