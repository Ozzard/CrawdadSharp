using CrawdadSharp;
using NUnit.Framework;
using System.Collections.Generic;

namespace TestCrawdadSharp
{
    public class Tests
    {
        [SetUp]
        public void Setup()
        {
        }

        /// <summary>
        /// Same test as Hannes Roest uses in his fork of the Crawdad code - 1 of 2.
        /// </summary>
        [Test]
        public void TestData1Of2()
        {
            float[] rts = { 1474.34f, 1477.11f, 1479.88f, 1482.64f, 1485.41f, 1488.19f, 1490.95f, 1493.72f, 1496.48f, 1499.25f, 1502.03f, 1504.8f, 1507.56f, 1510.33f, 1513.09f, 1515.87f, 1518.64f, 1521.42f };
            float[] intensities = { 3.26958f, 3.74189f, 3.31075f, 86.1901f, 3.47528f, 387.864f, 13281f, 6375.84f, 39852.6f, 2.66726f, 612.747f, 3.34313f, 793.12f, 3.29156f, 4.00586f, 4.1591f, 3.23035f, 3.90591f };

            CrawdadPeakFinder pf = new CrawdadPeakFinder();
            pf.SetChromatogram(rts, intensities);
            IList<CrawdadPeak> res = pf.CalcPeaks();

            Assert.IsTrue(res.Count == 1);
            Assert.IsTrue(res[0].StartIndex == 2);
            Assert.IsTrue(res[0].TimeIndex == 8);
            Assert.IsTrue(res[0].EndIndex == 13);
        }

        /// <summary>
        /// Same test as Hannes Roest uses in his fork of the Crawdad code - 2 of 2.
        /// </summary>
        [Test]
        public void TestData2Of2()
        {
            float[] rts = { 1473.55f, 1476.31f, 1479.08f, 1481.84f, 1484.61f, 1487.39f, 1490.15f, 1492.92f, 1495.69f, 1498.45f, 1501.23f, 1504f, 1506.76f, 1509.53f, 1512.29f, 1515.07f, 1517.84f, 1520.62f };
            float[] intensities = { 3.44054f, 2142.31f, 3.58763f, 3076.97f, 6663.55f, 45681f, 157694f, 122844f, 86034.7f, 85391.1f, 15992.8f, 2293.94f, 6934.85f, 2735.18f, 459.413f, 3.93863f, 3.36564f, 3.44005f };

            CrawdadPeakFinder pf = new CrawdadPeakFinder();
            pf.SetChromatogram(rts, intensities);
            IList<CrawdadPeak> res = pf.CalcPeaks();

            Assert.IsTrue(res.Count == 1);
            Assert.IsTrue(res[0].StartIndex == 2);
            Assert.IsTrue(res[0].TimeIndex == 6);
            Assert.IsTrue(res[0].EndIndex == 13);
        }
    }
}