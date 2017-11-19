using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CountDNABases
{
    public class StatBuilder
    {
        private bool mBuildHistogram;
        private int mHistogramBuckets;

        private int[] mHistogram;
        private int mHistogramLength;

        private String mTitle;

        private ulong mCount;
        private double mSum;
        private double mSquareSum;
        private double mMinValue;
        private double mMaxValue;

        private long mCountOneZero;
        private long mCountBothZero;

        public string Title { get { return mTitle; } }

        public long CountOneZero { get { return mCountOneZero; } }
        public long CountBothZero { get { return mCountBothZero; } }

        public StatBuilder(string title, bool buildHistogram, int histogramBuckets)
        {
            mBuildHistogram = buildHistogram;
            if (mBuildHistogram == true)
            {
                mHistogramBuckets = histogramBuckets;
                mHistogram = new int[mHistogramBuckets];
                mHistogramLength = mHistogram.Length;
                Array.Clear(mHistogram, 0, mHistogramLength);
            }

            mTitle = title;

            mCount = 0;
            mSum = 0;
            mSquareSum = 0;
            mMinValue = double.MaxValue;
            mMaxValue = double.MinValue;

            mCountOneZero = 0;
            mCountBothZero = 0;
        }

        public void IncreaseZeroCounts(double a, double b)
        {
            if ((a == 0) || (b == 0))
                mCountOneZero++;
            if ((a == 0) && (b == 0))
                mCountBothZero++;
        }

        public void IncreaseZeroCounts(int a, int b)
        {
            if ((a == 0) || (b == 0))
                mCountOneZero++;
            if ((a == 0) && (b == 0))
                mCountBothZero++;
        }

        public void AddValue(double val)
        {
            mCount++;
            mSum += val;
            mSquareSum += (val * val);
            //mMinValue = Math.Min(val, mMinValue);
            if (val < mMinValue)
                mMinValue = val;
            //mMaxValue = Math.Max(val, mMaxValue);
            if (val > mMaxValue)
                mMaxValue = val;

            if (mBuildHistogram == true)
            {
                // Build histogram, assume source 0..1 ==> 0..HISTOGRAM_BUCKETS
                int histogramBucket = (int)(val * (double)mHistogramBuckets);
                if (histogramBucket < 0)
                    histogramBucket = 0;
                if (histogramBucket >= mHistogramLength)
                    histogramBucket = mHistogramLength - 1;
                mHistogram[histogramBucket]++;
            }
        }

        public int[] GetHistogram()
        {
            if (mBuildHistogram == false)
                throw new Exception("mBuildHistogram == false");

            return mHistogram;
        }

        public double GetMaxValue() { return mMaxValue; }
        public double GetMinValue() { return mMinValue; }

        public ulong GetCount()
        {
            return mCount;
        }

        public double GetSum()
        {
            return mSum;
        }

        public double GetSquareSum()
        {
            return mSquareSum;
        }

        public double GetAverage()
        {
            return mSum / (double)mCount;
        }

        public double GetVariance()
        {
            double avg = GetAverage();

            return (mSquareSum / (double)mCount) - (avg * avg);
        }
    }
}
