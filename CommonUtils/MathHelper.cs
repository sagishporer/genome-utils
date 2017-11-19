using System;
using System.Collections.Generic;
using System.Text;

namespace SA
{
    public static class MathHelper
    {
        public static double Avg(List<int> d)
        {
            double sum = 0;
            for (int i = 0; i < d.Count; i++)
                sum += (double)d[i];

            return sum / (double)d.Count;
        }

        public static double Avg(List<double> d)
        {
            double sum = 0;
            for (int i = 0; i < d.Count; i++)
                sum += d[i];

            return sum / (double)d.Count;
        }

        public static double Var(List<double> d)
        {
            double sum = 0;
            double sumSqr = 0;

            for (int i = 0; i < d.Count; i++)
            {
                sum += d[i];
                sumSqr += (d[i] * d[i]);
            }

            double avg = sum / (double)d.Count;
            return (sumSqr / (double)d.Count) - avg * avg;
        }

        public static int Max(List<int> d)
        {
            int max = int.MinValue;
            for (int i = 0; i < d.Count; i++)
                max = Math.Max(max, d[i]);

            return max;
        }

        public static double Max(List<double> d)
        {
            double max = double.MinValue;
            for (int i = 0; i < d.Count; i++)
                max = Math.Max(max, d[i]);

            return max;
        }

        public static double Min(List<double> d)
        {
            double min = double.MaxValue;
            for (int i = 0; i < d.Count; i++)
                min = Math.Min(min, d[i]);

            return min;
        }

        /// <summary>
        /// L2 = Sqrt(Sigma(l^2))
        /// </summary>
        /// <param name="d"></param>
        /// <returns></returns>
        public static double L2(List<double> d)
        {
            double sumSqr = 0;
            for (int i = 0; i < d.Count; i++)
                sumSqr += (d[i] * d[i]);

            return Math.Sqrt(sumSqr);
        }
    }
}
