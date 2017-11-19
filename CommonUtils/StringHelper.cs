using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CommonUtils
{
    public class StringHelper
    {
        public static string NumberToString(int number)
        {
            if (number >= 1000 * 1000)
                return "" + number / (1000 * 1000) + "M";
            if (number >= 1000)
                return "" + number / 1000 + "K";
            return "" + number;
        }
    }
}
