using System;
using System.Collections.Generic;
using System.Text;

namespace SA
{
    public class Helper
    {
        public static string ArrayToString(List<int> array)
        {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < array.Count; i++)
            {
                if (i > 0)
                    sb.Append(" ");

                sb.Append(array[i]);
            }

            return sb.ToString();
        }
    }
}
