using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Security.Cryptography;

namespace SA
{
    public static class HashHelper
    {
        /// <summary>
        /// Code based on sample code from MSDN:
        /// http://msdn.microsoft.com/en-us/library/s02tk69a.aspx
        /// </summary>
        /// <param name="hashProvider"></param>
        /// <param name="input"></param>
        /// <returns></returns>
        public static string GetHash(HashAlgorithm hashProvider, string input)
        {

            // Convert the input string to a byte array and compute the hash. 
            byte[] data = hashProvider.ComputeHash(Encoding.UTF8.GetBytes(input));

            // Create a new Stringbuilder to collect the bytes 
            // and create a string.
            StringBuilder sBuilder = new StringBuilder(data.Length*2);

            // Loop through each byte of the hashed data  
            // and format each one as a hexadecimal string. 
            for (int i = 0; i < data.Length; i++)
                sBuilder.Append(data[i].ToString("X2"));

            // Return the hexadecimal string. 
            return sBuilder.ToString();
        }
    }
}
