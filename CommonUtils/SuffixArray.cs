using System;
using System.Collections.Generic;
using System.Text;
using System.Runtime.InteropServices;

namespace SA
{
    public class SuffixArray
    {
        int[] _input;
        int _inputRealSize;
        int _k;

        int[] _suffixArray;
        //int[] _suffixArrayLookup;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="input">Input data, values 1..k, the input must not contain zeros '0'. The input array must be longer than the actual input</param>
        /// <param name="n"></param>
        /// <param name="k"></param>
        public SuffixArray(int[] input, int n, int k)
        {
            if (input.Length < n + 10)
                throw new Exception("Input buffer size must contains at least 10 items beyond the real values");

            for (int i = 0; i < n; i++)
            {
                if (input[i] == 0)
                    throw new Exception(string.Format("The input must not contain zeros (0), problem in position: {0}", i));
                if ((input[i] < 1) || (input[i] > k))
                    throw new Exception(string.Format("The input must be 1..k, problem in position: {0}, value: {1}", i, input[i]));
            }

            // Reset the buffer items in the input;
            for (int i = n; i < input.Length; i++)
                input[i] = 0;

            _input = input;
            _inputRealSize = n;
            _k = k + 1;
        }

        public void CreateSuffixArray()
        {
            throw new NotImplementedException();
        }


        public int FindSequence(int[] sequence)
        {
            int startRangeMin = -1;
            int startRangeMax = -1;

            return FindSequence(sequence, ref startRangeMin, ref startRangeMax);
        }

        public int FindSequence(int[] sequence, ref int startRangeMin, ref int startRangeMax)
        {
            return FindSequence(sequence, 0, sequence.Length, ref startRangeMin, ref startRangeMax);
        }

        public int FindSequence(int[] sequence, int sequenceStart, int sequenceLength)
        {
            int startRangeMin = -1;
            int startRangeMax = -1;

            return FindSequence(sequence, 0, sequence.Length, ref startRangeMin, ref startRangeMax);
        }

        public int FindSequence(int[] sequence, int sequenceStart, int sequenceLength, ref int startRangeMin, ref int startRangeMax, int compareSequenceStringStartPosition = 0)
        {
            int fs = FindSequenceSA(sequence, sequenceStart, sequenceLength, ref startRangeMin, ref startRangeMax, compareSequenceStringStartPosition);
            if (fs < 0)
                return fs;

            return _suffixArray[fs];
        }

        public int FindSequenceExclude(int[] sequence, int sequenceStart, int sequenceLength, ref int startRangeMin, ref int startRangeMax, int compareSequenceStringStartPosition, int excludePosition)
        {
            int fs = FindSequenceSAExcludePosition(sequence, sequenceStart, sequenceLength, ref startRangeMin, ref startRangeMax, compareSequenceStringStartPosition, excludePosition);
            if (fs < 0)
                return fs;

            return _suffixArray[fs];
        }

        /// <summary>
        /// This function assume that the binary search was performed, and now all left is to return the ALL
        /// the values in range
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="sequenceStart"></param>
        /// <param name="sequenceLength"></param>
        /// <param name="rangeMin"></param>
        /// <param name="rangeMax"></param>
        /// <param name="excludePosition"></param>
        /// <returns></returns>
        public List<int> FindAllSequencesInRange(int[] sequence, int sequenceStart, int sequenceLength, ref int startRangeMin, ref int startRangeMax, bool skipCheckingMidPoint, int compareSequenceStringStartPosition, int excludePosition)
        {
            List<int> positions = new List<int>();

            // Binary search sequence in source
            int rangeMin = 0;
            if (startRangeMin >= 0)
                rangeMin = startRangeMin;

            int rangeMax = _inputRealSize;
            if (startRangeMax >= 0)
                rangeMax = startRangeMax;

            if (rangeMax == rangeMin)
                throw new ArgumentException("No range to search");

            int rangeMid = rangeMin + (rangeMax - rangeMin) / 2;
            int compValue;

            if (skipCheckingMidPoint == false)
            {
                compValue = CompareSequence(sequence, sequenceStart, sequenceLength, _suffixArray[rangeMid], compareSequenceStringStartPosition);
                if (compValue != 0)
                    throw new ArgumentException("Range middle does not contain the sequence");
            }

            if (_suffixArray[rangeMid] != excludePosition)
                positions.Add(_suffixArray[rangeMid]);

            // Search for sequences below the mid point
            for (int i = rangeMid - 1; i >= rangeMin; i--)
            {
                compValue = CompareSequence(sequence, sequenceStart, sequenceLength, _suffixArray[i], compareSequenceStringStartPosition);
                if (compValue != 0)
                {
                    startRangeMin = i + 1;
                    break;
                }

                if (_suffixArray[i] != excludePosition)
                    positions.Add(_suffixArray[i]);
            }

            // Search for sequences above the mid point
            for (int i = rangeMid + 1; i < rangeMax; i++)
            {
                compValue = CompareSequence(sequence, sequenceStart, sequenceLength, _suffixArray[i], compareSequenceStringStartPosition);
                if (compValue != 0)
                {
                    startRangeMax = i;
                    break;
                }

                if (_suffixArray[i] != excludePosition)
                    positions.Add(_suffixArray[i]);
            }

            return positions;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns>If sequence not found -1, else a position within the SuffixArray which contains
        /// the sequence, possibly not the first or last possition, but somewhere in the middle</returns>
        public int FindSequenceSA(int[] sequence, int sequenceStart, int sequenceLength, ref int startRangeMin, ref int startRangeMax, int compareSequenceStringStartPosition)
        {
            // Binary search sequence in source
            int rangeMin = 0;
            if (startRangeMin >= 0)
                rangeMin = startRangeMin;

            int rangeMax = _inputRealSize;
            if (startRangeMax >= 0)
                rangeMax = startRangeMax;

            //// Use the suffix array lookup table
            //if ((startRangeMax == -1) && (startRangeMin == -1))
            //{
            //    int sequenceLookupPos = GetLookupPosition(sequence, sequenceStart, sequenceLength);
            //    if (sequenceLookupPos >= 0)
            //    {
            //        rangeMin = _suffixArrayLookup[sequenceLookupPos];
            //        rangeMax = rangeMin;
            //        while (rangeMax == rangeMin)
            //        {
            //            if (_suffixArrayLookup.Length == sequenceLookupPos + 1)
            //                break;

            //            rangeMax = _suffixArrayLookup[++sequenceLookupPos];
            //        }
            //    }
            //}

            while (rangeMax != rangeMin)
            {
                int rangeMid = rangeMin + (rangeMax - rangeMin) / 2;
                int compValue = CompareSequence(sequence, sequenceStart, sequenceLength, _suffixArray[rangeMid], compareSequenceStringStartPosition);
                if (compValue == 0)
                {
                    startRangeMin = rangeMin;
                    startRangeMax = rangeMax;

                    return rangeMid;
                }

                if (compValue < 0)
                    rangeMax = rangeMid;
                else // compValue > 0
                {
                    // Not found
                    if (rangeMin == rangeMid)
                        return -1;

                    rangeMin = rangeMid;
                }
            }

            return -1;
        }

        public int FindSequenceSAExcludePosition(int[] sequence, int sequenceStart, int sequenceLength, ref int startRangeMin, ref int startRangeMax, int compareSequenceStringStartPosition, int excludePosition)
        {
            // Binary search sequence in source
            int rangeMin = 0;
            if (startRangeMin >= 0)
                rangeMin = startRangeMin;

            int rangeMax = _inputRealSize;
            if (startRangeMax >= 0)
                rangeMax = startRangeMax;

            while (rangeMax != rangeMin)
            {
                int rangeMid = rangeMin + (rangeMax - rangeMin) / 2;
                int compValue = CompareSequence(sequence, sequenceStart, sequenceLength, _suffixArray[rangeMid], compareSequenceStringStartPosition);
                if (compValue == 0)
                {
                    startRangeMin = rangeMin;
                    startRangeMax = rangeMax;

                    // If range-mid is the excluded item
                    if (_suffixArray[rangeMid] == excludePosition)
                    {
                        // Search below the excluded point for a match
                        int tempRangeMin1 = rangeMin;
                        int tempRangeMax1 = rangeMid;
                        int res1 = FindSequenceSA(sequence, sequenceStart, sequenceLength, ref tempRangeMin1, ref tempRangeMax1, compareSequenceStringStartPosition);

                        // Search above the excluded point for a match
                        int tempRangeMin2 = rangeMid + 1;
                        int tempRangeMax2 = rangeMax;
                        int res2 = FindSequenceSA(sequence, sequenceStart, sequenceLength, ref tempRangeMin2, ref tempRangeMax2, compareSequenceStringStartPosition);

                        // No match other than the excluded point
                        if ((res1 == -1) && (res2 == -1))
                            return -1;
                        // A single match - res1
                        else if ((res1 != -1) && (res2 == -1))
                        {
                            startRangeMin = tempRangeMin1;
                            startRangeMax = tempRangeMax1;
                            return res1;
                        }
                        // A single match - res2
                        else if ((res1 == -1) && (res2 != -1))
                        {
                            startRangeMin = tempRangeMin2;
                            startRangeMax = tempRangeMax2;
                            return res2;
                        }
                        else // res1 != -1 && res2 != -1 ==> Nothing to do, continue with the full range, return res1 position
                            return res1;
                    }

                    return rangeMid;
                }

                if (compValue < 0)
                    rangeMax = rangeMid;
                else // compValue > 0
                {
                    // Not found
                    if (rangeMin == rangeMid)
                        return -1;

                    rangeMin = rangeMid;
                }
            }

            return -1;
        }

        [DllImport("msvcrt.dll")]
        static extern int memcmp(IntPtr ptr1, IntPtr ptr2, int count);

        // Check if 'sequence' can be found in 'source' at the given start position.
        private int CompareSequence(int[] sequence, int sequenceStart, int sequenceLength, int startPoint, int compareSequenceStringStartPosition)
        {
            // WARNING: This is x64 specific code
            unsafe
            {
                fixed (int* pSequence = &sequence[sequenceStart + compareSequenceStringStartPosition])
                fixed (int* pIntput = &_input[startPoint + compareSequenceStringStartPosition])
                {
                    IntPtr intPtrSequence = new IntPtr((void*)pSequence);
                    IntPtr intPtrInput = new IntPtr((void*)pIntput);
                    return memcmp(intPtrSequence, intPtrInput, (sequenceLength - compareSequenceStringStartPosition) * 4);
                }
            }

            /*
            for (int i = compareSequenceStringStartPosition; i < sequenceLength; i++)
            {
                // Equal so far, but the Sequence cotinues, and the Input is over => Sequence > Input
                if (_inputRealSize <= startPoint + i)
                    return 1;

                if (sequence[sequenceStart + i] < _input[startPoint + i])
                    return -1;

                if (sequence[sequenceStart + i] > _input[startPoint + i])
                    return 1;
            }

            return 0;
            */
        }
    }
}
