using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace SA
{
    /// <summary>
    /// T - A
    /// G - C
    /// </summary>
    public abstract class DNA : IDisposable  //: IComparable<DNA>, IComparable, IEquatable<DNA>
    {
        public abstract int Length
        {
            get;
        }

        public abstract List<int> FindSubDna(MemoryDNA other);

        public static DNA LoadChromFile(string fileName)
        {
            if (string.IsNullOrEmpty(fileName))
                throw new Exception("No DNA file selected");

            if (File.Exists(fileName) == false)
                throw new Exception("Selected DNA file does not exist");

            List<char> dna = new List<char>();
            using (StreamReader sr = new StreamReader(fileName))
            {
                // Skip the first line
                string line;

                while ((line = sr.ReadLine()) != null)
                {
                    if (line.StartsWith(">"))
                        continue;

                    for (int i = 0; i < line.Length; i++)
                    {
                        switch (line[i])
                        {
                            case 'a':
                            case 'A':
                                dna.Add('A');
                                break;
                            case 't':
                            case 'T':
                                dna.Add('T');
                                break;
                            case 'c':
                            case 'C':
                                dna.Add('C');
                                break;
                            case 'g':
                            case 'G':
                                dna.Add('G');
                                break;
                            default:
                                dna.Add('N');
                                break;
                        }
                    }
                    /*
                    line = line.ToUpper();

                    List<char> replaceChars = new List<char>();
                    for (int i = 0; i < line.Length; i++)
                        if ((line[i] != 'A') && (line[i] != 'T') && (line[i] != 'C') && (line[i] != 'G'))
                            replaceChars.Add(line[i]);

                    for (int i = 0; i < replaceChars.Count; i++)
                        line = line.Replace(replaceChars[i].ToString(), "N");

                    dna.AddRange(line.ToCharArray());
                    */
                }
            }

            return new MemoryDNA(dna);
        }

        public static DNA GenerateRandom(int length)
        {
            char[] bases = new char[] { 'A', 'C', 'G', 'T' };

            List<char> dna = new List<char>();
            Random rand = new Random((int)DateTime.Now.Ticks);

            for (int i = 0; i < length; i++)
                dna.Add(bases[rand.Next(4)]);

            return new MemoryDNA(dna);
        }

        #region Find Inverted Reversed K-mer
        public List<KmerRankData> GetInvertedReverseKmerRanks(int windowSize, int step)
        {
            if (step <= 0)
                step = 1;

            List<KmerRankData> windowsKmerRanks = new List<KmerRankData>();

            for (int startPos = 0; startPos + windowSize <= this.Length; startPos += step)
            {
                List<RepeatData> invertedRepeatList = MemoryDNA.GetInvertedReverseKmerRankDataList(this.GetSubSequence(startPos, windowSize));
                List<int> rankList = new List<int>();
                for (int i = 0; i < invertedRepeatList.Count; i++)
                    rankList.Add(invertedRepeatList[i].ArmLength);

                KmerRankData kmerRankData = new KmerRankData();
                kmerRankData.WindowRank = MathHelper.Avg(rankList);
                kmerRankData.LongestSequence = MathHelper.Max(rankList);

                windowsKmerRanks.Add(kmerRankData);
            }

            return windowsKmerRanks;
        }

        #endregion

        #region Find Inverted K-mer
        public List<KmerRankData> GetInvertedKmerRanks(int windowSize, int step)
        {
            if (step <= 0)
                step = 1;

            List<KmerRankData> windowsKmerRanks = new List<KmerRankData>();

            for (int startPos = 0; startPos + windowSize <= this.Length; startPos += step)
            {
                List<int> rankList = MemoryDNA.GetInvertedKmerRankList(this.GetSubSequence(startPos, windowSize));

                KmerRankData kmerRankData = new KmerRankData();
                kmerRankData.WindowRank = MathHelper.Avg(rankList);
                kmerRankData.LongestSequence = MathHelper.Max(rankList);

                windowsKmerRanks.Add(kmerRankData);
            }

            return windowsKmerRanks;
        }

        #endregion

        #region DNA Int & Char manipulation functions
        public static int GetIntValue(char ch)
        {
            switch (ch)
            {
                case 'T':
                case 't':
                    return 1;
                case 'G':
                case 'g':
                    return 2;
                case 'C':
                case 'c':
                    return 3;
                case 'A':
                case 'a':
                    return 4;

                case 'N':
                case 'n':
                    return 5;
                case 'U':
                case 'u':
                    return 6;
            }

            throw new Exception(string.Format("Invalid DNA building block: {0}", ch));
        }

        public static char GetCharValue(int i)
        {
            switch (i)
            {
                case 1:
                    return 'T';
                case 2:
                    return 'G';
                case 3:
                    return 'C';
                case 4:
                    return 'A';

                // ------------------
                case 5:
                    return 'N';
                case 6:
                    return 'U';
            }

            throw new Exception("Invalid DNA building block");
        }

        protected static int GetComplement(int i)
        {
            switch (i)
            {
                case 1:
                    return 4;
                case 2:
                    return 3;
                case 3:
                    return 2;
                case 4:
                    return 1;

                case 5:
                    return 6;
                case 6:
                    return 5;
            }

            throw new Exception("Unknown value: " + i.ToString());
        }
        #endregion

        public abstract MemoryDNA GetSubSequence(int startPos, int length);

        #region IDisposable Members

        public void Dispose()
        {
        }

        #endregion
    }

    public class MemoryDNA : DNA, ICloneable //, IEquatable<MemoryDNA>
    {
        private static readonly int DNA_VALUE_U = GetIntValue('U');
        private static readonly int DNA_VALUE_N = GetIntValue('N');

        int[] _dna;

        public override int Length
        {
            get { return _dna.Length; }
        }

        public char this[int idx]
        {
            get
            {
                if ((idx < 0) || (idx >= this.Length))
                    throw new ArgumentOutOfRangeException("idx");

                return GetCharValue(_dna[idx]);
            }
        }

        #region C-tors
        public MemoryDNA(string dna)
            : this(dna.ToCharArray()) { }

        public MemoryDNA(char[] dna)
            : this(new List<char>(dna))
        {
        }

        public MemoryDNA(List<char> dna)
        {
            _dna = new int[dna.Count];
            int i = 0;
            try
            {
                for (i = 0; i < dna.Count; i++)
                    _dna[i] = GetIntValue(dna[i]);
            }
            catch (Exception ex)
            {
                throw new Exception(string.Format("Position: {0}, Exception: {1}", i, ex.Message), ex);
            }
        }

        private MemoryDNA(int[] dna)
        {
            _dna = dna;
        }
        #endregion

        /// <summary>
        /// Convert all T to A
        /// Convert all G to C
        /// </summary>
        public void ConvertToBinaryDNA()
        {
            int valA = GetIntValue('A');
            int valT = GetIntValue('T');
            int valG = GetIntValue('G');
            int valC = GetIntValue('C');

            for (int i = 0; i < _dna.Length; i++)
            {
                int pos = _dna[i];
                if (pos == valT)
                    _dna[i] = valA;
                else if (pos == valG)
                    _dna[i] = valC;
            }
        }

        public void WriteToFile(string path)
        {
            using (System.IO.StreamWriter sw = new System.IO.StreamWriter(path))
            {
                int chunkSize = 80;
                for (int startPos = 0; startPos < this.Length; startPos += chunkSize)
                {
                    int len = chunkSize;
                    if (startPos + len >= this.Length)
                        len = this.Length - startPos;
                    sw.WriteLine(this.GetSubSequence(startPos, len).ToString());
                }
            }
        }

        public static List<RepeatData> GetInvertedReverseKmerRankDataList(MemoryDNA sequence)
        {
            return GetInvertedReverseKmerRankDataList(sequence, sequence, true);
        }

        public static List<RepeatData> GetInvertedReverseKmerRankDataList(MemoryDNA sequenceToFindIn, MemoryDNA sequenceForSecondArm, bool selfSequenceSearch = false, int minimumSequenceLength = 50)
        {
            int[] window = null;
            using (MemoryDNA complement = sequenceForSecondArm.GetComplementReverse())
            {
                // Create a suffix array for the complement reversed DNA window
                window = new int[complement.Length + 10];
                Array.Clear(window, 0, window.Length);
                Array.Copy(complement._dna, 0, window, 0, complement.Length);
            }

            return ProcessKmerRankDataList(true, selfSequenceSearch, sequenceToFindIn, window, sequenceForSecondArm.Length, minimumSequenceLength);

            /*
            Dictionary<RepeatData, bool> maxKmerLengthFound = new Dictionary<RepeatData, bool>();

            // 6 = A, C, G, T, N&U
            SuffixArray sa = new SuffixArray(window, sequenceForSecondArm.Length, 6);
            sa.CreateSuffixArray();

            int lastKmerLength = 0;
            for (int i = 0; i < sequenceToFindIn.Length; i++)
            {
                // Optimization: 
                // kmer length can be more than 1 sorter than the last one - reduce 1 every time
                // and start checking from here.
                int kmerLength = Math.Max(0, lastKmerLength - 1);

                int pos = -1;
                int lastFoundPos = -1;
                int startRangeMin = -1;
                int startRangeMax = -1;

                do
                {
                    lastFoundPos = pos; 
                    kmerLength++;

                    if (kmerLength <= sequenceToFindIn.Length - i)
                    {                        
                        pos = sa.FindSequence(sequenceToFindIn._dna, i, kmerLength, ref startRangeMin, ref startRangeMax);
                    }
                    else
                        pos = -1;
                }
                while (pos >= 0);

                lastKmerLength = kmerLength - 1;

                // lastFoundPos == -1 if the current IR is sub-IR of the previous IR.
                if ((lastFoundPos >= 0) && (lastKmerLength >= minimumSequenceLength))
                {
                    int firstArmPosition = i;
                    int secondArmPosition = sequenceForSecondArm.Length - lastFoundPos - lastKmerLength;

                    // Swap - only if searching in self
                    if ((firstArmPosition > secondArmPosition) && (sequenceToFindIn == sequenceForSecondArm))
                    {
                        int temp = firstArmPosition;
                        firstArmPosition = secondArmPosition;
                        secondArmPosition = temp;
                    }

                    RepeatData invertedRepeatData = new RepeatData(firstArmPosition, secondArmPosition, lastKmerLength);
                    if (maxKmerLengthFound.ContainsKey(invertedRepeatData) == false)
                        maxKmerLengthFound.Add(invertedRepeatData, true);
                }

                if (i % 1000000 == 0)
                    System.Console.WriteLine("{0} {1}", DateTime.Now, i);
            }

            return new List<RepeatData>(maxKmerLengthFound.Keys);
            */
        }

        public static List<RepeatData> GetReverseRepeatsKmerRankDataList(MemoryDNA sequence)
        {
            return GetReverseRepeatsKmerRankDataList(sequence, sequence, true);
        }

        public static List<RepeatData> GetReverseRepeatsKmerRankDataList(MemoryDNA sequenceToFindIn, MemoryDNA sequenceForSecondArm, bool selfSequenceSearch = false, int minimumSequenceLength = 50)
        {
            int[] window = null;
            using (MemoryDNA complement = sequenceForSecondArm.GetReverse())
            {
                // Create a suffix array for the complement reversed DNA window
                window = new int[complement.Length + 10];
                Array.Clear(window, 0, window.Length);
                Array.Copy(complement._dna, 0, window, 0, complement.Length);
            }

            return ProcessKmerRankDataList(true, selfSequenceSearch, sequenceToFindIn, window, sequenceForSecondArm.Length, minimumSequenceLength);

            /*
            Dictionary<RepeatData, bool> maxKmerLengthFound = new Dictionary<RepeatData, bool>();

            // 6 = A, C, G, T, N&U
            SuffixArray sa = new SuffixArray(window, sequenceForSecondArm.Length, 6);
            sa.CreateSuffixArray();

            int dnaValueU = GetIntValue('U');
            int dnaValueN = GetIntValue('N');

            int lastKmerLength = 0;
            for (int i = 0; i < sequenceToFindIn.Length; i++)
            {
                // Optimization: 
                // kmer length can be more than 1 sorter than the last one - reduce 1 every time
                // and start checking from here.
                int kmerLength = Math.Max(0, lastKmerLength - 1);

                int pos = -1;
                int lastFoundPos = -1;
                int startRangeMin = -1;
                int startRangeMax = -1;

                do
                {
                    lastFoundPos = pos;
                    kmerLength++;

                    if (kmerLength > sequenceToFindIn.Length - i)
                        pos = -1;
                    else if ((sequenceToFindIn._dna[i + kmerLength - 1] == dnaValueN) || (sequenceToFindIn._dna[i + kmerLength - 1] == dnaValueU))
                        pos = -1;
                    else
                        pos = sa.FindSequence(sequenceToFindIn._dna, i, kmerLength, ref startRangeMin, ref startRangeMax);
                }
                while (pos >= 0);

                lastKmerLength = kmerLength - 1;

                // lastFoundPos == -1 if the current IR is sub-IR of the previous IR.
                if ((lastFoundPos >= 0) && (lastKmerLength >= minimumSequenceLength))
                {
                    int firstArmPosition = i;
                    int secondArmPosition = sequenceForSecondArm.Length - lastFoundPos - lastKmerLength;

                    // Swap - only if searching in self
                    if ((firstArmPosition > secondArmPosition) && (sequenceToFindIn == sequenceForSecondArm))
                    {
                        int temp = firstArmPosition;
                        firstArmPosition = secondArmPosition;
                        secondArmPosition = temp;
                    }

                    RepeatData invertedRepeatData = new RepeatData(firstArmPosition, secondArmPosition, lastKmerLength);
                    if (maxKmerLengthFound.ContainsKey(invertedRepeatData) == false)
                        maxKmerLengthFound.Add(invertedRepeatData, true);
                }

                if (i % 1000000 == 0)
                    System.Console.WriteLine("{0} {1}", DateTime.Now, i);
            }

            return new List<RepeatData>(maxKmerLengthFound.Keys);
            */
        }

        public static List<RepeatData> GetInvertedKmerRankDataList(MemoryDNA sequence)
        {
            return GetInvertedKmerRankDataList(sequence, sequence, true);
        }

        public static List<RepeatData> GetInvertedKmerRankDataList(MemoryDNA sequenceToFindIn, MemoryDNA sequenceForSecondArm, bool selfSequenceSearch = false, int minimumSequenceLength = 50)
        {
            int[] window = null;
            using (MemoryDNA complement = sequenceForSecondArm.GetComplement())
            {
                // Create a suffix array for the complement reversed DNA window
                window = new int[complement.Length + 10];
                Array.Clear(window, 0, window.Length);
                Array.Copy(complement._dna, 0, window, 0, complement.Length);
            }

            return ProcessKmerRankDataList(false, selfSequenceSearch, sequenceToFindIn, window, sequenceForSecondArm.Length, minimumSequenceLength);

            /*
            Dictionary<RepeatData, bool> maxKmerLengthFound = new Dictionary<RepeatData, bool>();

            // 6 = A, C, G, T, N&U
            SuffixArray sa = new SuffixArray(window, sequenceForSecondArm.Length, 6);
            sa.CreateSuffixArray();

            int dnaValueU = GetIntValue('U');
            int dnaValueN = GetIntValue('N');

            int lastKmerLength = 0;
            for (int i = 0; i < sequenceToFindIn.Length; i++)
            {
                // Optimization: 
                // kmer length can be more than 1 sorter than the last one - reduce 1 every time
                // and start checking from here.
                int kmerLength = Math.Max(0, lastKmerLength - 1);

                int pos = -1;
                int lastFoundPos = -1;
                int startRangeMin = -1;
                int startRangeMax = -1;

                do
                {
                    lastFoundPos = pos;
                    kmerLength++;

                    if (kmerLength > sequenceToFindIn.Length - i)
                        pos = -1;
                    else if ((sequenceToFindIn._dna[i + kmerLength - 1] == dnaValueN) || (sequenceToFindIn._dna[i + kmerLength - 1] == dnaValueU))
                        pos = -1;
                    else
                        pos = sa.FindSequenceExclude(sequenceToFindIn._dna, i, kmerLength, ref startRangeMin, ref startRangeMax, i);
                }
                while (pos >= 0);

                lastKmerLength = kmerLength - 1;

                // lastFoundPos == -1 if the current IR is sub-IR of the previous IR.
                if ((lastFoundPos >= 0) && (lastKmerLength >= minimumSequenceLength))
                {
                    int firstArmPosition = i;
                    int secondArmPosition = lastFoundPos;

                    // Swap - only if searching in self
                    if ((firstArmPosition > secondArmPosition) && (sequenceToFindIn == sequenceForSecondArm))
                    {
                        int temp = firstArmPosition;
                        firstArmPosition = secondArmPosition;
                        secondArmPosition = temp;
                    }

                    RepeatData invertedRepeatData = new RepeatData(firstArmPosition, secondArmPosition, lastKmerLength);
                    if (maxKmerLengthFound.ContainsKey(invertedRepeatData) == false)
                        maxKmerLengthFound.Add(invertedRepeatData, true);
                }

                if (i % 1000000 == 0)
                    System.Console.WriteLine("{0} {1}", DateTime.Now, i);
            }

            return new List<RepeatData>(maxKmerLengthFound.Keys);
            */
        }

        public static List<RepeatData> GetRepeatKmerRankDataList(MemoryDNA sequence)
        {
            return GetRepeatKmerRankDataList(sequence, sequence, true);
        }

        public static List<RepeatData> GetRepeatKmerRankDataList(MemoryDNA sequenceToFindIn, MemoryDNA sequenceForSecondArm, bool selfSequenceSearch = false, int minimumSequenceLength = 50)
        {
            int[] window = null;

            // Create a suffix array for the DNA window
            window = new int[sequenceForSecondArm.Length + 10];
            Array.Clear(window, 0, window.Length);
            Array.Copy(sequenceForSecondArm._dna, 0, window, 0, sequenceForSecondArm.Length);

            return ProcessKmerRankDataList(false, selfSequenceSearch, sequenceToFindIn, window, sequenceForSecondArm.Length, minimumSequenceLength);
            /*
            Dictionary<RepeatData, bool> maxKmerLengthFound = new Dictionary<RepeatData, bool>();

            // 6 = A, C, G, T, N&U
            SuffixArray sa = new SuffixArray(window, sequenceForSecondArm.Length, 6);
            sa.CreateSuffixArray();

            int dnaValueU = GetIntValue('U');
            int dnaValueN = GetIntValue('N');

            int lastKmerLength = 0;
            for (int i = 0; i < sequenceToFindIn.Length; i++)
            {
                // Optimization: 
                // kmer length can be more than 1 sorter than the last one - reduce 1 every time
                // and start checking from here.
                int kmerLength = Math.Max(0, lastKmerLength - 1);

                int pos = -1;
                int lastFoundPos = -1;
                int startRangeMin = -1;
                int startRangeMax = -1;

                do
                {
                    lastFoundPos = pos; 
                    kmerLength++;

                    if (kmerLength > sequenceToFindIn.Length - i)
                        pos = -1;
                    else if ((sequenceToFindIn._dna[i + kmerLength - 1] == dnaValueN) || (sequenceToFindIn._dna[i + kmerLength - 1] == dnaValueU))
                        pos = -1;
                    else
                        pos = sa.FindSequenceExclude(sequenceToFindIn._dna, i, kmerLength, ref startRangeMin, ref startRangeMax, i);
                }
                while (pos >= 0);

                lastKmerLength = kmerLength - 1;

                // lastFoundPos == -1 if the current IR is sub-IR of the previous IR.
                if ((lastFoundPos >= 0) && (lastKmerLength >= minimumSequenceLength))
                {
                    int firstArmPosition = i;
                    int secondArmPosition = lastFoundPos;

                    // Swap - only if searching in self
                    if ((firstArmPosition > secondArmPosition) && (sequenceToFindIn == sequenceForSecondArm))
                    {
                        int temp = firstArmPosition;
                        firstArmPosition = secondArmPosition;
                        secondArmPosition = temp;
                    }

                    RepeatData invertedRepeatData = new RepeatData(firstArmPosition, secondArmPosition, lastKmerLength);
                    if (maxKmerLengthFound.ContainsKey(invertedRepeatData) == false)
                        maxKmerLengthFound.Add(invertedRepeatData, true);
                }

                if (i % 1000000 == 0)
                    System.Console.WriteLine("{0} {1}", DateTime.Now, i);
            }

            return new List<RepeatData>(maxKmerLengthFound.Keys);
            */
        }

        /*
        // This function finds only the maximal(s) kmers per each location
        private static List<RepeatData> ProcessKmerRankDataList(bool invertedSearch, bool selfSequenceSeach, MemoryDNA sequenceToFindIn, int[] sequenceForSecondArmWindow, int sequenceForSecondArmLength, int minimumSequenceLength)
        {
            Dictionary<RepeatData, bool> maxKmerLengthFound = new Dictionary<RepeatData, bool>();

            // 6 = A, C, G, T, N&U
            SuffixArray sa = new SuffixArray(sequenceForSecondArmWindow, sequenceForSecondArmLength, 6);
            sa.CreateSuffixArray();

            int dnaValueU = GetIntValue('U');
            int dnaValueN = GetIntValue('N');

            int lastKmerLength = 0;
            for (int i = 0; i < sequenceToFindIn.Length; i++)
            {
                // Optimization: 
                // kmer length can be more than 1 sorter than the last one - reduce 1 every time
                // and start checking from here.
                int kmerLength = Math.Max(minimumSequenceLength - 1, lastKmerLength - 1);

                int pos = -1;
                int lastFoundPos = -1;
                int startRangeMin = -1;
                int startRangeMax = -1;

                do
                {
                    lastFoundPos = pos;
                    kmerLength++;

                    if (kmerLength > sequenceToFindIn.Length - i)
                        pos = -1;
                    else if ((sequenceToFindIn._dna[i + kmerLength - 1] == dnaValueN) || (sequenceToFindIn._dna[i + kmerLength - 1] == dnaValueU))
                        pos = -1;
                    else
                    {
                        if ((selfSequenceSeach == true) && (invertedSearch == false))
                        {
                            // Prevent from finding the same position as a repeat incase of searching self 
                            pos = sa.FindSequenceExclude(sequenceToFindIn._dna, i, kmerLength, ref startRangeMin, ref startRangeMax, i);
                        }
                        else
                            pos = sa.FindSequence(sequenceToFindIn._dna, i, kmerLength, ref startRangeMin, ref startRangeMax);
                    }
                }
                while (pos >= 0);

                lastKmerLength = kmerLength - 1;

                // lastFoundPos == -1 if the current IR is sub-IR of the previous IR.
                if ((lastFoundPos >= 0) && (lastKmerLength >= minimumSequenceLength))
                {
                    int excludePosition = -1;
                    if ((selfSequenceSeach == true) && (invertedSearch == false))
                        excludePosition = i;

                    List<int> positions = sa.FindAllSequencesInRange(sequenceToFindIn._dna, i, lastKmerLength, ref startRangeMin, ref startRangeMax, excludePosition);

                    int firstArmPosition = i;
                    int secondArmPosition;

                    for (int foundPositions = 0; foundPositions < positions.Count; foundPositions++)
                    {
                        if (invertedSearch == true)
                            secondArmPosition = sequenceForSecondArmLength - positions[foundPositions] - lastKmerLength;
                        else
                            secondArmPosition = positions[foundPositions];

                        // Swap - only if searching in self
                        if ((firstArmPosition > secondArmPosition) && (selfSequenceSeach == true))
                        {
                            int temp = firstArmPosition;
                            firstArmPosition = secondArmPosition;
                            secondArmPosition = temp;
                        }

                        RepeatData invertedRepeatData = RepeatData.obtain(firstArmPosition, secondArmPosition, lastKmerLength);
                        if (maxKmerLengthFound.ContainsKey(invertedRepeatData) == false)
                            maxKmerLengthFound.Add(invertedRepeatData, true);
                    }
                }

                if (i % 1000000 == 0)
                    System.Console.WriteLine("{0} {1}", DateTime.Now, i);
            }

            List<RepeatData> list = new List<RepeatData>(maxKmerLengthFound.Keys);
            list.Sort();

            return list;
        }
        */

        private static List<RepeatData> ProcessKmerRankDataList(bool invertedSearch, bool selfSequenceSeach, MemoryDNA sequenceToFindIn, int[] sequenceForSecondArmWindow, int sequenceForSecondArmLength, int minimumSequenceLength)
        {
            Dictionary<RepeatData, bool> maxKmerLengthFound = new Dictionary<RepeatData, bool>();
            List<RepeatData> cachedRecentFoundKmers = new List<RepeatData>();

            System.Console.WriteLine("{0} Starting to build suffix array", DateTime.Now);

            // 6 = A, C, G, T, N&U
            SuffixArray sa = new SuffixArray(sequenceForSecondArmWindow, sequenceForSecondArmLength, 6);
            sa.CreateSuffixArray();

            for (int i = 0; i < sequenceToFindIn.Length; i++)
            {
                // Candidate mapping <OtherArmPosition, KmerLength>
                Dictionary<int, int> candidatesMaxKmerLengthFound = ProcessKmerFindCandidates(cachedRecentFoundKmers, invertedSearch, selfSequenceSeach, sequenceToFindIn, sequenceForSecondArmLength, minimumSequenceLength, sa, i);

                if (candidatesMaxKmerLengthFound != null)
                {
                    // Trim the cache used from trimming duplications - if the arm does not cover newly found arms - no point in being in the cache
                    while (cachedRecentFoundKmers.Count > 0)
                    {
                        RepeatData repeatData = cachedRecentFoundKmers[0];
                        if (repeatData.FirstArmPosition + repeatData.ArmLength >= i + minimumSequenceLength)
                            break;

                        cachedRecentFoundKmers.RemoveAt(0);
                        repeatData.recycle();
                    }

                    // Remove overlapping kmers this loop (e.g. the first sequence completly overlaps the second)
                    List<int> otherArmPositionsSorter = new List<int>(candidatesMaxKmerLengthFound.Keys);
                    otherArmPositionsSorter.Sort();

                    LinkedList<int> otherArmPositions = new LinkedList<int>(otherArmPositionsSorter);
                    LinkedListNode<int> otherArmPositionsNext = otherArmPositions.First;
                    while (otherArmPositionsNext != null)
                    {
                        int otherArmPosition1 = otherArmPositionsNext.Value;
                        int otherArmLength1 = candidatesMaxKmerLengthFound[otherArmPosition1];
                        if (otherArmPositionsNext.Next != null)
                        {
                            int otherArmPosition2 = otherArmPositionsNext.Next.Value;
                            int otherArmLength2 = candidatesMaxKmerLengthFound[otherArmPosition2];
                            if (otherArmPosition1 + otherArmLength1 >= otherArmPosition2 + otherArmLength2)
                            {
                                otherArmPositions.Remove(otherArmPositionsNext.Next);
                                continue;
                            }
                        }

                        bool removedItem = false;
                        for (int cachedRecentIndex = 0; cachedRecentIndex < cachedRecentFoundKmers.Count; cachedRecentIndex++)
                        {
                            RepeatData cachedRecentFoundKmer = cachedRecentFoundKmers[cachedRecentIndex];
                            if (
                                (cachedRecentFoundKmer.FirstArmPosition <= i) &&
                                (cachedRecentFoundKmer.FirstArmPosition + cachedRecentFoundKmer.ArmLength >= i + otherArmLength1) &&
                                (cachedRecentFoundKmer.SecondArmPosition <= otherArmPosition1) &&
                                (cachedRecentFoundKmer.SecondArmPosition + cachedRecentFoundKmer.ArmLength >= otherArmPosition1 + otherArmLength1)
                                )
                            {
                                LinkedListNode<int> itemToRemove = otherArmPositionsNext;
                                otherArmPositionsNext = otherArmPositionsNext.Next;
                                otherArmPositions.Remove(itemToRemove);
                                removedItem = true;
                                break;
                            }
                        }
                        if (removedItem == true)
                            continue;

                        otherArmPositionsNext = otherArmPositionsNext.Next;
                    }

                    for (otherArmPositionsNext = otherArmPositions.First; otherArmPositionsNext != null; otherArmPositionsNext = otherArmPositionsNext.Next)
                    {
                        int secondArmPosition = otherArmPositionsNext.Value;
                        int length = candidatesMaxKmerLengthFound[secondArmPosition];

                        int firstArmPosition = i;

                        // Swap - only if searching in self
                        if ((firstArmPosition > secondArmPosition) && (selfSequenceSeach == true))
                        {
                            int temp = firstArmPosition;
                            firstArmPosition = secondArmPosition;
                            secondArmPosition = temp;
                        }

                        cachedRecentFoundKmers.Add(RepeatData.obtain(i, otherArmPositionsNext.Value, length));

                        RepeatData invertedRepeatData = RepeatData.obtain(firstArmPosition, secondArmPosition, length);
                        if (maxKmerLengthFound.ContainsKey(invertedRepeatData) == false)
                            maxKmerLengthFound.Add(invertedRepeatData, true);
                    }
                }

                if (i % 1000000 == 0)
                    System.Console.WriteLine("{0} {1}", DateTime.Now, i);
            }

            List<RepeatData> list = new List<RepeatData>(maxKmerLengthFound.Keys);
            list.Sort();

            return list;
        }

        private static Dictionary<int, int> ProcessKmerFindCandidates(List<RepeatData> cachedRecentFoundKmers, bool invertedSearch, bool selfSequenceSeach, MemoryDNA sequenceToFindIn, int sequenceForSecondArmLength, int minimumSequenceLength, SuffixArray sa, int i)
        {
            int excludePosition = -1;
            if ((selfSequenceSeach == true) && (invertedSearch == false))
                excludePosition = i;

            int kmerLength = minimumSequenceLength - 1;

            // Candidate mapping <OtherArmPosition, KmerLength>
            Dictionary<int, int> candidatesMaxKmerLengthFound = null;

            // Validate only valid bases
            if (kmerLength > sequenceToFindIn.Length - i)
                return null;
            for (int j = 0; j < kmerLength; j++)
                if ((sequenceToFindIn._dna[i + j] == DNA_VALUE_N) || (sequenceToFindIn._dna[i + j] == DNA_VALUE_U))
                    return null;

            int pos = -1;
            int lastFoundPos = -1;
            int startRangeMin = -1;
            int startRangeMax = -1;

            // Starting position for the 'CompareSequence' function, assuming all the position below
            // the starting position are equal
            int compareSequenceStartPosition = 0;

            do
            {
                lastFoundPos = pos;
                kmerLength++;

                if (kmerLength > sequenceToFindIn.Length - i)
                    pos = -1;
                else if ((sequenceToFindIn._dna[i + kmerLength - 1] == DNA_VALUE_N) || (sequenceToFindIn._dna[i + kmerLength - 1] == DNA_VALUE_U))
                    pos = -1;
                else
                {
                    // Prevent from finding the same position as a repeat incase of searching self 
                    pos = sa.FindSequenceExclude(sequenceToFindIn._dna, i, kmerLength, ref startRangeMin, ref startRangeMax, compareSequenceStartPosition, excludePosition);

                    if ((pos >= 0) && (kmerLength >= minimumSequenceLength))
                    {
                        if (candidatesMaxKmerLengthFound == null)
                            candidatesMaxKmerLengthFound = new Dictionary<int, int>();

                        List<int> positions = sa.FindAllSequencesInRange(sequenceToFindIn._dna, i, kmerLength, ref startRangeMin, ref startRangeMax, true, compareSequenceStartPosition, excludePosition);
                        // After 'FindAllSequencesInRange' it is known that all the items left in 
                        // the SuffixArray in the search range, are equal in their first 'kmerLength'
                        // items
                        compareSequenceStartPosition = kmerLength;

                        int firstArmPosition = i;

                        for (int foundPositions = 0; foundPositions < positions.Count; foundPositions++)
                        {
                            int otherArmPosition;
                            if (invertedSearch == true)
                                otherArmPosition = sequenceForSecondArmLength - positions[foundPositions] - kmerLength;
                            else
                                otherArmPosition = positions[foundPositions];

                            // Super early pruning
                            // Perform this only on the first kmer batch found
                            if (kmerLength == minimumSequenceLength)
                            {
                                bool removedItem = false;
                                for (int cachedRecentIndex = 0; cachedRecentIndex < cachedRecentFoundKmers.Count; cachedRecentIndex++)
                                {
                                    RepeatData cachedRecentFoundKmer = cachedRecentFoundKmers[cachedRecentIndex];
                                    if (
                                        (cachedRecentFoundKmer.FirstArmPosition <= i) &&
                                        (cachedRecentFoundKmer.FirstArmPosition + cachedRecentFoundKmer.ArmLength >= i + kmerLength) &&
                                        (cachedRecentFoundKmer.SecondArmPosition <= otherArmPosition) &&
                                        (cachedRecentFoundKmer.SecondArmPosition + cachedRecentFoundKmer.ArmLength >= otherArmPosition + kmerLength) &&
                                        // This is unique to this triming - the new pair found is a direct alignment with an existing pair
                                        (cachedRecentFoundKmer.FirstArmPosition - i == cachedRecentFoundKmer.SecondArmPosition - otherArmPosition)
                                        )
                                    {
                                        removedItem = true;
                                        break;
                                    }
                                }
                                // Skip entering the item to the candidates
                                if (removedItem == true)
                                    continue;
                            }
                            ///////

                            candidatesMaxKmerLengthFound[otherArmPosition] = kmerLength;
                        }

                        // If the super early pruning removed all items - there's nothing to do here, exist
                        if (candidatesMaxKmerLengthFound.Count == 0)
                            return null;
                    }
                }
            }
            while (pos >= 0);
            return candidatesMaxKmerLengthFound;
        }


        public SuffixArray GetSuffixArray()
        {
            int[] window = new int[this.Length + 10];
            Array.Clear(window, 0, window.Length);
            Array.Copy(this._dna, 0, window, 0, this.Length);

            // 6 = A, C, G, T, N&U
            SuffixArray sa = new SuffixArray(window, this.Length, 6);
            sa.CreateSuffixArray();

            return sa;
        }

        public static List<int> GetInvertedKmerRankList(MemoryDNA sequence)
        {
            List<int> maxKmerLengthFound = new List<int>();

            // Create a suffix array for the complement reversed DNA window
            MemoryDNA complement = sequence.GetComplement();

            // Create a suffix array for the complement reversed DNA window
            int[] window = new int[sequence.Length + 10];
            Array.Clear(window, 0, window.Length);
            Array.Copy(complement._dna, 0, window, 0, sequence.Length);

            // 6 = A, C, G, T, N&U
            SuffixArray sa = new SuffixArray(window, sequence.Length, 6);
            sa.CreateSuffixArray();

            for (int i = 0; i < sequence.Length; i++)
            {
                int kmerLength = 0;

                int pos = 0;
                int startRangeMin = -1;
                int startRangeMax = -1;

                while (pos >= 0)
                {
                    kmerLength++;
                    if (kmerLength <= sequence.Length - i)
                        pos = sa.FindSequence(sequence._dna, i, kmerLength, ref startRangeMin, ref startRangeMax);
                    else
                        pos = -1;
                }

                maxKmerLengthFound.Add(kmerLength - 1);
            }

            return maxKmerLengthFound;
        }

        public static int FindLongestCommonSubsequence(SuffixArray sa, MemoryDNA sequence)
        {
            int maxLengthFound = 0;
            for (int i = 0; i < sequence.Length; i++)
            {
                int kmerLength = maxLengthFound;

                int pos = 0;
                int startRangeMin = -1;
                int startRangeMax = -1;

                while (pos >= 0)
                {
                    kmerLength++;
                    if (kmerLength <= sequence.Length - i)
                        pos = sa.FindSequence(sequence._dna, i, kmerLength, ref startRangeMin, ref startRangeMax);
                    else
                        pos = -1;
                }

                maxLengthFound = Math.Max(maxLengthFound, kmerLength - 1);
            }

            return maxLengthFound;
        }

        /// <summary>
        /// Get the Mirror reversed DNS.
        /// E.g. for AATTCCGG => CCGGAATT
        /// </summary>
        /// <returns></returns>
        public MemoryDNA GetComplementReverse()
        {
            int[] complementDna = new int[this.Length];
            for (int i = 0; i < this.Length; i++)
                complementDna[this.Length - 1 - i] = GetComplement(_dna[i]);

            return new MemoryDNA(complementDna);
        }

        public void PerformInplacePartialInvertedRepeat(int startPos, int length)
        {
            if (length < 0)
                throw new ArgumentException("length < 0");
            if ((startPos < 0) || (startPos >= _dna.Length))
                throw new ArgumentException("startPos out of range");
            if (startPos + length > _dna.Length)
                throw new ArgumentException("startPos+Length out of range");

            int endPos = startPos + length - 1;
            while (startPos <= endPos)
            {
                int startPosValue = _dna[startPos];
                _dna[startPos] = GetComplement(_dna[endPos]);
                if (startPos != endPos)
                    _dna[endPos] = GetComplement(startPosValue);

                startPos++;
                endPos--;
            }
        }

        /// <summary>
        /// Get the Mirror reversed DNS.
        /// E.g. for AATTCCGG => GGCCTTAA
        /// </summary>
        /// <returns></returns>
        public MemoryDNA GetReverse()
        {
            int[] reversedDna = new int[this.Length];
            for (int i = 0; i < this.Length; i++)
                reversedDna[this.Length - 1 - i] = _dna[i];

            return new MemoryDNA(reversedDna);
        }

        /// <summary>
        /// Get the Mirror reversed DNS.
        /// E.g. for AATTCCGG => TTAAGGCC
        /// </summary>
        /// <returns></returns>
        public MemoryDNA GetComplement()
        {
            int[] complementDna = new int[this.Length];
            for (int i = 0; i < this.Length; i++)
                complementDna[i] = GetComplement(_dna[i]);

            return new MemoryDNA(complementDna);
        }

        public override List<int> FindSubDna(MemoryDNA otherMemory)
        {
            List<int> foundPositions = new List<int>();

            for (int i = 0; i < this._dna.Length; i++)
            {
                if (PositionEquals(i, otherMemory) == true)
                    foundPositions.Add(i);
            }

            return foundPositions;
        }

        public override MemoryDNA GetSubSequence(int startPos, int length)
        {
            int[] newDna = new int[length];

            Array.Copy(this._dna, startPos, newDna, 0, length);

            return new MemoryDNA(newDna);
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < this.Length; i++)
                sb.Append(this[i]);

            return sb.ToString();
        }

        #region DNA Comparison Members

        public int CompareTo(object obj)
        {
            DNA other = obj as DNA;
            if (other == null)
                return 1;

            return CompareTo(other);
        }

        public int CompareTo(MemoryDNA other)
        {
            if (this.Length < other.Length)
                return -1;

            if (this.Length > other.Length)
                return 1;

            for (int i = 0; i < this.Length; i++)
            {
                if (this._dna[i] < other._dna[i])
                    return -1;

                if (this._dna[i] > other._dna[i])
                    return 1;
            }

            return 0;
        }

        public bool PositionEquals(int startPos, MemoryDNA other)
        {
            if (this.Length < startPos + other.Length)
                return false;

            for (int i = 0; i < other.Length; i++)
                if (this._dna[startPos + i] != other._dna[i])
                    return false;

            return true;
        }

        public bool Equals(DNA other)
        {
            return PositionEquals(0, other as MemoryDNA);
        }

        #endregion

        public object Clone()
        {
            return new MemoryDNA((int[])_dna.Clone());
        }
    }

    public class FileDNA : DNA
    {
        const int FileBlockSize = 16 * 1024;
        string _fileName;
        int _length;

        // <DNA Index, File Index>
        List<DnaFileMarker> _dnaFileMarkers;

        public override int Length
        {
            get { return _length; }
        }

        public FileDNA(string fileName)
        {
            _fileName = fileName;
            if (File.Exists(fileName) == false)
                throw new FileNotFoundException(fileName);

            // Build DNA file markers
            _dnaFileMarkers = new List<DnaFileMarker>();

            using (FileStream fs = new FileStream(_fileName, FileMode.Open, FileAccess.Read))
            {
                byte[] buffer = new byte[FileBlockSize];

                // Read and skip first line
                while (fs.ReadByte() != 10) ;

                _dnaFileMarkers.Add(new DnaFileMarker(fs.Position, 0));

                int readSize;
                int dnaPosition = 0;
                while ((readSize = fs.Read(buffer, 0, buffer.Length)) > 0)
                {
                    for (int i = 0; i < readSize; i++)
                    {
                        if (buffer[i] != 10)
                            dnaPosition++;
                    }

                    _dnaFileMarkers.Add(new DnaFileMarker(fs.Position, dnaPosition));
                }

                _dnaFileMarkers.Sort();
                _length = dnaPosition;
            }

            _fileReadBuffer = new List<char>();
            _fileReadBufferStartPos = 0;
        }

        int _fileReadBufferStartPos;
        List<char> _fileReadBuffer;

        public override MemoryDNA GetSubSequence(int startPos, int length)
        {
            lock (this)
            {
                if ((_fileReadBufferStartPos > startPos) ||
                    (_fileReadBufferStartPos + _fileReadBuffer.Count < startPos + length))
                {
                    int beforeIndex = _dnaFileMarkers.BinarySearch(new DnaFileMarker(-1, startPos));
                    if (beforeIndex < 0)
                        beforeIndex = ~beforeIndex - 1;

                    using (FileStream fs = new FileStream(_fileName, FileMode.Open, FileAccess.Read))
                    {
                        fs.Seek(_dnaFileMarkers[beforeIndex].FilePosition, SeekOrigin.Begin);
                        _fileReadBufferStartPos = _dnaFileMarkers[beforeIndex].DnaPosition;
                        _fileReadBuffer = LoadFromStream(fs, Math.Max(FileBlockSize * 10, length * 10));
                    }

                }

                try
                {
                    return new MemoryDNA(_fileReadBuffer.GetRange(startPos - _fileReadBufferStartPos, length));
                }
                catch
                {
                    System.Console.WriteLine("{0} {1}, {2} {3}", _fileReadBufferStartPos, _fileReadBuffer.Count,
                        startPos, length);

                    throw;
                }
            }
        }

        List<char> LoadFromStream(FileStream fs, int desiredBufferSize)
        {
            List<char> dna = new List<char>();

            using (StreamReader sr = new StreamReader(fs))
            {
                string line;

                while ((line = sr.ReadLine()) != null)
                {
                    if (line.StartsWith(">"))
                        continue;

                    line = line.ToUpper();

                    List<char> replaceChars = new List<char>();
                    for (int i = 0; i < line.Length; i++)
                        if ((line[i] != 'A') && (line[i] != 'T') && (line[i] != 'C') && (line[i] != 'G'))
                            replaceChars.Add(line[i]);

                    for (int i = 0; i < replaceChars.Count; i++)
                        line = line.Replace(replaceChars[i].ToString(), "N");

                    dna.AddRange(line.ToCharArray());

                    if (dna.Count >= desiredBufferSize)
                        return dna;
                }
            }

            return dna;
        }

        class DnaFileMarker : IComparable<DnaFileMarker>
        {
            long _filePosition;
            int _dnaPosition;

            public long FilePosition
            {
                get { return _filePosition; }
            }

            public int DnaPosition
            {
                get { return _dnaPosition; }
            }

            public DnaFileMarker(long filePosition, int dnaPosition)
            {
                _filePosition = filePosition;
                _dnaPosition = dnaPosition;
            }

            #region IComparable<DnaFileMarker> Members

            public int CompareTo(DnaFileMarker other)
            {
                return this.DnaPosition.CompareTo(other.DnaPosition);
            }

            #endregion
        }

        public override List<int> FindSubDna(MemoryDNA otherMemory)
        {
            List<int> foundPositions = new List<int>();

            for (int i = 0; i < this.Length - otherMemory.Length; i++)
            {
                MemoryDNA thisDnaSubsequence = this.GetSubSequence(i, otherMemory.Length);
                if (thisDnaSubsequence.Equals(otherMemory) == true)
                    foundPositions.Add(i);
            }

            return foundPositions;
        }
    }

    public class KmerRankData
    {
        public double WindowRank;
        public int LongestSequence;
    }

    public class RepeatData : IComparable<RepeatData>, IEquatable<RepeatData>
    {
        int _firstArmPosition;
        int _secondArmPosition;
        int _armLength;

        // Object pool structure
        private RepeatData mNext;
        private static Object gPoolSync = new Object();
        private static RepeatData gPool;
        private static int gPoolSize = 0;
        private const int MAX_POOL_SIZE = 1000;
        //

        public int FirstArmPosition
        {
            get { return _firstArmPosition; }
        }

        public int SecondArmPosition
        {
            get { return _secondArmPosition; }
        }

        public int ArmLength
        {
            get { return _armLength; }
        }

        private RepeatData(int firstArmPosition, int secondArmPosition, int armLength)
        {
            _firstArmPosition = firstArmPosition;
            _secondArmPosition = secondArmPosition;

            _armLength = armLength;
        }

        #region Object Pool
        public static RepeatData obtain(int firstArmPosition, int secondArmPosition, int armLength)
        {
            lock (gPoolSync)
            {
                if (gPool != null)
                {
                    RepeatData m = gPool;
                    m._firstArmPosition = firstArmPosition;
                    m._secondArmPosition = secondArmPosition;
                    m._armLength = armLength;

                    gPool = m.mNext;
                    m.mNext = null;
                    gPoolSize--;
                    return m;
                }
            }

            return new RepeatData(firstArmPosition, secondArmPosition, armLength);
        }

        public void recycle()
        {
            lock (gPoolSync)
            {
                if (gPoolSize < MAX_POOL_SIZE)
                {
                    mNext = gPool;
                    gPool = this;
                    gPoolSize++;
                }
            }
        }
        #endregion

        #region IComparable<InvertedRepeatData> Members

        public int CompareTo(RepeatData other)
        {
            if (this.ArmLength < other.ArmLength)
                return 1;

            if (this.ArmLength > other.ArmLength)
                return -1;

            if (this.FirstArmPosition > other.FirstArmPosition)
                return 1;

            if (this.FirstArmPosition < other.FirstArmPosition)
                return -1;

            if (this.SecondArmPosition > other.SecondArmPosition)
                return 1;

            if (this.SecondArmPosition < other.SecondArmPosition)
                return -1;

            return 0;
        }

        #endregion

        public override int GetHashCode()
        {
            return this.FirstArmPosition;
        }

        public override bool Equals(object obj)
        {
            RepeatData other = obj as RepeatData;
            if (other == null)
                return false;

            return this.Equals(other);
        }

        public override string ToString()
        {
            return string.Format("Len: {0}, Arm1: {1}, Arm2: {2}", this.ArmLength, this.FirstArmPosition, this.SecondArmPosition);
        }

        #region IEquatable<InvertedRepeatData> Members

        public bool Equals(RepeatData other)
        {
            return (this.FirstArmPosition == other.FirstArmPosition) && (this.SecondArmPosition == other.SecondArmPosition);
        }

        #endregion
    }

}
