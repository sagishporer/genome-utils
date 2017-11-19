using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using SA;
using CommonUtils;

namespace CountDNABases
{
    /// <summary>
    /// Given a CoverageRegions.txt in the following format
    /// Chromosome Pos1 Pos2
    /// 
    /// Count, for each chromosome, the number of bases that appear in those regions. It also count sequences (how many AA, CG, etc appear)
    /// </summary>    
    class Program
    {
        const string INPUT_COVERAGE_FILE = @"E:\Genome\HG38Cover.txt";
        const string INPUT_FILE_PATTERN = @"{0}.fa";

        //const string INPUT_FILE_FOLDER = @"C:\Temp\Genome\CCC\frog.xenTro3.fa.masked\";
        const string INPUT_FILE_SEARCH_PATTERN = "*.fa.*";

        const bool JOIN_STAT_FOR_MULTIPLE_FILES = false;
        const bool EXPORT_FULL_DATA = true;
        const bool PROCESS_COVERAGE_FILE = false;
        const bool APPEND_SPECIMEN_NAME = false;

        const int HISTOGRAM_BUCKETS = 250;
        const bool PRINT_HISTOGRAM = false;

        const int MAX_SEQUENCE_LENGTH = 12;

        static readonly char[] BASES = { 'A', 'T', 'G', 'C' };
        static Random gRandom = new Random((int)DateTime.Now.Ticks);

        const bool SELECT_FOLDER_IN_UI = false;

        [STAThread]
        static void Main(string[] args)
        {
            if ((args == null) || (args.Length != 1))
            {
                Console.WriteLine("Usage:");
                Console.WriteLine("CountDNABase [File/Directory]");
                return;
            }

            if (Directory.Exists(args[0]) == false)
            {
                Console.WriteLine("Folder not found: {0}", args[0]);
                return;
            }

            ProcessFolder(args[0]);
        }

        static void BuildLowComplexityChromosomes(string originalFolder, string maskedFolder, string diffFolder)
        {
            if (Directory.Exists(diffFolder) == false)
                Directory.CreateDirectory(diffFolder);

            string[] sourceFiles = Directory.GetFiles(originalFolder, "*.fa");
            foreach (string sourceFile in sourceFiles)
            {
                Console.WriteLine("{0} - Processing file: {1}", DateTime.Now, sourceFile);

                FileInfo fileInfo = new FileInfo(sourceFile);
                string maskedFile = Path.Combine(maskedFolder, fileInfo.Name);
                if (File.Exists(maskedFile) == false)
                {
                    Console.WriteLine("{0} - No masked file found", DateTime.Now);
                    continue;
                }

                MemoryDNA source = (MemoryDNA)DNA.LoadChromFile(sourceFile);
                MemoryDNA masked = (MemoryDNA)DNA.LoadChromFile(maskedFile);

                if (source.Length != masked.Length)
                {
                    Console.WriteLine("{0} - Source & masked chromosomes are not in the same length");
                    continue;
                }

                string outputFile = Path.Combine(diffFolder, fileInfo.Name);
                using (StreamWriter sw = new StreamWriter(outputFile))
                {
                    StringBuilder sb = new StringBuilder();
                    bool lowComplexityRegion = false;

                    for (int i = 0; i < source.Length; i++)
                    {
                        if (masked[i] == 'N')
                        {
                            char ch = source[i];
                            sb.Append(ch);
                            if (ch == 'N')
                                lowComplexityRegion = false;
                            else
                                lowComplexityRegion = true;
                        }
                        else if (lowComplexityRegion == true)
                        {
                            sb.Append('N');
                            lowComplexityRegion = false;
                        }

                        if (sb.Length >= 80)
                        {
                            sw.WriteLine(sb.ToString());
                            sb.Clear();
                        }
                    }

                    if (sb.Length > 0)
                        sw.WriteLine(sb.ToString());
                }
            }
        }

        static void CountKmersAndStatsinChromosome(string folderName, string chromosomeFile, int maxLength, long windowSize, double minimumValidBasesPercent)
        {
            DNAMemoryStreamReader dnaStreamReader = new DNAMemoryStreamReader((MemoryDNA)MemoryDNA.LoadChromFile(chromosomeFile));
            long dnaLength = dnaStreamReader.Length;
            StringBuilder output = new StringBuilder();

            for (int kmerLength = 1; kmerLength <= maxLength; kmerLength++)
            {
                System.Console.WriteLine("{0} - Processing file: {1}, Len: {2}", DateTime.Now, chromosomeFile, kmerLength);

                ulong patternMask = 0;
                for (int j = 0; j < kmerLength; j++)
                    patternMask = (patternMask << 2) + 0x3;

                // Prepare statistics object for specific length
                VLargeArray<StatBuilder> stats = VLargeArray<StatBuilder>.AllocateArray((ulong)Math.Pow(4, kmerLength));
                VLargeArray<StatBuilder> statsN1Nir = VLargeArray<StatBuilder>.AllocateArray((ulong)Math.Pow(4, kmerLength));
                VLargeArray<StatBuilder> statsAbsN1Nir = VLargeArray<StatBuilder>.AllocateArray((ulong)Math.Pow(4, kmerLength));

                List<VLargeArray<StatBuilder>> statList = new List<VLargeArray<StatBuilder>>()
                {
                    stats,
                    statsN1Nir,
                    statsAbsN1Nir
                };

                for (ulong i = 0; i < stats.Length; i++)
                {
                    stats[i] = new StatBuilder("Count: " + GetPatternStringFromValue(i, kmerLength), PRINT_HISTOGRAM, HISTOGRAM_BUCKETS);
                    statsN1Nir[i] = new StatBuilder("(N1-Nir)/(N1+Nir): " + GetPatternStringFromValue(i, kmerLength), PRINT_HISTOGRAM, HISTOGRAM_BUCKETS);
                    statsAbsN1Nir[i] = new StatBuilder("|N1-Nir|/(N1+Nir): " + GetPatternStringFromValue(i, kmerLength), PRINT_HISTOGRAM, HISTOGRAM_BUCKETS);
                }

                // Build statistics for specific kmer length
                for (long dnaPos = 0; dnaPos < dnaLength; dnaPos += windowSize)
                {
                    List<Pair> pairs = new List<Pair>() { new Pair(dnaPos, dnaPos + windowSize) };
                    dnaStreamReader.Seek(0L);
                    VLargeArray<int> kmerCount = BuildKmerCount(dnaStreamReader, pairs, kmerLength);

                    // Validate minimum valid base-pairs (ACGT)
                    long sum = 0;
                    for (ulong i = 0; i < kmerCount.Length; i++)
                        sum += kmerCount[i];
                    if ((double)sum / (double)windowSize < minimumValidBasesPercent)
                        continue;

                    for (ulong i = 0; i < kmerCount.Length; i++)
                    {
                        ulong irPattern = ConvertPattern_BuildIR(kmerLength, patternMask, i);

                        int patternCountI = kmerCount[i];
                        int irCount = kmerCount[irPattern];

                        stats[i].AddValue(patternCountI);
                        if ((irCount > 0) || (patternCountI > 0))
                        {
                            statsN1Nir[i].AddValue((double)(patternCountI - irCount) / (double)(irCount + patternCountI));
                            statsAbsN1Nir[i].AddValue(Math.Abs((double)(patternCountI - irCount) / (double)(irCount + patternCountI)));
                        }
                    }
                }

                string simpleFolderName = new DirectoryInfo(folderName).Name;
                // Print statistics
                foreach (VLargeArray<StatBuilder> statistics in statList)
                    for (ulong i = 0; i < statistics.Length; i++)
                    {
                        if (APPEND_SPECIMEN_NAME == true)
                        {
                            output.Append(simpleFolderName);
                            output.Append('\t');
                        }

                        output.Append(statistics[i].Title);
                        output.Append('\t');
                        output.Append(kmerLength);
                        output.Append('\t');
                        output.Append(statistics[i].GetCount());
                        output.Append('\t');
                        output.Append(statistics[i].GetMinValue());
                        output.Append('\t');
                        output.Append(statistics[i].GetMaxValue());
                        output.Append('\t');
                        output.Append(statistics[i].GetAverage());
                        output.Append('\t');
                        output.Append(statistics[i].GetVariance());

                        output.AppendLine();
                    }
            }

            using (StreamWriter sw = new StreamWriter(chromosomeFile + ".kmerCount.Windows." + windowSize))
            {
                sw.Write(output);
            }
        }

        static void CountKmersInFolder(string folder, int maxLength)
        {
            // Build dataset
            string[] files = Directory.GetFiles(folder, "*.fa");
            SortedDictionary<string, List<VLargeArray<int>>> patternsCount = new SortedDictionary<string, List<VLargeArray<int>>>();
            foreach (string file in files)
            {
                System.Console.WriteLine("{0} - Processing file: {1}", DateTime.Now, file);

                List<VLargeArray<int>> count = new List<VLargeArray<int>>();
                DNAMemoryStreamReader dnaStreamReader = new DNAMemoryStreamReader((MemoryDNA)MemoryDNA.LoadChromFile(file));

                for (int kmerLength = 1; kmerLength <= maxLength; kmerLength++)
                {
                    System.Console.WriteLine("{0} - Processing file: {1}, Len: {2}", DateTime.Now, file, kmerLength);

                    dnaStreamReader.Seek(0L);
                    count.Add(BuildKmerCount(dnaStreamReader, kmerLength));
                }

                string fileName = new FileInfo(file).Name;
                fileName = fileName.Substring(0, fileName.Length - 3);
                patternsCount.Add(fileName, count);
            }

            // Create output
            for (int kmerLength = 1; kmerLength <= maxLength; kmerLength++)
                using (StreamWriter sw = new StreamWriter(Path.Combine(folder, string.Format("KmerCount-{0}.txt", kmerLength))))
                {
                    System.Console.WriteLine("{0} - Creating output: {1}/{2}", DateTime.Now, kmerLength, maxLength);

                    List<VLargeArray<int>> patternsCounts = new List<VLargeArray<int>>();

                    // Header row
                    sw.Write('\t');
                    foreach (KeyValuePair<string, List<VLargeArray<int>>> kvp in patternsCount)
                    {
                        sw.Write(kvp.Key);
                        sw.Write('\t');
                        patternsCounts.Add(kvp.Value[kmerLength - 1]);
                    }
                    sw.WriteLine();

                    // Write rows
                    for (ulong i = 0; i < patternsCounts[0].Length; i++)
                    {
                        sw.Write(GetPatternStringFromValue(i, kmerLength));
                        sw.Write('\t');
                        for (int j = 0; j < patternsCounts.Count; j++)
                        {
                            sw.Write(patternsCounts[j][i]);
                            sw.Write('\t');
                        }
                        sw.WriteLine();
                    }
                }
        }

        /// <summary>
        /// Search for Kmer which are significantly less frequent than expected. For example
        /// CG known to appear way less then is should be.
        /// </summary>
        /// <param name="dnaFile"></param>
        static void CompareCountBetweenKmerLengths(string dnaFile, int maxPatternLength)
        {
            Dictionary<string, bool> knownPatterns = new Dictionary<string, bool>();

            SA.MemoryDNA memoryDna = (SA.MemoryDNA)SA.MemoryDNA.LoadChromFile(dnaFile);
            DNAMemoryStreamReader dnaStreamReader = new DNAMemoryStreamReader(memoryDna);

            List<Pair> pairs = new List<Pair>();
            pairs.Add(new Pair(0, memoryDna.Length));

            VLargeArray<int> patternCountPlusOne = BuildKmerCount(dnaStreamReader, pairs, 1);
            VLargeArray<int> patternCount;

            for (int patternLength = 1; patternLength <= maxPatternLength; patternLength++)
            {
                System.Console.WriteLine("Length: {0}", patternLength);

                patternCount = patternCountPlusOne;
                patternCountPlusOne = BuildKmerCount(dnaStreamReader, pairs, patternLength + 1);
                // Compare
                for (ulong patternCountPos = 0; patternCountPos < patternCount.Length; patternCountPos++)
                {
                    for (ulong newBase = 0; newBase < 4; newBase++)
                    {
                        ulong patternCountPlusOnePos = (patternCountPos << 2) + newBase;

                        int val = patternCount[patternCountPos];
                        int valPlusOne = patternCountPlusOne[patternCountPlusOnePos];

                        if (valPlusOne * 10 < val)
                        {
                            // Filter patterns which are extensions of patterns found
                            string patternPlusOneString = GetPatternStringFromValue(patternCountPlusOnePos, patternLength + 1);
                            bool foundDup = false;
                            for (int i = 1; i < patternPlusOneString.Length; i++)
                                if (knownPatterns.ContainsKey(patternPlusOneString.Substring(i)))
                                {
                                    foundDup = true;
                                    break;
                                }

                            if (foundDup == true)
                                break;

                            System.Console.WriteLine("Found! {0}: {1}, {2}: {3}",
                                GetPatternStringFromValue(patternCountPos, patternLength),
                                val,
                                patternPlusOneString,
                                valPlusOne
                                );

                            knownPatterns.Add(
                                patternPlusOneString,
                                true
                                );
                        }
                    }
                }
            }
        }


        static string StringExtender(string initString, int targetLength, Boolean extendWithIRs, int maxExtendInterval)
        {
            int count = 0;
            StringBuilder sb = new StringBuilder();
            while (sb.Length < targetLength)
            {
                int pos1 = gRandom.Next(initString.Length);
                int pos2 = gRandom.Next(initString.Length);
                if (pos1 > pos2)
                {
                    int tmp = pos1;
                    pos1 = pos2;
                    pos2 = tmp;
                }

                if (pos2 - pos1 > maxExtendInterval)
                    continue;

                string subscring = initString.Substring(pos1, pos2 - pos1 + 1);
                count++;
                if ((extendWithIRs == true) && ((count % 2) == 0))
                    subscring = (new SA.MemoryDNA(subscring)).GetComplementReverse().ToString();

                sb.Append(subscring);
            }

            return sb.ToString(0, targetLength);
        }

        static void ChromosomeFolding(string[] files, int chromosomeExtendLength, bool extendWithIRs, int[] numberOfFoldings, int[] foldingkength)
        {
            StringBuilder outputString = new StringBuilder();
            StringBuilder sbStatistics = new StringBuilder();
            StringBuilder sbStatisticsShort = new StringBuilder();

            foreach (string file in files)
            {
                //for (int maxExtendInterval = 1000; maxExtendInterval <= 16 * 1000; maxExtendInterval *= 2)
                int maxExtendInterval = int.MaxValue;
                {
                    FileInfo info = new FileInfo(file);
                    string name = info.Name;
                    if (chromosomeExtendLength > 0)
                        name += "_" + StringHelper.NumberToString(chromosomeExtendLength);
                    if (extendWithIRs == true)
                        name += "_WithIR";
                    if (maxExtendInterval < int.MaxValue)
                        name += "_ExtendMaxBlock" + StringHelper.NumberToString(maxExtendInterval);

                    ChromosomeFolding(outputString, sbStatistics, sbStatisticsShort, name, file, chromosomeExtendLength, maxExtendInterval, extendWithIRs, numberOfFoldings, foldingkength);
                }
            }

            FileInfo fileInfo = new FileInfo(files[0]);
            WriteOutput(fileInfo.DirectoryName, outputString, sbStatistics, sbStatisticsShort);
        }

        static void ChromosomeFolding(StringBuilder outputString, StringBuilder sbStatistics, StringBuilder sbStatisticsShort, string chromosome, string fileName, int chromosomeExtendLength, int maxExtendInterval, bool extendWithIRs, int[] numberOfFoldings, int[] foldingkength)
        {
            DateTime startTime = DateTime.Now;

            SA.MemoryDNA dna = (SA.MemoryDNA)SA.MemoryDNA.LoadChromFile(fileName);
            if (chromosomeExtendLength > 0)
                dna = new SA.MemoryDNA(StringExtender(dna.ToString(), chromosomeExtendLength, extendWithIRs, maxExtendInterval));

            dna.WriteToFile(@"c:\temp\chromosome_Folds0.fa");

            ProcessesChromosomeFile2_Stream(outputString, sbStatistics, sbStatisticsShort, new DNAMemoryStreamReader(dna), fileName, chromosome, null);
            int folds = 0;
            for (int j = 0; j < numberOfFoldings.Length; j++)
            {
                for (int i = 0; i < numberOfFoldings[j]; i++)
                {
                    folds++;
                    dna = PerformDNAFold(dna, foldingkength[j]);

//                    if (folds % 10 == 0)
//                        ProcessesChromosomeFile2_Stream(outputString, sbStatistics, new DNAMemoryStreamReader(dna), chromosome + "_Folds" + folds, null);
                }

                ProcessesChromosomeFile2_Stream(outputString, sbStatistics, sbStatisticsShort, new DNAMemoryStreamReader(dna), fileName, chromosome + "_Folds" + folds, null);
                dna.WriteToFile(@"c:\temp\chromosome_Folds" + folds + ".fa");
            }

            Console.WriteLine("Runtime: {0}", new TimeSpan(DateTime.Now.Ticks - startTime.Ticks));
        }

        static SA.MemoryDNA PerformDNAFold(SA.MemoryDNA dna)
        {
            return PerformDNAFold(dna, -1);
        }

        static SA.MemoryDNA PerformDNAFold(SA.MemoryDNA dna, int maxFoldingLength)
        {
            int pos1 = gRandom.Next(dna.Length);
            int pos2;
            if (maxFoldingLength < 0)
                pos2 = gRandom.Next(dna.Length);
            else
            {
                pos2 = pos1 + gRandom.Next(maxFoldingLength);
                while (pos2 >= dna.Length)
                    pos2 = pos1 + gRandom.Next(maxFoldingLength);
            }

            if (pos1 > pos2)
            {
                int tmp = pos1;
                pos1 = pos2;
                pos2 = tmp;
            }

            //string newDnaStr =
            //    dna.GetSubSequence(0, pos1).ToString() +
            //    dna.GetSubSequence(pos1, pos2 - pos1 + 1).GetComplementReverse().ToString() +
            //    dna.GetSubSequence(pos2 + 1, dna.Length - (pos2 + 1)).ToString();
            //
            //return new SA.MemoryDNA(newDnaStr);

            dna.PerformInplacePartialInvertedRepeat(pos1, pos2 - pos1 + 1);
            return dna;
        }

        static void ProcessFileMultipleRanges(string inputFile, List<Pair> ranges, string fileNameSuffix)
        {
            DateTime startTime = DateTime.Now;

            StringBuilder outputString = new StringBuilder();
            StringBuilder sbStatistics = new StringBuilder();
            StringBuilder sbStatisticsShort = new StringBuilder();

            FileInfo fileInfo = new FileInfo(inputFile);

            foreach (Pair range in ranges)
            {
                string name = fileInfo.Name + "_" + range.Start + "-" + range.End;

                ProcessesChromosomeFile2(outputString, sbStatistics, sbStatisticsShort, inputFile, name, new List<Pair>() { range });
            }

            WriteOutput(fileInfo.Directory.FullName, outputString, sbStatistics, sbStatisticsShort, fileNameSuffix);

            Console.WriteLine("Runtime: {0}", new TimeSpan(DateTime.Now.Ticks - startTime.Ticks));
        }

        static void ProcessFolder(string inputFolder)
        {
            DateTime startTime = DateTime.Now;

            StringBuilder outputString = new StringBuilder();
            StringBuilder sbStatistics = new StringBuilder();
            StringBuilder sbStatisticsShort = new StringBuilder();

            if (JOIN_STAT_FOR_MULTIPLE_FILES == true)
            {
                string[] files = Directory.GetFiles(inputFolder, INPUT_FILE_SEARCH_PATTERN);
                Array.Sort(files);
                DirectoryInfo info = new DirectoryInfo(inputFolder);
                string title = info.Name;
                ProcessesChromosomeMultipleFiles(outputString, sbStatistics, sbStatisticsShort, files, title);
            }
            else if (PROCESS_COVERAGE_FILE == true)
            {
                Dictionary<string, List<Pair>> coverageFile = loadCoverageFile();
                foreach (KeyValuePair<string, List<Pair>> kvp in coverageFile)
                {
                    Console.WriteLine("Processing: {0}", kvp.Key);

                    List<Pair> pairs = kvp.Value;
                    string fileName = Path.Combine(inputFolder, string.Format(INPUT_FILE_PATTERN, kvp.Key));
                    ProcessesChromosomeFile2(outputString, sbStatistics, sbStatisticsShort, fileName, kvp.Key, pairs);
                }
            }
            else // Process entire chromosomes
            {
                string[] files = Directory.GetFiles(inputFolder, INPUT_FILE_SEARCH_PATTERN);
                Array.Sort(files);

                foreach (string fileName in files)
                {
                    FileInfo fileInfo = new FileInfo(fileName);
                    string chromosome = fileInfo.Name.Substring(0, fileInfo.Name.Length - fileInfo.Extension.Length);
                    ProcessesChromosomeFile2(outputString, sbStatistics, sbStatisticsShort, fileName, chromosome, null);
                }
            }

            WriteOutput(inputFolder, outputString, sbStatistics, sbStatisticsShort);

            Console.WriteLine("Runtime: {0}", new TimeSpan(DateTime.Now.Ticks - startTime.Ticks));
        }

        private static void WriteOutput(string inputFolder, StringBuilder outputString, StringBuilder sbStatistics, StringBuilder sbStatisticsShort, string fileNameSuffix = "")
        {
            DirectoryInfo di = new DirectoryInfo(inputFolder);
            string outputFile = Path.Combine(inputFolder, "chargraff." + di.Name + fileNameSuffix);

            if (EXPORT_FULL_DATA == true)
            {
                using (StreamWriter sw = new StreamWriter(outputFile + ".txt"))
                {
                    sw.Write(outputString);
                }
            }

            using (StreamWriter sw = new StreamWriter(outputFile + ".stat.txt"))
            {
                sw.Write(sbStatistics);
            }

            using (StreamWriter sw = new StreamWriter(outputFile + ".stat.short.txt"))
            {
                sw.Write(sbStatisticsShort);
            }

        }

        private static void ProcessesChromosomeMultipleFiles(StringBuilder outputString, StringBuilder sbStatistics, StringBuilder sbStatisticsShort, string[] fileNames, string title)
        {
            Console.WriteLine("{0} - Starting to process {1} - loading...", DateTime.Now, title);

            for (int patternLength = 1; patternLength <= MAX_SEQUENCE_LENGTH; patternLength++)
            {
                Console.WriteLine("{0} - Starting to process {1} length {2}", DateTime.Now, title, patternLength);
                //int[] patternCount = new int[(int)Math.Pow(4, patternLength)];
                //Array.Clear(patternCount, 0, patternCount.Length);
                VLargeArray<int> patternCount = VLargeArray<int>.AllocateArray((ulong)Math.Pow(4, patternLength));

                ulong patternMask = 0;
                for (int j = 0; j < patternLength; j++)
                    patternMask = (patternMask << 2) + 0x3;

                foreach (string fileName in fileNames)
                {
                    DNAStreamReader dnaStreamReader = new DNAFileStreamReader(fileName);

                    // Read first few chars
                    ulong patternVal = 0;
                    ulong patternInvalid = 0;
                    for (int j = 0; j < patternLength - 1; j++)
                        ProcessesChromosomeFile2_ReadCharToMask(dnaStreamReader, patternMask, ref patternVal, ref patternInvalid);

                    while (true)
                    {
                        bool foundChar = ProcessesChromosomeFile2_ReadCharToMask(dnaStreamReader, patternMask, ref patternVal, ref patternInvalid);
                        if (foundChar == false)
                            break;

                        if (patternInvalid > 0)
                            continue;

                        patternCount[patternVal]++;
                    }
                }

                // Print full output
                if (EXPORT_FULL_DATA == true)
                {
                    for (ulong i = 0; i < patternCount.Length; i++)
                    {
                        outputString.AppendLine(string.Format("{0} {1} {2} {3}", title, patternLength, GetPatternStringFromValue((ulong)i, patternLength), patternCount[i]));
                    }
                }

                GenerateStatistics(sbStatistics, sbStatisticsShort, null, title, patternLength, patternCount);
            }
        }

        private static void ProcessesChromosome_CompareTwoFiles(List<CompareGenome> compares)
        {
            DateTime startTime = DateTime.Now;

            StringBuilder outputString = new StringBuilder();
            StringBuilder sbStatistics = new StringBuilder();
            StringBuilder sbStatisticsShort = new StringBuilder();

            foreach (CompareGenome comp in compares)
                ProcessesChromosome_CompareTwoFiles(outputString, sbStatistics, sbStatisticsShort, comp.Path1, comp.Path2, comp.Title);

            FileInfo fileInfo = new FileInfo(compares[0].Path1);
            WriteOutput(fileInfo.DirectoryName, outputString, sbStatistics, sbStatisticsShort);

            Console.WriteLine("Runtime: {0}", new TimeSpan(DateTime.Now.Ticks - startTime.Ticks));
        }

        private static void ProcessesChromosome_CompareTwoFiles(StringBuilder outputString, StringBuilder sbStatistics, StringBuilder sbStatisticsShort, string fileName1, string fileName2, string title)
        {
            Console.WriteLine("{0} - Starting to process {1} - loading...", DateTime.Now, title);
            DNAStreamReader dnaStreamReader1 = new DNAFileStreamReader(fileName1);
            DNAStreamReader dnaStreamReader2 = new DNAFileStreamReader(fileName2);

            for (int patternLength = 1; patternLength <= MAX_SEQUENCE_LENGTH; patternLength++)
            {
                Console.WriteLine("{0} - Starting to process {1} length {2}", DateTime.Now, title, patternLength);
                VLargeArray<double> patternCount1 = BuildKmerPercent(dnaStreamReader1, patternLength);
                VLargeArray<double> patternCount2 = BuildKmerPercent(dnaStreamReader2, patternLength);

                GenerateStatisticsComp(sbStatistics, title, patternLength, patternCount1, patternCount2);
            }
        }

        private static void ProcessesChromosomeFile2(StringBuilder outputString, StringBuilder sbStatistics, StringBuilder sbStatisticsShort, string fileName, string chromosome, List<Pair> pairs)
        {
            Console.WriteLine("{0} - Starting to process {1} - loading...", DateTime.Now, chromosome);
            DNAStreamReader dnaStreamReader = new DNAFileStreamReader(fileName);

            ProcessesChromosomeFile2_Stream(outputString, sbStatistics, sbStatisticsShort, dnaStreamReader, fileName, chromosome, pairs);
        }

        private static void ProcessesChromosomeFile2_Stream(StringBuilder outputString, StringBuilder sbStatistics, StringBuilder sbStatisticsShort, DNAStreamReader dnaStreamReader, string fileName, string chromosome, List<Pair> pairs)
        {
            string simpleFolderName = new FileInfo(fileName).Directory.Name;

            if (pairs == null)
            {
                pairs = new List<Pair>();
                pairs.Add(new Pair(0, dnaStreamReader.Length));
            }

            for (int patternLength = 1; patternLength <= MAX_SEQUENCE_LENGTH; patternLength++)
            {
                Console.WriteLine("{0} - Starting to process {1} length {2}", DateTime.Now, chromosome, patternLength);
                //int[] patternCount = new int[(int)Math.Pow(4, patternLength)];
                //Array.Clear(patternCount, 0, patternCount.Length);
                VLargeArray<int> patternCount = BuildKmerCount(dnaStreamReader, pairs, patternLength);

                // Print full output
                if (EXPORT_FULL_DATA == true)
                {
                    ulong patternMask = 0;
                    for (int j = 0; j < patternLength; j++)
                        patternMask = (patternMask << 2) + 0x3;                

                    for (ulong i = 0; i < patternCount.Length; i++)
                    {
                        ulong irPattern = ConvertPattern_BuildIR(patternLength, patternMask, i);
                        ulong reversedPattern = ConvertPattern_BuildReversed(patternLength, patternMask, i);
                        ulong invertedPattern = ConvertPattern_BuildInverted(patternLength, patternMask, i);

                        if (APPEND_SPECIMEN_NAME == true)
                        {
                            outputString.Append(simpleFolderName);
                            outputString.Append('\t');
                        }

                        outputString.AppendLine(string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}", 
                            chromosome, patternLength, 
                            GetPatternStringFromValue(i, patternLength), patternCount[i], 
                            GetPatternStringFromValue(irPattern, patternLength), patternCount[irPattern],
                            GetPatternStringFromValue(reversedPattern, patternLength), patternCount[reversedPattern],
                            GetPatternStringFromValue(invertedPattern, patternLength), patternCount[invertedPattern]
                            ));
                    }
                }

                GenerateStatistics(sbStatistics, sbStatisticsShort, simpleFolderName, chromosome, patternLength, patternCount);
            }
        }

        private static VLargeArray<double> BuildKmerPercent(DNAStreamReader dnaStreamReader, int patternLength)
        {
            VLargeArray<int> count = BuildKmerCount(dnaStreamReader, patternLength);
            long streamLength = 0;
            for (ulong i = 0; i < count.Length; i++)
                streamLength += count[i];

            VLargeArray<double> countPercent = VLargeArray<double>.AllocateArray(count.Length);
            for (ulong i = 0; i < countPercent.Length; i++)
                countPercent[i] = (double)count[i] / (double)streamLength;

            return countPercent;
        }

        private static VLargeArray<int> BuildKmerCount(DNAStreamReader dnaStreamReader, int patternLength)
        {
            List<Pair> pairs = new List<Pair>();
            pairs.Add(new Pair(0, dnaStreamReader.Length));

            return BuildKmerCount(dnaStreamReader, pairs, patternLength);
        }

        private static VLargeArray<int> BuildKmerCount(DNAStreamReader dnaStreamReader, List<Pair> pairs, int patternLength)
        {
            VLargeArray<int> patternCount = VLargeArray<int>.AllocateArray((ulong)Math.Pow(4, patternLength));

            ulong patternMask = 0;
            for (int j = 0; j < patternLength; j++)
                patternMask = (patternMask << 2) + 0x3;

            foreach (Pair pair in pairs)
            {
                dnaStreamReader.Seek(pair.Start);

                // Read first few chars
                ulong patternVal = 0;
                ulong patternInvalid = 0;
                for (int j = 0; j < patternLength - 1; j++)
                    ProcessesChromosomeFile2_ReadCharToMask(dnaStreamReader, patternMask, ref patternVal, ref patternInvalid);

                long pairEnd = Math.Min(pair.End, dnaStreamReader.Length);
                for (long i = pair.Start; i < pairEnd - (patternLength - 1); i++)
                {
                    ProcessesChromosomeFile2_ReadCharToMask(dnaStreamReader, patternMask, ref patternVal, ref patternInvalid);
                    if (patternInvalid > 0)
                        continue;

                    patternCount[patternVal]++;
                }
            }
            return patternCount;
        }

        private static void GenerateStatisticsComp(StringBuilder sbStatistics, string title, int patternLength, VLargeArray<double> patternPercent1, VLargeArray<double> patternPercent2)
        {
            StatBuilder n1n2AbsDivStat = new StatBuilder("Abs(N1-N2)/(N1+N2)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // (N1-N2)/(N1+N2)
            //StatBuilder n1n2DivStat = new StatBuilder("(N1-N2)/(N1+N2)"); // (N1-N2)/(N1+N2)
            StatBuilder n1n2irAbsDivStat = new StatBuilder("Abs(N1-N2ir)/(N1+N2ir)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // (N1-N2)/(N1+N2)
            StatBuilder n1n2AbsStat = new StatBuilder("Abs(N1 - N2)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS);
            StatBuilder n1n2irRootAbsDivStat = new StatBuilder("Abs(N1-N2ir)/Sqrt(N1+N2ir)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // (N1-N2)/(N1+N2)

            StatBuilder[] statBuilders = new StatBuilder[] { 
                //n1n2DivStat, 
                n1n2AbsDivStat, 
                n1n2irAbsDivStat,
                n1n2AbsStat,
                n1n2irRootAbsDivStat
            };

            ulong patternMask = 0;
            for (int j = 0; j < patternLength; j++)
                patternMask = (patternMask << 2) + 0x3;

            for (ulong i = 0; i < (ulong)patternPercent1.Length; i++)
            {
                ulong irPattern = ConvertPattern_BuildIR(patternLength, patternMask, i);

                double p1 = patternPercent1[i];
                double p2 = patternPercent2[i];
                double p2ir = patternPercent2[irPattern];

                n1n2AbsDivStat.IncreaseZeroCounts(p1, p2);
                //n1n2DivStat.IncreaseZeroCounts(p1, p2);
                n1n2irAbsDivStat.IncreaseZeroCounts(p1, p2ir);
                n1n2irRootAbsDivStat.IncreaseZeroCounts(p1, p2ir);
                n1n2AbsStat.IncreaseZeroCounts(p1, p2);

                // (N1-N2)/(N1+N2)
                if ((p1 > 0) || (p2 > 0))
                {
                    n1n2AbsDivStat.AddValue(Math.Abs((p1 - p2) / (p1 + p2)));
                    //n1n2DivStat.AddValue((p1 - p2) / (p1 + p2));
                }

                if ((p1 > 0) || (p2ir > 0))
                {
                    n1n2irAbsDivStat.AddValue(Math.Abs((p1 - p2ir) / (p1 + p2ir)));
                    n1n2irRootAbsDivStat.AddValue(Math.Abs((p1 - p2ir) / Math.Sqrt(p1 + p2ir)));
                }

                n1n2AbsStat.AddValue(Math.Abs(p1 - p2));
            }

            foreach (StatBuilder statBuilder in statBuilders)
                GenerateStatistics_BuildStatRow(sbStatistics, null, title, patternLength, 1, statBuilder, PRINT_HISTOGRAM);
        }

        private static void GenerateStatistics(StringBuilder sbStatistics, StringBuilder sbStatisticsShort, string specimen, string title, int patternLength, VLargeArray<int> patternCount)
        {
            // Generate statistics
            string minPattern = null;
            double minPercent = double.MaxValue;
            int minValue = 0;
            string maxPattern = null;
            double maxPercent = double.MinValue;
            int maxValue = int.MaxValue;

            //StatBuilder patternDiffPercentStat = new StatBuilder(); // N1/N2
            //StatBuilder pairCountStat = new StatBuilder(); // N1+N2
            StatBuilder baseCountStat = new StatBuilder("N1 Count", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // N1
            StatBuilder n1n2DivStat = new StatBuilder("|N1-Nir|/(N1+Nir)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // (N1-N2)/(N1+N2)
            StatBuilder n1n2SqrtDivStat = new StatBuilder("|N1-Nir|/Sqrt(N1+Nir)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS);
            StatBuilder n1RandomN1DivStat = new StatBuilder("|N1 - Nrand|/(N1 + Nrand)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // (N1 - Nr)/(N1 + Nr)
            StatBuilder n1ReversedN1DivStat = new StatBuilder("|N1 - Nrev|/(N1 + Nrev)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // (N1 - Nr)/(N1 + Nr)
            StatBuilder n1ReversedN1SqrtDivStat = new StatBuilder("|N1 - Nrev|/Sqrt(N1 + Nrev)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // (N1 - Nr)/(N1 + Nr)
            StatBuilder n1InvertedN1DivStat = new StatBuilder("|N1 - Ninv|/(N1 + Ninv)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS);

            StatBuilder n1Minusn2Stat = new StatBuilder("Abs(N1-Nir)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // Abs(N1-N2)
            StatBuilder n1Plusn2Stat = new StatBuilder("(N1+Nir)", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // (N1+N2)

            StatBuilder n1n2DivNoCgStat = new StatBuilder("(N1-Nir)/(N1+Nir), No CG", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // (N1-N2)/(N1+N2)
            StatBuilder n1RandomN1DivNoCgStat = new StatBuilder("(N1 - Nrand)/(N1 + Nrand), No CG", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // (N1 - Nr)/(N1 + Nr)
            StatBuilder n1ReversedN1DivNoCgStat = new StatBuilder("(N1 - Nrev)/(N1 + Nrev), No CG", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS); // (N1 - Nr)/(N1 + Nr)
            StatBuilder n1InvertedN1DivNoCgStat = new StatBuilder("(N1 - Ninv)/(N1 + Ninv), No CG", PRINT_HISTOGRAM, HISTOGRAM_BUCKETS);
            
            //double countG = patternCount[(ulong)GetIntValue('G')];
            //double countC = patternCount[(ulong)GetIntValue('C')];
            //double countA = patternCount[(ulong)GetIntValue('A')];
            //double countT = patternCount[(ulong)GetIntValue('T')];
            //StatBuilder skewStat = new StatBuilder("(#G-#C)/(#G+#C) + (#T-#A)/(#T+#A)");
            //if ((countC + countG > 0)&&(countA + countT > 0))
            //    skewStat.AddValue(
            //        (countG - countC) / (countG + countC) +
            //        (countT - countA) / (countT + countA)
            //        );

            StatBuilder[] statBuilders = new StatBuilder[] { 
                n1n2DivStat, n1n2SqrtDivStat, n1RandomN1DivStat, n1ReversedN1DivStat, n1ReversedN1SqrtDivStat, n1InvertedN1DivStat,
                n1n2DivNoCgStat, n1RandomN1DivNoCgStat, n1ReversedN1DivNoCgStat, n1InvertedN1DivNoCgStat 
                , n1Minusn2Stat, n1Plusn2Stat
            //    , skewStat
            };

            ulong patternMask = 0;
            for (int j = 0; j < patternLength; j++)
                patternMask = (patternMask << 2) + 0x3;

            // We are processing most of the pairs: Pattern-IR twice, however, there are some patterns
            // Forwhich the Pattern == IR (for example 'AT'), so we can't eliminate and double by 2.
            for (ulong i = 0; i < (ulong)patternCount.Length; i++)
            {
                // Build bit-wise IR. Assumption: every A,C,G,T bits: A == ~T , C == ~G
                ulong irPattern = ConvertPattern_BuildIR(patternLength, patternMask, i);

                int patternCountI = patternCount[i];
                int irCount = patternCount[irPattern];

                string highPattern;
                int highValue;
                string lowPattern;
                int lowValue;

                if (irCount < patternCountI)
                {
                    lowValue = irCount;
                    lowPattern = GetPatternStringFromValue(irPattern, patternLength); ;
                    highValue = patternCountI;
                    highPattern = GetPatternStringFromValue(i, patternLength);
                }
                else
                {
                    highValue = irCount;
                    highPattern = GetPatternStringFromValue(irPattern, patternLength);
                    lowValue = patternCountI;
                    lowPattern = GetPatternStringFromValue(i, patternLength);
                }

                if (highValue > 0)
                {
                    double percent = (double)lowValue / (double)highValue;

                    if ((percent < minPercent)||
                        ((percent == maxPercent) && (minValue < highValue)))
                    {
                        minPercent = percent;
                        minValue = highValue;
                        minPattern = highPattern;
                    }

                    if ((percent > maxPercent)||
                        ((percent == maxPercent)&&(maxValue < highValue)))
                    {
                        maxPercent = percent;
                        maxValue = highValue;
                        maxPattern = highPattern;
                    }

                    //patternDiffPercentStat.AddValue(percent);
                }

                n1Minusn2Stat.AddValue(Math.Abs((double)(irCount - patternCountI)));
                n1Plusn2Stat.AddValue((double)(irCount + patternCountI));

                // IR (N1-N2)/(N1+N2)
                if ((irCount > 0) || (patternCountI > 0))
                {
                    n1n2DivStat.AddValue(Math.Abs((double)(irCount - patternCountI) / (double)(irCount + patternCountI)));
                    n1n2SqrtDivStat.AddValue(Math.Abs((double)(irCount - patternCountI) / Math.Sqrt((double)(irCount + patternCountI))));
                }

                // Inverted (N1-N2)/(N1+N2)
                ulong patternInverted = ConvertPattern_BuildInverted(patternLength, patternMask, i);
                int invertedN2value = patternCount[patternInverted];
                if ((invertedN2value > 0) || (patternCountI > 0))
                    n1InvertedN1DivStat.AddValue(Math.Abs((double)(invertedN2value - patternCountI) / (double)(invertedN2value + patternCountI)));

                // Random (N1-N2)/(N1+N2)
                ulong randomN2 = (ulong)(gRandom.NextDouble() * (double)patternCount.Length);
                int randomN2value = patternCount[randomN2];
                if ((randomN2value > 0) || (patternCountI > 0))
                    n1RandomN1DivStat.AddValue(Math.Abs((double)(randomN2value - patternCountI) / (double)(randomN2value + patternCountI)));

                // Reversed (N1 - N2)/(N1 + N2)
                ulong patternReversed = ConvertPattern_BuildReversed(patternLength, patternMask, i);
                int reversedN2value = patternCount[patternReversed];
                if ((reversedN2value > 0) || (patternCountI > 0))
                {
                    n1ReversedN1DivStat.AddValue(Math.Abs((double)(reversedN2value - patternCountI) / (double)(reversedN2value + patternCountI)));
                    n1ReversedN1SqrtDivStat.AddValue(Math.Abs((double)(reversedN2value - patternCountI) / Math.Sqrt((double)(reversedN2value + patternCountI))));
                }
                n1n2DivStat.IncreaseZeroCounts(irCount, patternCountI);
                n1n2SqrtDivStat.IncreaseZeroCounts(irCount, patternCountI);
                n1InvertedN1DivStat.IncreaseZeroCounts(invertedN2value, patternCountI);
                n1RandomN1DivStat.IncreaseZeroCounts(randomN2value, patternCountI);
                n1ReversedN1DivStat.IncreaseZeroCounts(reversedN2value, patternCountI);
                n1ReversedN1SqrtDivStat.IncreaseZeroCounts(reversedN2value, patternCountI);

                if (CheckIfContainsPatternCGorGC(i) == false)
                {
                    if ((irCount > 0) || (patternCountI > 0))
                        n1n2DivNoCgStat.AddValue(Math.Abs((double)(irCount - patternCountI) / (double)(irCount + patternCountI)));

                    if ((invertedN2value > 0) || (patternCountI > 0))
                        n1InvertedN1DivNoCgStat.AddValue(Math.Abs((double)(invertedN2value - patternCountI) / (double)(invertedN2value + patternCountI)));

                    if (CheckIfContainsPatternCGorGC(randomN2) == false)
                        if ((randomN2value > 0) || (patternCountI > 0))
                            n1RandomN1DivNoCgStat.AddValue(Math.Abs((double)(randomN2value - patternCountI) / (double)(randomN2value + patternCountI)));

                    if ((reversedN2value > 0) || (patternCountI > 0))
                        n1ReversedN1DivNoCgStat.AddValue(Math.Abs((double)(reversedN2value - patternCountI) / (double)(reversedN2value + patternCountI)));

                    n1n2DivNoCgStat.IncreaseZeroCounts(irCount, patternCountI);
                    n1InvertedN1DivNoCgStat.IncreaseZeroCounts(invertedN2value, patternCountI);
                    n1RandomN1DivNoCgStat.IncreaseZeroCounts(randomN2value, patternCountI);
                    n1ReversedN1DivNoCgStat.IncreaseZeroCounts(reversedN2value, patternCountI);
                }

                //pairCountStat.AddValue(lowValue + highValue);
                baseCountStat.AddValue(patternCountI);
            }

            foreach (StatBuilder statBuilder in statBuilders)
                GenerateStatistics_BuildStatRow(sbStatistics, specimen, title, patternLength, baseCountStat.GetSum(), statBuilder, PRINT_HISTOGRAM);

            GenerateStatistics_BuildStatRow(sbStatisticsShort, specimen, title, patternLength, baseCountStat.GetSum(), n1n2DivStat, false);
        }

        private static void GenerateStatistics_BuildStatRow(
            StringBuilder sbStatistics, string specimen, string title, int patternLength,
            double patternCount, StatBuilder stat, Boolean printHistogram)
        {
            if (specimen != null)
            {
                sbStatistics.Append(specimen);
                sbStatistics.Append('\t');
            }

            sbStatistics.Append(title);
            sbStatistics.Append('\t');
            sbStatistics.Append(patternLength);

            // Processed chromosome length
            sbStatistics.Append('\t');
            sbStatistics.Append(patternCount);

            /*
            sbStatistics.Append('\t');
            sbStatistics.Append(minPattern);
            sbStatistics.Append('\t');
            sbStatistics.Append(minPercent);
            sbStatistics.Append('\t');
            sbStatistics.Append(minValue);
            sbStatistics.Append('\t');
            sbStatistics.Append(maxPattern);
            sbStatistics.Append('\t');
            sbStatistics.Append(maxPercent);
            sbStatistics.Append('\t');
            sbStatistics.Append(maxValue);
            */

            // (N1-N2)/(N1+N2)
            sbStatistics.Append('\t');
            sbStatistics.Append(stat.Title);

            sbStatistics.Append('\t');
            sbStatistics.Append(stat.CountOneZero);
            sbStatistics.Append('\t');
            sbStatistics.Append(stat.CountBothZero);

            sbStatistics.Append('\t');
            sbStatistics.Append(stat.GetCount());
            sbStatistics.Append('\t');
            sbStatistics.Append(stat.GetAverage());
            sbStatistics.Append('\t');
            sbStatistics.Append(stat.GetVariance());

            if (printHistogram == true)
            {
                foreach (int histogramVal in stat.GetHistogram())
                {
                    sbStatistics.Append('\t');
                    sbStatistics.Append((double)histogramVal / (double)stat.GetCount());
                }
            }

            sbStatistics.AppendLine();
        }

        private static ulong ConvertPattern_BuildIR(int patternLength, ulong patternMask, ulong pattern)
        {
            ulong irPattern = 0;
            ulong tmpPattern = pattern;
            for (int j = 0; j < patternLength; j++)
            {
                irPattern = (irPattern << 2) + (tmpPattern & 0x3);
                tmpPattern = tmpPattern >> 2;
            }
            irPattern = (~irPattern) & patternMask;

            return irPattern;
        }

        private static ulong ConvertPattern_BuildInverted(int patternLength, ulong patternMask, ulong pattern)
        {
            return (~pattern) & patternMask;
        }

        private static ulong ConvertPattern_BuildReversed(int patternLength, ulong patternMask, ulong pattern)
        {
            ulong irPattern = 0;
            ulong tmpPattern = pattern;

            for (int j = 0; j < patternLength; j++)
            {
                irPattern = (irPattern << 2) + (tmpPattern & 0x3);
                tmpPattern = tmpPattern >> 2;
            }

            return irPattern;
        }

        private static bool ProcessesChromosomeFile2_ReadCharToMask(DNAStreamReader dnaStreamReader, ulong patternMask, ref ulong patternVal, ref ulong patternInvalid)
        {
            char ch = dnaStreamReader.Read();
            if (ch == DNAStreamReader.EOF)
                return false;

            ulong chVal = (ulong)GetIntValue(ch);
            if (chVal > 3)
            {
                patternVal = (patternVal << 2);
                patternInvalid = (patternInvalid << 2) + 0x1;
            }
            else
            {
                patternVal = (patternVal << 2) + chVal;
                patternInvalid = (patternInvalid << 2);
            }

            patternVal = patternVal & patternMask;
            patternInvalid = patternInvalid & patternMask;

            return true;
        }

        private static void ProcessesChromosomeFile(StringBuilder outputString, StringBuilder sbStatistics, string fileName, string chromosome, List<Pair> pairs)
        {
            Console.WriteLine("{0} - Starting to process {1} - loading...", DateTime.Now, chromosome);

            SA.MemoryDNA memoryDNA = (SA.MemoryDNA)SA.MemoryDNA.LoadChromFile(fileName);
            char[] dna = new char[memoryDNA.Length];
            for (int i = 0; i < dna.Length; i++)
            {
                char c = memoryDNA[i];
                if ((c != 'A') && (c != 'T') && (c != 'C') && (c != 'G'))
                    c = ' ';
                dna[i] = c;
            }
            memoryDNA = null;

            if (pairs == null)
            {
                pairs = new List<Pair>();
                pairs.Add(new Pair(0, dna.Length));
            }

            for (int patternLength = 1; patternLength <= MAX_SEQUENCE_LENGTH; patternLength++)
            {
                Console.WriteLine("{0} - Starting to process {1} length {2}", DateTime.Now, chromosome, patternLength);
                Dictionary<string, IntContainer> patternCount = new Dictionary<string, IntContainer>();

                foreach (Pair pair in pairs)
                {
                    long pairEnd = Math.Min(pair.End, dna.Length);
                    for (long i = pair.Start; i < pairEnd - (patternLength - 1); i++)
                    {
                        bool foundInvalidBase = false;
                        for (int j = 0; j < patternLength; j++)
                            if (dna[i + j] == ' ')
                            {
                                foundInvalidBase = true;
                                break;
                            }

                        if (foundInvalidBase == true)
                            continue;

                        string pattern = new String(dna, (int)i, patternLength);

                        IntContainer count;
                        if (patternCount.TryGetValue(pattern, out count) == false)
                        {
                            count = new IntContainer();
                            patternCount.Add(pattern, count);
                        }
                        count.Value++;
                    }
                }

                List<string> sortedKeys = new List<string>(patternCount.Keys);
                sortedKeys.Sort();

                // Print full output
                if (EXPORT_FULL_DATA == true)
                {
                    for (int i = 0; i < sortedKeys.Count; i++)
                        outputString.AppendLine(string.Format("{0} {1} {2} {3}", chromosome, patternLength, sortedKeys[i], patternCount[sortedKeys[i]].Value));
                }

                // Generate statistics
                string minPattern = null;
                double minPercent = double.MaxValue;
                int minValue = 0;
                string maxPattern = null;
                double maxPercent = double.MinValue;
                int maxValue = int.MaxValue;

                long n = 0;
                double percentSum = 0;
                double percentSequareSum = 0;

                // We are processing most of the pairs: Pattern-IR twice, however, there are some patterns
                // Forwhich the Pattern == IR (for example 'AT'), so we can't eliminate and double by 2.
                for (int k = 0; k < sortedKeys.Count; k++)
                {
                    string ir = (new SA.MemoryDNA(sortedKeys[k])).GetComplementReverse().ToString();

                    IntContainer irCount;
                    if (patternCount.TryGetValue(ir, out irCount) == false)
                        irCount = new IntContainer();

                    string highPattern;
                    int highValue;
                    string lowPattern;
                    int lowValue;

                    if (irCount.Value < patternCount[sortedKeys[k]].Value)
                    {
                        lowValue = irCount.Value;
                        lowPattern = ir;
                        highValue = patternCount[sortedKeys[k]].Value;
                        highPattern = sortedKeys[k];
                    }
                    else
                    {
                        highValue = irCount.Value;
                        highPattern = ir;
                        lowValue = patternCount[sortedKeys[k]].Value;
                        lowPattern = sortedKeys[k];
                    }

                    double percent = (double)lowValue / (double)highValue;
                    if (percent < minPercent)
                    {
                        minPercent = percent;
                        minValue = lowValue;
                        minPattern = lowPattern;
                    }
                    if (percent > maxPercent)
                    {
                        maxPercent = percent;
                        maxValue = lowValue;
                        maxPattern = lowPattern;
                    }

                    n++;
                    percentSum += percent;
                    percentSequareSum += (percent * percent);
                }

                double avg = percentSum / (double)n;
                double var = (percentSequareSum / (double)n) - (avg * avg);

                sbStatistics.Append(chromosome);
                sbStatistics.Append('\t');
                sbStatistics.Append(patternLength);
                sbStatistics.Append('\t');
                sbStatistics.Append(n);
                sbStatistics.Append('\t');
                sbStatistics.Append(avg);
                sbStatistics.Append('\t');
                sbStatistics.Append(var);
                sbStatistics.Append('\t');
                sbStatistics.Append(minPattern);
                sbStatistics.Append('\t');
                sbStatistics.Append(minPercent);
                sbStatistics.Append('\t');
                sbStatistics.Append(minValue);
                sbStatistics.Append('\t');
                sbStatistics.Append(maxPattern);
                sbStatistics.Append('\t');
                sbStatistics.Append(maxPercent);
                sbStatistics.Append('\t');
                sbStatistics.Append(maxValue);

                sbStatistics.AppendLine();
            }
        }

        private static Dictionary<string, List<Pair>> loadCoverageFile()
        {
            Dictionary<string, List<Pair>> coverageFile = new Dictionary<string, List<Pair>>();
            using (StreamReader sr = new StreamReader(INPUT_COVERAGE_FILE))
            {
                List<Pair> pairs = new List<Pair>();
                string currentChr = null;
                string line;
                while ((line = sr.ReadLine()) != null)
                {
                    string[] parts = line.Split(' ');
                    if (parts.Length != 3)
                        continue;

                    if (parts[0] != currentChr)
                    {
                        if (pairs.Count > 0)
                            coverageFile.Add(currentChr, pairs);

                        pairs = new List<Pair>();
                        currentChr = parts[0];
                    }

                    pairs.Add(new Pair(long.Parse(parts[1]), long.Parse(parts[2])));
                }

                if (currentChr != null)
                    coverageFile.Add(currentChr, pairs);
            }

            return coverageFile;
        }

        private static List<string> buildAllCombinations(char[] letters, int length)
        {
            List<string> combinations = new List<string>();
            buildAllCombinations_Recursive(combinations, letters, length, "");

            return combinations;
        }

        private static void buildAllCombinations_Recursive(List<string> patterns, char[] letters, int length, string currentPattern)
        {
            if (currentPattern.Length == length)
            {
                patterns.Add(currentPattern);
                return;
            }

            for (int i = 0; i < letters.Length; i++)
                buildAllCombinations_Recursive(patterns, letters, length, currentPattern + letters[i].ToString());
        }

        // Bit pattern to string
        private static String GetPatternStringFromValue(ulong val, int patternLength)
        {
            StringBuilder sb = new StringBuilder(patternLength);            
            for (int j = 0; j < patternLength; j++)
            {
                char ch = GetCharValue((int)(val & 0x3));
                sb.Insert(0, ch);
                val = val >> 2;
            }

            return sb.ToString();
        }

        private static char GetCharValue(int val)
        {
            switch (val)
            {
                case 0: return 'A';
                case 1: return 'C';
                case 2: return 'G';
                case 3: return 'T';
            }

            throw new Exception("Unsupported value: " + val);
        }

        private static int GetIntValue(char ch)
        {
            switch (ch)
            {
                case 'A':
                case 'a':
                    return 0;
                case 'C':
                case 'c':
                    return 1;
                case 'G':
                case 'g':
                    return 2;
                case 'T':
                case 't':
                    return 3;
                default:
                    return 666;
            }
        }

        private static readonly ulong PATTERN_CG = (ulong)(GetIntValue('C') << 2 | GetIntValue('G'));
        private static readonly ulong PATTERN_GC = (ulong)(GetIntValue('G') << 2 | GetIntValue('C'));

        private static bool CheckIfContainsPatternCGorGC(ulong pattern)
        {
            return 
                CheckIfContainsPattern(pattern, PATTERN_CG) ||
                CheckIfContainsPattern(pattern, PATTERN_GC);
        }

        private static bool CheckIfContainsPattern(ulong pattern, ulong searchPattern)
        {
            if (searchPattern == 0)
                throw new Exception("Search pattern can not be 0");

            while (pattern > 0)
            {
                if ((pattern & searchPattern) == searchPattern)
                    return true;

                pattern >>= 2;
            }

            return false;
        }
    }

    class IntContainer
    {
        public int Value;

        public IntContainer()
        {
            Value = 0;
        }
    }

    struct Pair
    {
        public long Start;
        public long End;

        public Pair(long start, long end)
        {
            Start = start;
            End = end;
        }
    }

    class CompareGenome
    {
        public string Path1;
        public string Path2;
        public string Title;

        public CompareGenome(string path1, string path2, string title)
        {
            this.Path1 = path1;
            this.Path2 = path2;
            this.Title = title;
        }
    }
}
