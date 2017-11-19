using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace SA
{
    public abstract class DNAStreamReader : IDisposable
    {
        public const char EOF = '\0';

        public abstract long Length
        {
            get;
        }

        public abstract void Seek(long position);

        public abstract char Read();

        public virtual void Dispose()
        {
        }
    }

    public class DNAMemoryStreamReader : DNAStreamReader
    {
        private SA.MemoryDNA mMemoryDNA;
        private int mPosition;

        public override long Length
        {
            get { return mMemoryDNA.Length; }
        }

        public DNAMemoryStreamReader(SA.MemoryDNA memoryDNA)
        {
            mMemoryDNA = memoryDNA;
            mPosition = 0;
        }

        public override void Seek(long position)
        {
            mPosition = (int)position;
        }

        public override char Read()
        {
            if (mPosition >= this.Length)
                return EOF;
            char ch = mMemoryDNA[mPosition];

            mPosition++;

            return ch;
        }
    }

    public class DNAFileStreamReader : DNAStreamReader
    {
        private string mPath;
        private StreamReader mStreamReader;
        private long mLength;
        private long mCurrentPosition;

        string mLineBuffer;
        int mLineBufferNextPosition;

        public override long Length
        {
            get 
            {
                if (mLength == -1)
                {
                    long oldPos = mCurrentPosition;
                    mLength = 0;
                    char ch;
                    while ((ch = Read()) != EOF)
                        mLength++;
                    Seek(oldPos);
                }

                return mLength; 
            }
        }

        public DNAFileStreamReader(string path)
        {
            mPath = path;

            Reset();
            mLength = -1;
        }

        public override void Seek(long position)
        {
            Reset();
            for (long i = 0; i < position; i++)
                Read();
        }

        public void Reset()
        {
            mCurrentPosition = 0;
            mStreamReader = new StreamReader(mPath);
            ReadLine();
        }

        private void ReadLine()
        {
            while ((mLineBuffer = mStreamReader.ReadLine()) != null)
            {
                if (mLineBuffer.StartsWith(">") == false)
                    break;
            }

            mLineBufferNextPosition = 0;
        }

        public override char Read()
        {
            if (mLineBuffer == null)
                return EOF;

            while (mLineBufferNextPosition >= mLineBuffer.Length)
            {
                ReadLine();
                if (mLineBuffer == null)
                    return EOF;
            }

            char ch = mLineBuffer[mLineBufferNextPosition];
            mLineBufferNextPosition++;
            mCurrentPosition++;

            return ch;
        }

        public void Close()
        {
            mStreamReader.Close();
        }

        public override void Dispose()
        {
            base.Dispose();

            if (mStreamReader != null)
            {
                mStreamReader.Dispose();
                mStreamReader = null;
            }
        }
    }
}
