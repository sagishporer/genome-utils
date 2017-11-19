using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SA
{
    public abstract class VLargeArray<T>
    {
        public abstract ulong Length
        {
            get;
        }

        public abstract T this[ulong position]
        {
            get;
            set;
        }

        public static VLargeArray<T> AllocateArray(ulong length)
        {
            if ((length >> VLargeArrayLong<T>.ARRAY_SPLIT_2_POW) > 0)
                return new VLargeArrayLong<T>(length);
            else
                return new VLargeArraySmall<T>(length);
        }
    }

    class VLargeArraySmall<T> : VLargeArray<T>
    {
        private ulong mLength;
        private T[] mData;

        public override ulong Length 
        { 
            get { return mLength; } 
        }

        public override T this[ulong position]
        {
            get { return mData[position]; }
            set { mData[position] = value; }
        }

        public VLargeArraySmall(ulong length)
        {
            mLength = length;
            mData = new T[length];
            Array.Clear(mData, 0, (int)length);
        }
    }

    class VLargeArrayLong<T> : VLargeArray<T>
    {
        public const int ARRAY_SPLIT_2_POW = 31;

        private ulong mLength;
        private ulong mArrayMask;

        private T[][] mData;

        public override ulong Length
        {
            get { return mLength; }
        }

        public override T this[ulong position]
        {
            get { return Get(position); }
            set { Set(position, value); }
        }

        public VLargeArrayLong(ulong length)
        {
            mArrayMask = 0;
            for (int i = 0; i < ARRAY_SPLIT_2_POW; i++)
                mArrayMask = (mArrayMask << 1) + 1;

            mLength = length;
            ulong numberOfArrays = (mLength >> ARRAY_SPLIT_2_POW) + 1;
            mData = new T[numberOfArrays][];
            long lengthLeftToAllocate = (long)mLength;
            long allocationBlock = (1 << ARRAY_SPLIT_2_POW);
            for (ulong i = 0; i < numberOfArrays; i++)
            {
                if (allocationBlock > lengthLeftToAllocate)
                    mData[i] = new T[lengthLeftToAllocate];
                else
                    mData[i] = new T[allocationBlock];

                Array.Clear(mData[i], 0, mData[i].Length);

                lengthLeftToAllocate -= mData[i].Length;
            }

            if (lengthLeftToAllocate != 0)
                throw new Exception("lengthLeftToAllocate != 0");
        }

        public void Set(ulong position, T val)
        {
            ulong array = (position >> ARRAY_SPLIT_2_POW);

            mData[array][position & mArrayMask] = val;
        }

        public T Get(ulong position)
        {
            ulong array = (position >> ARRAY_SPLIT_2_POW);

            return mData[array][position & mArrayMask];
        }
    }
}
