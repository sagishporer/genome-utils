using System;
using System.Collections.Generic;
using System.Text;

namespace SA
{
    public class BitIntArray
    {
        const int BitsInArrayEntry = 64;

        ulong[] _array;
        int _arrayLength;
        ulong[] _arrayItemClearMask;

        /// <summary>
        /// The word space for the array (0..k-1)
        /// </summary>
        int _k;

        int _bitsPerItem;
        int _itemsPerArrayEntry;
        ulong _singleItemBitMask;

        public int Length
        {
            get { return _arrayLength; }
        }

        public int this[int idx]
        {
            get
            {
                CheckIndexRange(idx);

                int arrayAbsPosition = idx / _itemsPerArrayEntry;
                int arrayRelPosition = idx % _itemsPerArrayEntry;

                ulong shiftedItem = _array[arrayAbsPosition] >> (arrayRelPosition * _bitsPerItem);
                return (int)(shiftedItem & _singleItemBitMask);
            }
            set
            {
                CheckIndexRange(idx);
                CheckItemValue(value);

                int arrayAbsPosition = idx / _itemsPerArrayEntry;
                int arrayRelPosition = idx % _itemsPerArrayEntry;

                ulong targetShiftedItem = ((ulong)value) << (arrayRelPosition * _bitsPerItem);
                _array[arrayAbsPosition] &= _arrayItemClearMask[arrayRelPosition];
                _array[arrayAbsPosition] ^= targetShiftedItem;
            }
        }

        public BitIntArray(int length, int k)
        {
            k = k + 1;

            _bitsPerItem = (int)Math.Ceiling(Math.Log(k, 2));
            _singleItemBitMask = 0;
            for (int i = 0; i < _bitsPerItem; i++)
                _singleItemBitMask += (uint)Math.Pow(2, i);            

            _itemsPerArrayEntry = BitsInArrayEntry / _bitsPerItem;
            _arrayItemClearMask = new ulong[_itemsPerArrayEntry];
            for (int i = 0; i < _itemsPerArrayEntry; i++)
                _arrayItemClearMask[i] = ~(_singleItemBitMask << (_bitsPerItem * i));

            _k = k;
            _arrayLength = length;

            _array = new ulong[(int)Math.Ceiling((double)_arrayLength / _itemsPerArrayEntry)];
            Array.Clear(_array, 0, _array.Length);
        }

        void CheckIndexRange(int idx)
        {
            if ((idx < 0) || (idx >= _arrayLength))
                throw new ArgumentOutOfRangeException("index == " + idx.ToString());
        }

        void CheckItemValue(int item)
        {
            if ((item < 0) || (item >= _k))
                throw new ArgumentOutOfRangeException("item value == " + item.ToString());
        }

        public int[] ToArray()
        {
            int[] target = new int[_arrayLength];
            for (int i = 0; i < _arrayLength; i++)
                target[i] = this[i];

            return target;
        }

        public static BitIntArray FromArray(int[] source, int k)
        {
            BitIntArray target = new BitIntArray(k, source.Length);
            for (int i = 0; i < source.Length; i++)
                target[i] = source[i];

            return target;
        }
    }
}
