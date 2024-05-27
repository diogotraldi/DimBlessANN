typedef unsigned __int64 my_int;
#define MYINT_SIZE 64

#include <vector>
#include <map>
#include <algorithm>
#include "MovingAverage.h"
#include "Walsh.h"
#include "MyPseudoSort.h"

using namespace std;

namespace DimBless
{
	class SearchSpaceANN
	{
	private:
		const my_int BIT0 = ((my_int)1 << 0); const my_int BIT1 = ((my_int)1 << 1); const my_int BIT2 = ((my_int)1 << 2); const my_int BIT3 = ((my_int)1 << 3); const my_int BIT4 = ((my_int)1 << 4); const my_int BIT5 = ((my_int)1 << 5); const my_int BIT6 = ((my_int)1 << 6); const my_int BIT7 = ((my_int)1 << 7); const my_int BIT8 = ((my_int)1 << 8); const my_int BIT9 = ((my_int)1 << 9); const my_int BIT10 = ((my_int)1 << 10); const my_int BIT11 = ((my_int)1 << 11); const my_int BIT12 = ((my_int)1 << 12); const my_int BIT13 = ((my_int)1 << 13); const my_int BIT14 = ((my_int)1 << 14); const my_int BIT15 = ((my_int)1 << 15); const my_int BIT16 = ((my_int)1 << 16); const my_int BIT17 = ((my_int)1 << 17); const my_int BIT18 = ((my_int)1 << 18); const my_int BIT19 = ((my_int)1 << 19); const my_int BIT20 = ((my_int)1 << 20); const my_int BIT21 = ((my_int)1 << 21); const my_int BIT22 = ((my_int)1 << 22); const my_int BIT23 = ((my_int)1 << 23); const my_int BIT24 = ((my_int)1 << 24); const my_int BIT25 = ((my_int)1 << 25); const my_int BIT26 = ((my_int)1 << 26); const my_int BIT27 = ((my_int)1 << 27); const my_int BIT28 = ((my_int)1 << 28); const my_int BIT29 = ((my_int)1 << 29); const my_int BIT30 = ((my_int)1 << 30); const my_int BIT31 = ((my_int)1 << 31); const my_int BIT32 = ((my_int)1 << 32); const my_int BIT33 = ((my_int)1 << 33); const my_int BIT34 = ((my_int)1 << 34); const my_int BIT35 = ((my_int)1 << 35); const my_int BIT36 = ((my_int)1 << 36); const my_int BIT37 = ((my_int)1 << 37); const my_int BIT38 = ((my_int)1 << 38); const my_int BIT39 = ((my_int)1 << 39); const my_int BIT40 = ((my_int)1 << 40); const my_int BIT41 = ((my_int)1 << 41); const my_int BIT42 = ((my_int)1 << 42); const my_int BIT43 = ((my_int)1 << 43); const my_int BIT44 = ((my_int)1 << 44); const my_int BIT45 = ((my_int)1 << 45); const my_int BIT46 = ((my_int)1 << 46); const my_int BIT47 = ((my_int)1 << 47); const my_int BIT48 = ((my_int)1 << 48); const my_int BIT49 = ((my_int)1 << 49); const my_int BIT50 = ((my_int)1 << 50); const my_int BIT51 = ((my_int)1 << 51); const my_int BIT52 = ((my_int)1 << 52); const my_int BIT53 = ((my_int)1 << 53); const my_int BIT54 = ((my_int)1 << 54); const my_int BIT55 = ((my_int)1 << 55); const my_int BIT56 = ((my_int)1 << 56); const my_int BIT57 = ((my_int)1 << 57); const my_int BIT58 = ((my_int)1 << 58); const my_int BIT59 = ((my_int)1 << 59); const my_int BIT60 = ((my_int)1 << 60); const my_int BIT61 = ((my_int)1 << 61); const my_int BIT62 = ((my_int)1 << 62); const my_int BIT63 = ((my_int)1 << 63);

		const my_int BIT_INT_HALF_H = 0xffffffff00000000;
		const my_int BIT_INT_HALF_L = 0x00000000ffffffff;

		const my_int BIT_SHORT_QUARTER4 = 0xffff000000000000;
		const my_int BIT_SHORT_QUARTER3 = 0x0000ffff00000000;
		const my_int BIT_SHORT_QUARTER2 = 0x00000000ffff0000;
		const my_int BIT_SHORT_QUARTER1 = 0x000000000000ffff;

		const my_int BIT_SHORT_OCT8 = 0xff00000000000000;
		const my_int BIT_SHORT_OCT7 = 0x00ff000000000000;
		const my_int BIT_SHORT_OCT6 = 0x0000ff0000000000;
		const my_int BIT_SHORT_OCT5 = 0x000000ff00000000;
		const my_int BIT_SHORT_OCT4 = 0x00000000ff000000;
		const my_int BIT_SHORT_OCT3 = 0x0000000000ff0000;
		const my_int BIT_SHORT_OCT2 = 0x000000000000ff00;
		const my_int BIT_SHORT_OCT1 = 0x00000000000000ff;

		double limiarStepUp, limiarStepDown;
		int colMin, colStepCheck;
		int stepMonitor;
		long seqCol = 0;

		int _spaceDim = 0; //dimension before transformation. Outside vector dimension.
		int _spaceDimW = 0; //dimension after transformation
		int _spaceLen = 0; //total number of vector of The Space
		int _numOfSpaces = 0;
		int _bitsPerSample = 0; //memory choose for each 

		vector<__int64> m_ID; //keep a values indicating the ID of each array position
		vector<__int64> m_IDFound; //keep the found ids.

		int TABLE_LENGTH = 0;
		my_int *** TABLE_BIT = nullptr;
		my_int * TABLE_RESULT = nullptr;// [156250]; //keep memory safe for 10 million position
		my_int ** TABLE_BIT_TOP = nullptr;// [156250]; //keep memory safe for 10 million position
		my_int ** TABLE_BIT_BOTTOM = nullptr;
		my_int ResetBottomHelper = 0;
		//my_int TABLE_RESULT[15625];//  //keep memory safe for 10 million position

		vector<double> _average;
		vector<double> _stdv;
		vector<double> _stdvQuotient;

		float * _tempF1 = nullptr;
		float * _tempF2 = nullptr;
		float * _tempF3 = nullptr;
		double * _tempD = nullptr;

		int * _tempI1 = nullptr;
		int * _tempI2 = nullptr;

		struct MySortData
		{
			float number;
			int index;
		};
		std::vector<MySortData> mySortItems;

		Walsh _walsh;

		const int totalBuckets = 2000;
		MyPseudoSort * _pseudoSort = nullptr;
		MovingAverage * _movingAverage = nullptr;


		inline float * GetTrasnformW(float * vect)
		{
			float * vectNorm = _tempF1;
			float * vectOrdered = _tempF2;

			//normalize vector
			Norm(vect, vectNorm, _spaceDim);

			float * pW = VectW;
			for (int iSpace = 0; iSpace < _numOfSpaces; iSpace++)
			{
				Reorder(_randomIndex[iSpace], vectNorm, vectOrdered);
				_walsh.fwt(vectOrdered);

				for (int iDim = 1; iDim < _spaceDim; iDim++)
				{
					*pW++ = vectOrdered[iDim];
				}

			}
			return VectW;
		}

		inline void Reorder(vector<int> &index, float * from, float * to)
		{
			size_t len = index.size();

			for (int i = 0; i < len; i++)
				to[i] = from[index[i]];
		}

		inline void Norm(float * from, float * to, int length)
		{
			double sum = 0;
			for (int i = 0; i < length; i++)
			{
				double v = (double)from[i];

				_tempD[i] = v;

				sum += v;
			}
			double average = sum / length;
			double energy = 0;
			for (int i = 0; i < length; i++)
			{
				double v =  _tempD[i] - average;
				energy += (v * v);
				_tempD[i] = v;
			}

			double rms = sqrt(energy); //energy;
			double rmsAdjust = (rms == 0) ? 0 : (1 / rms);
			for (int i = 0; i < length; i++)
			{
				to[i] = (float)(_tempD[i] * rmsAdjust);
			}
		}

		inline int * GetSortedDimension(float * vectW, int topSelection, int arrayLength, int IdExlude)
		{
			int * index = OrderW;
			float * orderValues = _tempF1;

			for (int iDim = 0; iDim < arrayLength; iDim++)
			{
				float value = vectW[iDim];

				orderValues[iDim] = -(float)abs((value - _average[iDim]) * _stdvQuotient[iDim]);
				index[iDim] = (value >= _average[iDim]) ? -(iDim + 1) : iDim + 1;

				if (iDim == IdExlude)
					orderValues[iDim] = 1; //Dimension should be at the last position. 
			}

			PartialSort(orderValues, index, topSelection, arrayLength);

			return index;
		}

		inline bool IsZeroVector(float * vect, int len)
		{
			for (int i = 0; i < len; i++)
			{
				if (vect[i] != 0)
					return false;
			}
			return true;
		}

		inline void PartialSort(float * values, int * index, int topSelection, int arrayLength)
		{
			struct by_number
			{
				bool operator()(MySortData const &left, MySortData const &right)
				{
					return left.number < right.number;
				}
			};

			for (int i = 0; i < arrayLength; i++)
			{
				mySortItems[i].number = values[i];
				mySortItems[i].index = index[i];
			}

			std::partial_sort(mySortItems.begin(), mySortItems.begin() + topSelection, mySortItems.end(), by_number());

			for (int i = 0; i < topSelection; i++)
			{
				index[i] = mySortItems[i].index;
			}
		}

		inline void BitTableLIMITS(bool isTop, my_int * and, int &start, int &stop)
		{
			if (isTop)
			{
				for (int i = start; i < TABLE_LENGTH; i++)
				{
					if (and[i] != 0)
						break;

					start++;
				}

				for (int i = stop - 1; i >= start; i--)
				{
					if (and[i] != 0)
						break;

					stop--;
				}
			}
			else
			{
				for (int i = start; i < TABLE_LENGTH; i++)
				{
					if (~and[i] != 0)
						break;

					start++;
				}

				for (int i = stop - 1; i >= start; i--)
				{
					if (~and[i] != 0)
						break;

					stop--;
				}
			}
		}

		inline int GetMemoryIndex(float value, int iDimSignal, double limiar)
		{
			int iDim = abs(iDimSignal) - 1;
			float * mem = &(_roundMem[iDim][0]);
			if (iDimSignal > 0)
			{
				double limiarTop = value + limiar;

				if (mem[_bitsPerSample - 1] < limiarTop)
					return -1;

				for (int m = _bitsPerSample - 2; m >= 0; m--)
				{
					if (mem[m] < limiarTop)
						return m + 1;
				}
				return 0;
			}
			else
			{
				double limiarBottom = value - limiar;

				if (limiarBottom < mem[0])
					return -1;

				for (int m = 1; m < _bitsPerSample; m++)
				{
					if (mem[m] > limiarBottom)
						return m - 1;
				}
				return _bitsPerSample - 1;
			}
		}

		inline int GetSpottedCount(int start, int stop, int bruteTarget)
		{
			if (start == stop)
				return 0;

			int numOfFound = 0;

			int detailLevel = (4 * bruteTarget / (stop - start)) + 1;
			if (detailLevel > 3)
				detailLevel = 3;
			
			switch (detailLevel)
			{
			case 1:
				for (int j = start; j < stop; j += stepMonitor)
				{
					if (TABLE_RESULT[j] != 0)
						numOfFound++;
				}
				break;

			case 2:
				for (int j = start; j < stop; j += stepMonitor)
				{
					my_int y = TABLE_RESULT[j];
					if (y == 0)
						continue;

					if ((BIT_INT_HALF_L & y) > 0)
					{
						numOfFound++;
					}
					if ((y & BIT_INT_HALF_H) > 0)
					{
						numOfFound++;
					}
				}
				break;

			case 3:
				for (int j = start; j < stop; j += stepMonitor)
				{
					my_int y = TABLE_RESULT[j];
					if (y == 0)
						continue;

					if ((BIT_INT_HALF_L & y) > 0)
					{
						if ((BIT_SHORT_QUARTER1 & y) > 0)
							numOfFound++;

						if ((BIT_SHORT_QUARTER2 & y) > 0)
							numOfFound++;

					}
					if ((y & BIT_INT_HALF_H) > 0)
					{
						if ((BIT_SHORT_QUARTER3 & y) > 0)
							numOfFound++;

						if ((BIT_SHORT_QUARTER4 & y) > 0)
							numOfFound++;
					}
				}
				break;

			default:
				for (int j = start; j < stop; j += stepMonitor)
				{
					if (TABLE_RESULT[j] != 0)
						numOfFound++;
				}
				break;
			}

			return numOfFound*stepMonitor;
		}

		inline void SetBitTopN(my_int ** __restrict fromN, my_int * __restrict to, int start, int stop, int count)
		{
			to += start;
			my_int * __restrict from1 = fromN[0] + start;
			my_int * __restrict from2 = fromN[1] + start;
			my_int * __restrict from3 = fromN[2] + start;
			my_int * __restrict from4 = fromN[3] + start;
			my_int * __restrict from5 = fromN[4] + start;

			int len = stop - start;
			switch (count)
			{
			case 5:
				for (int i = len; i--; )
				{
					*to++ &= (*from1++) & (*from2++) & (*from3++) & (*from4++) & (*from5++);
				}
				break;

			case 4:
				for (int i = len; i--; )
				{
					*to++ &= (*from1++) & (*from2++) & (*from3++) & (*from4++);
				}

				break;

			case 3:
				for (int i = len; i--; )
				{
					*to++ &= (*from1++) & (*from2++) & (*from3++);
				}
				break;

			case 2:
				for (int i = len; i--; )
				{
					*to++ &= (*from1++) & (*from2++);
				}
				break;

			case 1:
				for (int i = len; i--; )
				{
					*to++ &= *from1++;
				}
				break;

			default:
				break;
			}

		}

		inline void SetBitBottomN(my_int ** __restrict fromN, my_int * __restrict to, int start, int stop, int count)
		{
			to += start;
			my_int * __restrict from1 = fromN[0] + start;
			my_int * __restrict from2 = fromN[1] + start;
			my_int * __restrict from3 = fromN[2] + start;
			my_int * __restrict from4 = fromN[3] + start;
			my_int * __restrict from5 = fromN[4] + start;

			int len = stop - start;
			switch (count)
			{
			case 5:
				for (int i = len; i--; )
				{
					*to++ &= (~*from1++) & (~*from2++) & (~*from3++) & (~*from4++) & (~*from5++);
				}
				break;

			case 4:
				for (int i = len; i--; )
				{
					*to++ &= (~*from1++) & (~*from2++) & (~*from3++) & (~*from4++);
				}
				break;

			case 3:
				for (int i = len; i--; )
				{
					*to++ &= (~*from1++) & (~*from2++) & (~*from3++);
				}
				break;

			case 2:
				for (int i = len; i--; )
				{
					*to++ &= (~*from1++) & (~*from2++);
				}
				break;

			case 1:
				for (int i = len; i--; )
				{
					*to++ &= ~*from1++;
				}
				break;

			default:
				break;
			}

		}

		inline void ResetBitTop(my_int * __restrict from, my_int * __restrict to, int start, int stop)
		{
			memcpy(to + start, from + start, (stop - start) * sizeof(my_int));
		}

		inline void ResetBitBottom(my_int * __restrict from, my_int * __restrict to, int start, int stop)
		{
			int len = stop - start;

			to += start;
			from += start;

			for (int i = len; i--; )
			{
				*to++ = ~(*from++);
			}

			to--;
			*to &= ResetBottomHelper;

			//int len = stop - start;

			//to += start;
			//from += start;

			//for (int i = len / 4; i--; )
			//{
			//	*to++ = ~(*from++);
			//	*to++ = ~(*from++);
			//	*to++ = ~(*from++);
			//	*to++ = ~(*from++);
			//}

			//for (int i = len % 4; i--; )
			//{
			//	*to++ = ~(*from++);
			//}
		}

		int LoadFoundBits(int start, int stop, __int64 * foundID)
		{
			int numOfFound = 0;

			__int64 * idFound = &(m_IDFound[0]);
			int baseId = start * MYINT_SIZE;
			for (int j = start; j < stop; j++, baseId += MYINT_SIZE)
			{
				my_int y = TABLE_RESULT[j];
				if (y == 0)
					continue;

				if ((BIT_INT_HALF_L & y) > 0)
				{
					if ((BIT_SHORT_QUARTER1 & y) > 0)
					{
						if ((BIT_SHORT_OCT1 & y) > 0)
						{
							if ((BIT0 & y) > 0) { *idFound++ = baseId + 0; numOfFound++; }
							if ((BIT1 & y) > 0) { *idFound++ = baseId + 1; numOfFound++; }
							if ((BIT2 & y) > 0) { *idFound++ = baseId + 2; numOfFound++; }
							if ((BIT3 & y) > 0) { *idFound++ = baseId + 3; numOfFound++; }
							if ((BIT4 & y) > 0) { *idFound++ = baseId + 4; numOfFound++; }
							if ((BIT5 & y) > 0) { *idFound++ = baseId + 5; numOfFound++; }
							if ((BIT6 & y) > 0) { *idFound++ = baseId + 6; numOfFound++; }
							if ((BIT7 & y) > 0) { *idFound++ = baseId + 7; numOfFound++; }
						}

						if ((BIT_SHORT_OCT2 & y) > 0)
						{
							if ((BIT8 & y) > 0) { *idFound++ = baseId + 8; numOfFound++; }
							if ((BIT9 & y) > 0) { *idFound++ = baseId + 9; numOfFound++; }
							if ((BIT10 & y) > 0) { *idFound++ = baseId + 10; numOfFound++; }
							if ((BIT11 & y) > 0) { *idFound++ = baseId + 11; numOfFound++; }
							if ((BIT12 & y) > 0) { *idFound++ = baseId + 12; numOfFound++; }
							if ((BIT13 & y) > 0) { *idFound++ = baseId + 13; numOfFound++; }
							if ((BIT14 & y) > 0) { *idFound++ = baseId + 14; numOfFound++; }
							if ((BIT15 & y) > 0) { *idFound++ = baseId + 15; numOfFound++; }
						}
					}

					if ((BIT_SHORT_QUARTER2 & y) > 0)
					{
						if ((BIT_SHORT_OCT3 & y) > 0)
						{
							if ((BIT16 & y) > 0) { *idFound++ = baseId + 16; numOfFound++; }
							if ((BIT17 & y) > 0) { *idFound++ = baseId + 17; numOfFound++; }
							if ((BIT18 & y) > 0) { *idFound++ = baseId + 18; numOfFound++; }
							if ((BIT19 & y) > 0) { *idFound++ = baseId + 19; numOfFound++; }
							if ((BIT20 & y) > 0) { *idFound++ = baseId + 20; numOfFound++; }
							if ((BIT21 & y) > 0) { *idFound++ = baseId + 21; numOfFound++; }
							if ((BIT22 & y) > 0) { *idFound++ = baseId + 22; numOfFound++; }
							if ((BIT23 & y) > 0) { *idFound++ = baseId + 23; numOfFound++; }
						}
						if ((BIT_SHORT_OCT4 & y) > 0)
						{
							if ((BIT24 & y) > 0) { *idFound++ = baseId + 24; numOfFound++; }
							if ((BIT25 & y) > 0) { *idFound++ = baseId + 25; numOfFound++; }
							if ((BIT26 & y) > 0) { *idFound++ = baseId + 26; numOfFound++; }
							if ((BIT27 & y) > 0) { *idFound++ = baseId + 27; numOfFound++; }
							if ((BIT28 & y) > 0) { *idFound++ = baseId + 28; numOfFound++; }
							if ((BIT29 & y) > 0) { *idFound++ = baseId + 29; numOfFound++; }
							if ((BIT30 & y) > 0) { *idFound++ = baseId + 30; numOfFound++; }
							if ((BIT31 & y) > 0) { *idFound++ = baseId + 31; numOfFound++; }
						}
					}
				}

				if ((y & BIT_INT_HALF_H) > 0)
				{
					if ((BIT_SHORT_QUARTER3 & y) > 0)
					{
						if ((BIT_SHORT_OCT5 & y) > 0)
						{
							if ((BIT32 & y) > 0) { *idFound++ = baseId + 32; numOfFound++; }
							if ((BIT33 & y) > 0) { *idFound++ = baseId + 33; numOfFound++; }
							if ((BIT34 & y) > 0) { *idFound++ = baseId + 34; numOfFound++; }
							if ((BIT35 & y) > 0) { *idFound++ = baseId + 35; numOfFound++; }
							if ((BIT36 & y) > 0) { *idFound++ = baseId + 36; numOfFound++; }
							if ((BIT37 & y) > 0) { *idFound++ = baseId + 37; numOfFound++; }
							if ((BIT38 & y) > 0) { *idFound++ = baseId + 38; numOfFound++; }
							if ((BIT39 & y) > 0) { *idFound++ = baseId + 39; numOfFound++; }
						}

						if ((BIT_SHORT_OCT6 & y) > 0)
						{
							if ((BIT40 & y) > 0) { *idFound++ = baseId + 40; numOfFound++; }
							if ((BIT41 & y) > 0) { *idFound++ = baseId + 41; numOfFound++; }
							if ((BIT42 & y) > 0) { *idFound++ = baseId + 42; numOfFound++; }
							if ((BIT43 & y) > 0) { *idFound++ = baseId + 43; numOfFound++; }
							if ((BIT44 & y) > 0) { *idFound++ = baseId + 44; numOfFound++; }
							if ((BIT45 & y) > 0) { *idFound++ = baseId + 45; numOfFound++; }
							if ((BIT46 & y) > 0) { *idFound++ = baseId + 46; numOfFound++; }
							if ((BIT47 & y) > 0) { *idFound++ = baseId + 47; numOfFound++; }
						}
					}

					if ((BIT_SHORT_QUARTER4 & y) > 0)
					{
						if ((BIT_SHORT_OCT7 & y) > 0)
						{
							if ((BIT48 & y) > 0) { *idFound++ = baseId + 48; numOfFound++; }
							if ((BIT49 & y) > 0) { *idFound++ = baseId + 49; numOfFound++; }
							if ((BIT50 & y) > 0) { *idFound++ = baseId + 50; numOfFound++; }
							if ((BIT51 & y) > 0) { *idFound++ = baseId + 51; numOfFound++; }
							if ((BIT52 & y) > 0) { *idFound++ = baseId + 52; numOfFound++; }
							if ((BIT53 & y) > 0) { *idFound++ = baseId + 53; numOfFound++; }
							if ((BIT54 & y) > 0) { *idFound++ = baseId + 54; numOfFound++; }
							if ((BIT55 & y) > 0) { *idFound++ = baseId + 55; numOfFound++; }
						}

						if ((BIT_SHORT_OCT8 & y) > 0)
						{
							if ((BIT56 & y) > 0) { *idFound++ = baseId + 56; numOfFound++; }
							if ((BIT57 & y) > 0) { *idFound++ = baseId + 57; numOfFound++; }
							if ((BIT58 & y) > 0) { *idFound++ = baseId + 58; numOfFound++; }
							if ((BIT59 & y) > 0) { *idFound++ = baseId + 59; numOfFound++; }
							if ((BIT60 & y) > 0) { *idFound++ = baseId + 60; numOfFound++; }
							if ((BIT61 & y) > 0) { *idFound++ = baseId + 61; numOfFound++; }
							if ((BIT62 & y) > 0) { *idFound++ = baseId + 62; numOfFound++; }
							if ((BIT63 & y) > 0) { *idFound++ = baseId + 63; numOfFound++; }
						}
					}
				}
			}

			__int64 * result = (foundID != nullptr) ? foundID : &(m_IDFound[0]);

			for (size_t i = 0; i < numOfFound; i++)
			{
				result[i] = m_ID[m_IDFound[i]];
			}

			return numOfFound;
		}

		void GetSTDV(vector<float> &vect, int length, double &average, double &stdv)
		{
			double sum = 0;
			for (size_t i = 0; i < length; i++)
			{
				sum += vect[i];
			}
			average = sum / length;

			double var = 0;
			for (size_t i = 0; i < length; i++)
			{
				double v = vect[i] - average;
				var += v * v;
			}

			stdv = sqrt(var / length);
		}

		int IndexOfMax(vector<double> &vect, int length)
		{
			double maxV = vect[0];
			int index = 0;
			for (int i = 1; i < length; i++)
			{
				double absolute = abs(vect[i]);
				if (maxV < absolute)
				{
					maxV = absolute;
					index = i;
				}
			}

			return index;
		}

		void SetMemoryTable(vector<float> &values, vector<int> &position, vector<float> &memory, my_int ** BIT_TABLE)
		{
			vector<int> sortedPosition(position);
			//std::sort(sortedPosition.begin(), sortedPosition.end(), [&values](float i1, size_t i2) {return values[i1] < values[i2]; });

			_pseudoSort->PeseudoSort(values, position, sortedPosition, memory);

			double step = (double)_spaceLen / (_bitsPerSample + 1);
			double stepSum = 0;

			vector<my_int> BIT_TABLE_DRAFT(TABLE_LENGTH, 0);

			int iStep = 0;
			for (int iMem = 0; iMem < _bitsPerSample; iMem++)
			{
				stepSum += step;
				for (; iStep < (int)(stepSum); iStep++)
				{
					int vectPos = sortedPosition[iStep];

					int bit = (vectPos % MYINT_SIZE);
					int pos = vectPos / MYINT_SIZE;

					BIT_TABLE_DRAFT[pos] |= ((my_int)1 << bit); //set bit
				}

				//copy values from BIT_TABLE_DRAFT to BIT_TABLE
				std::copy(BIT_TABLE_DRAFT.begin(), BIT_TABLE_DRAFT.end(), BIT_TABLE[iMem]);
			}
		}

		void SetKnownMemoryTable(vector<float> &values, vector<int> &position, vector<float> &memory, my_int ** BIT_TABLE)
		{
			vector<my_int> BIT_TABLE_DRAFT(TABLE_LENGTH, 0);

			for (int iMem = 0; iMem < _bitsPerSample; iMem++)
			{
				double memTop = memory[iMem];
				double memBottom = iMem == 0 ? -DBL_MAX : memory[iMem-1];
				for (int iStep = 0; iStep < _spaceLen; iStep++)
				{
					double value = values[iStep];

					if (value < memBottom || value >= memTop)
						continue;

					int vectPos = position[iStep];

					int bit = (vectPos % MYINT_SIZE);
					int pos = vectPos / MYINT_SIZE;

					BIT_TABLE_DRAFT[pos] |= ((my_int)1 << bit); //set bit
				}

				//copy values from BIT_TABLE_DRAFT to BIT_TABLE
				std::copy(BIT_TABLE_DRAFT.begin(), BIT_TABLE_DRAFT.end(), BIT_TABLE[iMem]);
			}
		}

		void InstanceSpaceMemory(int spaceLen, int bitsPerSample)
		{
			if (_spaceLen != spaceLen)
			{
				m_ID.resize(spaceLen);
				m_IDFound.resize(spaceLen);
			}

			ClearSpaceMemory();

			_spaceLen = spaceLen;

			int table_length = ((spaceLen % MYINT_SIZE == 0)) ? (spaceLen / MYINT_SIZE) : (spaceLen / MYINT_SIZE) + 1;
			_bitsPerSample = bitsPerSample;

			TABLE_LENGTH = table_length;
			TABLE_RESULT = new my_int[TABLE_LENGTH];
			TABLE_BIT_TOP = new my_int*[_spaceDimW];
			TABLE_BIT_BOTTOM = new my_int*[_spaceDimW];

			TABLE_BIT = new my_int**[_spaceDimW];
			for (int iDim = 0; iDim < _spaceDimW; iDim++)
			{
				TABLE_BIT[iDim] = new my_int*[bitsPerSample];
				for (int iMem = 0; iMem < bitsPerSample; iMem++)
					TABLE_BIT[iDim][iMem] = new my_int[TABLE_LENGTH];
			}

			ResetBottomHelper = 0;
			int uselessBit = _spaceLen % MYINT_SIZE;
			for (size_t bit = 0; bit < uselessBit; bit++)
			{
				ResetBottomHelper |= ((my_int)1 << bit);
			}

			keepColElim.reserve(1000);
			keepLimiarElim.reserve(1000);
			keepCountElim.reserve(1000);
		}

		void ClearGeneralMemory()
		{
			SafeDeleteArray(_tempF1);
			SafeDeleteArray(_tempF2);
			SafeDeleteArray(_tempF3);
			SafeDeleteArray(_tempI1);
			SafeDeleteArray(_tempI2);
			SafeDeleteArray(_tempD);

			SafeDeleteArray(VectW);
			SafeDeleteArray(OrderW);
			SafeDeleteArray(ReptitionControl);

			delete _pseudoSort;
			_pseudoSort = nullptr;

			delete _movingAverage;
			_movingAverage = nullptr;
		}

		void InstanceGeneralMemory(int spaceDim, int spaceDimW, int tempMemSize)
		{
			ClearGeneralMemory();

			_tempF1 = new float[tempMemSize];
			_tempF2 = new float[tempMemSize];
			_tempF3 = new float[tempMemSize];
			_tempI1 = new int[tempMemSize];
			_tempI2 = new int[tempMemSize];
			_tempD = new double[tempMemSize];

			VectW = new float[spaceDimW];
			OrderW = new int[spaceDimW];
			ReptitionControl = new int[spaceDimW];

			_pseudoSort = new MyPseudoSort(totalBuckets);
			_movingAverage = new MovingAverage(20);
		}

		template< class T > void SafeDeleteArray(T*& pVal)
		{
			if (pVal != NULL)
			{
				delete[] pVal;
				pVal = NULL;
			}
		}

	public:
		int iTop = 0;
		int iBot = 0;
		int _mainCol;
		int colNum;
		double limiar;
		float * VectW = nullptr;
		int * OrderW = nullptr;
		int * ReptitionControl = nullptr;

		vector<double> keepColElim;
		vector<double> keepLimiarElim;
		vector<double> keepCountElim;

		vector<vector<int>> _randomIndex; //randomize indedices before transformation
		vector<vector<float>> _roundMem;

		SearchSpaceANN(int spaceDim, int numOfSpaces): SearchSpaceANN(spaceDim, numOfSpaces, vector<vector<int>>())
		{
		}

		SearchSpaceANN(int spaceDim, int numOfSpaces, vector<vector<int>> &randomIndex)
		{
			_spaceDim = spaceDim;
			_numOfSpaces = numOfSpaces;
			_spaceDimW = (spaceDim - 1) * numOfSpaces;

			_walsh.Init(spaceDim);

			_average.resize(_spaceDimW);
			_stdv.resize(_spaceDimW);
			_stdvQuotient.resize(_spaceDimW);

			mySortItems.resize(_spaceDimW);

			if (randomIndex.size() > 0)
			{
				_randomIndex = randomIndex;
			}
			else
			{
				_randomIndex.resize(numOfSpaces, vector<int>(spaceDim));
				for (int iSpace = 0; iSpace < numOfSpaces; iSpace++)
				{
					for (int iDim = 0; iDim < spaceDim; iDim++)
					{
						_randomIndex[iSpace][iDim] = iDim;
					}

					random_shuffle(_randomIndex[iSpace].begin(), _randomIndex[iSpace].end());
				}
			}

			int tempMemSize = max(_spaceDimW, _spaceDim);

			InstanceGeneralMemory(spaceDim, _spaceDimW, tempMemSize);
		}

		~SearchSpaceANN()
		{
			_walsh.Close();
			Clear();
		}

		int Size()
		{
			return _spaceLen;
		}

		void IndexSpace(map<__int64, float*> &vectSpace, int bitsPerSample)
		{
			vector<vector<float>> emptyRoundMem;
			IndexSpace(vectSpace, bitsPerSample, emptyRoundMem, -1);
		}

		void IndexSpace(map<__int64, float*> &vectSpace, int bitsPerSample, vector<vector<float>> &roundMem, int mainCol)
		{
			InstanceSpaceMemory((int)vectSpace.size(), bitsPerSample);
			if (_spaceLen == 0)
				return;

			bool knownMemory = (roundMem.size() > 0);
			if (knownMemory)
				_roundMem = roundMem;
			else
				_roundMem.resize(_spaceDimW, vector<float>(bitsPerSample));

			vector<vector<float>> values(_spaceDim, vector<float>(_spaceLen));
			vector<int> reorderedVect(_spaceLen);

			for (int iSpace = 0; iSpace < _numOfSpaces; iSpace++)
			{
				int pos = 0;
				float * vectW = _tempF1; //getting an already allocated memory. It's temporary variable
				for (map<__int64, float*>::iterator it = vectSpace.begin(); it != vectSpace.end(); it++)
				{
					Reorder(_randomIndex[iSpace], it->second, vectW);
					Norm(vectW, vectW, _spaceDim);

					_walsh.fwt(vectW);

					for (int v = 1; v < _spaceDim; v++)
					{
						values[v][pos] = vectW[v];
					}
					pos++;
				}

				for (int i = 1, w = iSpace*(_spaceDim - 1); i < _spaceDim; i++, w++)
				{
					double average, stdv = 0;
					GetSTDV(values[i], _spaceLen, average, stdv);

					_average[w] = average;
					_stdv[w] = stdv;
					_stdvQuotient[w] = 1 / stdv;
				}

				if (iSpace == 0)
				{
					_mainCol = mainCol >= 0 ? mainCol : IndexOfMax(_stdv, (_spaceDim - 1));

					//load all space vectors IDs into a vector 
					vector<__int64> spaceId;
					for (map<__int64, float*>::iterator it = vectSpace.begin(); it != vectSpace.end(); it++)
						spaceId.push_back(it->first);

					//finding out new indexes
					int * sortedIndex = new int[_spaceLen];
					for (int i = 0; i < _spaceLen; i++)
						sortedIndex[i] = i;

					float * v = &values[_mainCol + 1][0];
					sort(sortedIndex, sortedIndex + _spaceLen, [&](int a, int b) { return v[a] < v[b]; });

					for (int iOrder = 0; iOrder < _spaceLen; iOrder++)
					{
						reorderedVect[sortedIndex[iOrder]] = iOrder; //write the order (new position) on the right position
						m_ID[iOrder] = spaceId[sortedIndex[iOrder]];
					}
				}


				for (int i = 1, w = iSpace * (_spaceDim - 1); i < _spaceDim; i++, w++)
				{
					if (knownMemory)
						SetKnownMemoryTable(values[i], reorderedVect, _roundMem[w], TABLE_BIT[w]);
					else
						SetMemoryTable(values[i], reorderedVect, _roundMem[w], TABLE_BIT[w]);
				}
			}
		}

		void SetDynamicParameters(int colNum, int colMin, int colStepCheck, float limiarStepUp, float limiarStepDown, int stepMonitor)
		{
			this->colNum = colNum;
			this->colMin = colMin;
			this->colStepCheck = colStepCheck;
			this->limiarStepUp = (1.0 + limiarStepUp);
			this->limiarStepDown = (1.0 - limiarStepDown);
			this->stepMonitor = stepMonitor;
		}

		int SearchDynamic(float * vect, int bruteTarget)
		{
			return SearchDynamic(vect, false, nullptr, bruteTarget);
		}

		int SearchDynamic(float * vect, bool debug, __int64 * FoundID, int bruteTarget)
		{
			if (IsZeroVector(vect, _spaceLen))
				return 0;

			float * vectW = GetTrasnformW(vect);
			int * order = GetSortedDimension(vectW, colNum, _spaceDimW, _mainCol);

			if (debug)
			{
				keepColElim.clear();
				keepLimiarElim.clear();
				keepCountElim.clear();
			}

			limiar = _movingAverage->size() == 0 ? _stdv[_mainCol] : _movingAverage->GetAverage();
			long totalCol = 0;
			int countSpotted = _spaceLen;
			int start = 0; //this is the portion sampled to estimate limiar
			int stop = TABLE_LENGTH;//this is the portion sampled to estimate limiar
			int iCol = 0;
			int countSaturated = 0;
			bool jobDone = false;
			bool resetTable = true; //variable indicates to assign table rather than 
			bool limiarChanged = true;
			double maxLimiar = limiar;

			while (!jobDone)
			{
				if (resetTable)
				{
					seqCol = 0;
					countSpotted = _spaceLen;

					for (int i = 0; i < colNum; i++)
						ReptitionControl[i] = INT_MIN;

					countSaturated = 0;

					start = 0; //this is the portion sampled to estimate limiar
					stop = TABLE_LENGTH;//this is the portion sampled to estimate limiar

					iCol = 0;
				}

				if (limiarChanged || resetTable)
				{
					int iMemMain = GetMemoryIndex(vectW[_mainCol], -(_mainCol + 1), limiar);
					if (iMemMain >= 0)
						BitTableLIMITS(false, TABLE_BIT[_mainCol][iMemMain], start, stop);

					iMemMain = GetMemoryIndex(vectW[_mainCol], (_mainCol + 1), limiar);
					if (iMemMain >= 0)
						BitTableLIMITS(true, TABLE_BIT[_mainCol][iMemMain], start, stop);
				}

				int iCheck = 0;
				int iTop = 0;
				int iBot = 0;
				while (iCheck < colStepCheck)
				{
					if (iCol >= colNum)
					{
						iCol = 0;
						limiar *= 0.99; //always decrease limiar after run over all columns
					}

					int iDim = abs(order[iCol]) - 1;
					int iMem = GetMemoryIndex(vectW[iDim], order[iCol], limiar);

					if (iMem != ReptitionControl[iCol])
					{
						bool isTop = (order[iCol] > 0);
						if (iMem >= 0 && iMem < _bitsPerSample)
						{
							iCheck++;
							if (resetTable)
							{
								if (isTop)
									ResetBitTop(TABLE_BIT[iDim][iMem], TABLE_RESULT, start, stop);
								else
									ResetBitBottom(TABLE_BIT[iDim][iMem], TABLE_RESULT, start, stop);

								resetTable = false;
							}
							else
							{
								if (isTop)
									TABLE_BIT_TOP[iTop++] = TABLE_BIT[iDim][iMem];
								else
									TABLE_BIT_BOTTOM[iBot++] = TABLE_BIT[iDim][iMem];
							}
						}

						ReptitionControl[iCol] = iMem;
						if (isTop)
						{
							if (iMem == 0)
								countSaturated++;
						}
						else
						{
							if (iMem == _bitsPerSample - 1)
								countSaturated++;
						}

						if (countSaturated >= colNum)
						{
							jobDone = true;
							break;
						}
					}

					iCol++;
				}

				for (int i = 0; i < iTop;)
				{
					int count = iTop - i;

					int div;
					for (div = 1; count / div >= 5; div++);
					count /= div;

					SetBitTopN(TABLE_BIT_TOP + i, TABLE_RESULT, start, stop, count);
					i += count;
				}

				for (int i = 0; i < iBot;)
				{
					int count = iBot - i;

					int div;
					for (div = 1; count / div >= 5; div++);
					count /= div;

					SetBitBottomN(TABLE_BIT_BOTTOM + i, TABLE_RESULT, start, stop, count);
					i += count;
				}

				totalCol += iCheck;
				seqCol += iCheck;

				int lastSpotted = countSpotted;
				countSpotted = GetSpottedCount(start, stop, bruteTarget);

				if (debug)
				{
					keepColElim.push_back(totalCol);
					keepLimiarElim.push_back(limiar);
					keepCountElim.push_back(countSpotted);
				}

				if (seqCol >= colMin && countSpotted <= bruteTarget)
				{
					jobDone = true;
					continue;
				}

				if (seqCol < colMin && countSpotted < bruteTarget)
				{
					if (limiar*limiarStepUp <= maxLimiar)
					{
						//It may 
						jobDone = true;
						continue;
					}

					limiar *= limiarStepUp;
					maxLimiar = limiar;
					limiarChanged = true;

					resetTable = true;
					continue;
				}

				long colLasting = seqCol >= colNum ? 1 : (colNum - seqCol) / colStepCheck + 1;
				double layerCount = ((double)countSpotted - bruteTarget * 0.9) / (lastSpotted - countSpotted) / colLasting;
				if (layerCount > 1)
				{
					limiar *= limiarStepDown;
					limiarChanged = true;
					continue;
				}
			}

			_movingAverage->AddValue(limiar);

			int numOfFound = LoadFoundBits(start, stop, FoundID);

			return numOfFound;
		}

		int SearchStatic(float * vect, double limiar, int colNum)
		{
			if (IsZeroVector(vect, _spaceLen))
				return 0;

			float * vectW = GetTrasnformW(vect);
			int * orderW = GetSortedDimension(vectW, colNum, _spaceDimW, _mainCol);

			//find first available memory e reset table
			for (int iCol = 0; iCol < colNum; iCol++)
			{
				int iDim = abs(orderW[iCol]) - 1;
				int iMem = GetMemoryIndex(vectW[iDim], orderW[iCol], limiar);

				ReptitionControl[iCol] = iMem;
			}

			return SearchStaticW(vectW, orderW, ReptitionControl, limiar, colNum, nullptr);
		}

		int SearchStaticW(float * vectW, int *orderW, int * memoryIndexes, double limiar, int colNum, __int64 * foundID)
		{
			if (IsZeroVector(vectW, colNum))
				return 0;

			int start = 0; //this is the portion sampled to estimate limiar
			int stop = TABLE_LENGTH;//this is the portion sampled to estimate limiar

			int iMemMain = GetMemoryIndex(vectW[_mainCol], -(_mainCol + 1), limiar);
			if (iMemMain >= 0)
				BitTableLIMITS(false, TABLE_BIT[_mainCol][iMemMain], start, stop);

			iMemMain = GetMemoryIndex(vectW[_mainCol], (_mainCol + 1), limiar);
			if (iMemMain >= 0)
				BitTableLIMITS(true, TABLE_BIT[_mainCol][iMemMain], start, stop);

			int iCol = 0;
			int iDim, iMem;
			//find first available memory e reset table
			for (; iCol < colNum; iCol++)
			{
				iDim = abs(orderW[iCol]) - 1;
				iMem = memoryIndexes[iCol];

				if (iMem >= 0 && iMem < _bitsPerSample)
				{
					if (orderW[iCol] > 0)
						ResetBitTop(TABLE_BIT[iDim][iMem], TABLE_RESULT, start, stop);
					else
						ResetBitBottom(TABLE_BIT[iDim][iMem], TABLE_RESULT, start, stop);

					iCol++;
					break;
				}
			}

			//Assing TABLE RESULT. (Reset TABLE Result) 
			iTop = 0;
			iBot = 0;
			for (; iCol < colNum; iCol++)
			{
				iMem = memoryIndexes[iCol];
				iDim = abs(orderW[iCol]) - 1;

				if (iMem >= 0 && iMem < _bitsPerSample)
				{
					if (orderW[iCol] > 0)
						TABLE_BIT_TOP[iTop++] = TABLE_BIT[iDim][iMem];
					else
						TABLE_BIT_BOTTOM[iBot++] = TABLE_BIT[iDim][iMem];
				}
			}

			for (int i = 0; i < iTop;)
			{
				int count = iTop - i;

				int div;
				for (div = 1; count / div >= 5; div++);
				count /= div;

				SetBitTopN(TABLE_BIT_TOP + i, TABLE_RESULT, start, stop, count);
				i += count;
			}

			for (int i = 0; i < iBot;)
			{
				int count = iBot - i;

				int div;
				for (div = 1; count / div >= 5; div++);
				count /= div;

				SetBitBottomN(TABLE_BIT_BOTTOM + i, TABLE_RESULT, start, stop, count);
				i += count;
			}

			int numOfFound = LoadFoundBits(start, stop, foundID);

			return numOfFound;
		}

		__int64 * GetFoundID()
		{
			return &(m_IDFound[0]);
		}

		vector<__int64> GetSortedSpaceID()
		{
			vector<__int64> vectID(_spaceLen);

			for (int i = 0; i < _spaceLen; i++)
			{
				vectID[i] = m_ID[i];
			}

			return vectID;
		}

		void Clear()
		{
			ClearSpaceMemory();
			ClearGeneralMemory();
		}

		void ClearSpaceMemory()
		{
			_spaceLen = 0;
			TABLE_LENGTH = 0;
			ResetBottomHelper = 0;

			SafeDeleteArray(TABLE_RESULT);
			SafeDeleteArray(TABLE_BIT_TOP);
			SafeDeleteArray(TABLE_BIT_BOTTOM);

			if (TABLE_BIT != nullptr)
			{
				for (int iDim = 0; iDim < _spaceDimW; iDim++)
					for (int iMem = 0; iMem < _bitsPerSample; iMem++)
						SafeDeleteArray(TABLE_BIT[iDim][iMem]);

				for (int iDim = 0; iDim < _spaceDimW; iDim++)
					SafeDeleteArray(TABLE_BIT[iDim]);

				SafeDeleteArray(TABLE_BIT);
			}
		}

	};
}