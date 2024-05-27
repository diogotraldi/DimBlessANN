#include <vector>

using namespace std;

namespace DimBless
{
	class MyPseudoSort
	{	
	private:
		//these variables help Pseudo sort algorithm
		int totalBuckets;
		vector<int> histogram;
		vector<vector<int>> matrixPosition;
		vector<vector<double>> matrixValues;


		struct data
		{
			double number;
			int index;
		};
		std::vector<data> items;


	public:
		MyPseudoSort(int total_buckets)
		{
			totalBuckets = total_buckets;
			histogram.resize(totalBuckets, 0);

			matrixPosition.resize(totalBuckets, vector<int>());
			matrixValues.resize(totalBuckets, vector<double>());
		}

		void PeseudoSort(vector<float> &values, vector<int> &position, vector<int> &sortedPosition, vector<float> &divisionValues)
		{
			int numOfDivision = (int)divisionValues.size();
			int length = (int)sortedPosition.size();

			double max = values[0], min = values[0];
			for (int i = 1; i < length; i++)
			{
				float value = values[i];

				if (value > max)
					max = value;

				if (value < min)
					min = value;
			}

			std::fill(histogram.begin(), histogram.end(), 0);

			double bucketSize = (max - min) / totalBuckets;
			for (int i = 0; i < length; i++)
			{
				double value = values[i];

				int bucketIndex = 0;
				if (bucketSize > 0.0)
				{
					bucketIndex = (int)((values[i] - min) / bucketSize);
					if (bucketIndex == totalBuckets)
						bucketIndex--;
				}
				int len = histogram[bucketIndex];

				if (matrixValues[bucketIndex].size() <= len)
				{
					matrixValues[bucketIndex].resize(2 * (len + 1));
					matrixPosition[bucketIndex].resize(2 * (len + 1));
				}

				histogram[bucketIndex] = len + 1;
				matrixValues[bucketIndex][len] = value;
				matrixPosition[bucketIndex][len] = position[i];
			}


			double step = (double)length / (numOfDivision + 1);
			double nextStep = step;
			int iDiv = 0;
			int count = 0;
			for (int i = 0; i < totalBuckets; i++)
			{
				int len = histogram[i];
				if (iDiv < numOfDivision && count + len > nextStep)
				{
					MySort(matrixValues[i], matrixPosition[i], len);
				}

				for (int j = 0; j < len; j++)
				{
					sortedPosition[count] = matrixPosition[i][j];
					count++;
					if (iDiv < numOfDivision && count > nextStep)
					{
						divisionValues[iDiv++] = (float)matrixValues[i][j];
						nextStep += step;
					}
				}
			}

		}

		void MySort(vector<double> &values, vector<int> &index, int arrayLength)
		{
			struct by_number
			{
				bool operator()(data const &left, data const &right)
				{
					return left.number < right.number;
				}
			};

			if (items.size() != arrayLength)
				items.resize(arrayLength);

			for (int i = 0; i<arrayLength; i++)
			{
				items[i].number = values[i];
				items[i].index = index[i];
			}

			std::sort(items.begin(), items.end(), by_number());

			for (int i = 0; i<arrayLength; i++)
			{
				index[i] = items[i].index;
				values[i] = items[i].number;
			}
		}


	};
}