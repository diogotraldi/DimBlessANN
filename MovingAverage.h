#include <deque>

using namespace std;

namespace DimBless
{
class MovingAverage 
{
private:
    deque<double> q;
    int queue_size;
public:
    /** Initialize your data structure here. */
    MovingAverage(int size)
	{
		queue_size = size;
    }
    
    void AddValue(double val)
	{
		if (q.size() < queue_size)
		{
			q.push_back(val);
		}
        else if(q.size() == queue_size)
        {
            q.pop_front();
        }
    }

	double GetAverage() 
	{
		double sum = 0;

		size_t len = q.size();
		if (len == 0)
			return 0;

		for (size_t i = 0; i < len; i++)
		{
			sum += q[i];
		}

		return sum / len;
	}

	double StandardDeviation()
	{
		size_t len = q.size();

		double average = GetAverage();

		double variance = 0;
		for (size_t i = 0; i < len; i++)
		{
			double diff = q[i] - average;
			variance += (diff * diff);
		}

		return sqrt((variance / len));
	}

	size_t size()
	{
		return q.size();
	}
};

}