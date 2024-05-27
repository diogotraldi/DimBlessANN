
namespace DimBless
{
	public class Walsh
	{
		int n;
		int i;
		int j;
		int j2;
		int jd;
		int js;
		int l;
		int m;
		int n2;
		int nx;
		int ny;
		int nz;
		int nzi;
		int nzn;
		float *y;

	public:
		void Init(int n);
		void Close();

		void fwt(float x[]);

	protected:
		inline int i4_log_2(int i);
		inline int i4_power(int i, int j);
		inline void r8vec_copy(int n, float a1[], float a2[]);
	};
}