#include <deque>

using namespace System;
using namespace System::Collections::Generic;

namespace DimBless
{
	public ref struct EliminationLog
	{
		List<double>^ ColElim = gcnew List<double>();
		List<double>^ LimiarElim = gcnew List<double>();
		List<double>^ CountElim = gcnew List<double>();
		List<double>^ RealColNum = gcnew List<double>();
	};
}