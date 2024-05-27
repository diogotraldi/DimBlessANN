// DimBlessANN.h

#pragma once
#include "SearchSpaceANN.h"
#include "KeepElimInfo.h"
#include <map>
using namespace System;
using namespace System::Collections::Generic;

namespace DimBless {

	public ref class DimBlessANN
	{

	private:
		SearchSpaceANN * searchSpaceDynamic;
		SearchSpaceANN * searchSpaceStatic;

		double _dynamicSpaceRatio;

		float * _searchVect;
		int _spaceDim;

		std::map<__int64, float *> *_spaceDynamic = nullptr;
		std::map<__int64, float *> *_spaceStatic = nullptr;

		int _foundIdLength = 0;
		__int64 * FoundID;	

		float * DimBlessANN::GetNative(cli::array<float> ^ array)
		{
			//return (float*)&array[0];
			try 
			{
				pin_ptr<float> array_pin = &array[0];
				return (float*)array_pin;
			}
			catch (...) 
			{
				return 0;
			}
		}

		void SetDynamicRatio()
		{
			int staticLen = searchSpaceStatic->Size();
			int dynamicLen = searchSpaceDynamic->Size();

			_dynamicSpaceRatio = (staticLen == 0) ? 1 : (double)dynamicLen / (dynamicLen + staticLen);

			if (_foundIdLength != staticLen + dynamicLen)
			{
				if (FoundID != nullptr)
					delete[] FoundID;

				FoundID = new __int64[staticLen + dynamicLen];
			}
		}

	public:
		DimBlessANN(int spaceDim, int numOfSpaces)
		{
			_spaceDim = spaceDim;
			_searchVect = new float[spaceDim];

			_spaceDynamic = new std::map<__int64, float *>();
			_spaceStatic = new std::map<__int64, float *>();

			searchSpaceStatic = new SearchSpaceANN(spaceDim, numOfSpaces);
			searchSpaceDynamic = new SearchSpaceANN(spaceDim, numOfSpaces, searchSpaceStatic->_randomIndex);
		}

		~DimBlessANN()
		{
			ClearDynamicVectors();
			searchSpaceDynamic->Clear(); 

			ClearStaticVectors();
			searchSpaceStatic->Clear();

			if (_searchVect != nullptr)
				delete[] _searchVect;

			if (FoundID != nullptr)
				delete[] FoundID;

			delete _spaceDynamic;
			delete _spaceStatic;
			delete searchSpaceDynamic;
			delete searchSpaceStatic;
		}

		void ClearSpace()
		{
			ClearDynamicVectors();
			searchSpaceDynamic->ClearSpaceMemory();

			ClearStaticVectors();
			searchSpaceStatic->ClearSpaceMemory();
		}

		void ClearDynamicVectors()
		{
			std::map<__int64, float *> &spaceDynamic = *_spaceDynamic;
			for (auto it = spaceDynamic.begin(); it != spaceDynamic.end(); it++)
			{
				delete[] it->second;
			}
			spaceDynamic.clear();
		}

		void ClearStaticVectors()
		{
			std::map<__int64, float *> &spaceStatic = *_spaceStatic;
			for (auto it = spaceStatic.begin(); it != spaceStatic.end(); it++)
			{
				delete[] it->second;
			}
			spaceStatic.clear();
		}

		void AddDynamic(cli::array<float>^ vect, Int64 id)
		{
			int len = vect->Length;

			float * vectNative = new float[len];
			for (int i = 0; i < len; i++)
			{
				vectNative[i] = vect[i];
			}

			std::map<__int64, float *> &spaceDynamic = *_spaceDynamic;
			spaceDynamic[(__int64)id] = vectNative;
		}

		void AddStatic(cli::array<float>^ vect, Int64 id)
		{
			int len = vect->Length;

			float * vectNative = new float[len];
			for (int i = 0; i < len; i++)
			{
				vectNative[i] = vect[i];
			}

			std::map<__int64, float *> &spaceStatic = *_spaceStatic;
			spaceStatic[(__int64)id] = vectNative;
		}

		int StaticVectorCount()
		{
			return (int)_spaceStatic->size();
		}

		int StaticSpaceLen()
		{
			return searchSpaceStatic->Size();
		}

		int DynamicSpaceLen()
		{
			return searchSpaceDynamic->Size();
		}

		int DynamicVectorCount()
		{
			return (int)_spaceDynamic->size();
		}

		void IndexSearchSpace(int bitsPerSample)
		{
			IndexSearchSpaceDynamic(bitsPerSample);
			IndexSearchSpaceStatic(bitsPerSample);
		}

		void IndexSearchSpaceDynamic(int bitsPerSample)
		{
			std::map<__int64, float *> &spaceDynamic = *_spaceDynamic;
			std::map<__int64, float *> &spaceStatic = *_spaceStatic;

			searchSpaceDynamic->IndexSpace(spaceDynamic, bitsPerSample);

			SetDynamicRatio();
		}

		void IndexSearchSpaceStatic(int bitsPerSample)
		{
			std::map<__int64, float *> &spaceStatic = *_spaceStatic;
			std::map<__int64, float *> &spaceDynamic = *_spaceDynamic;

			searchSpaceStatic->IndexSpace(spaceStatic, bitsPerSample, searchSpaceDynamic->_roundMem, searchSpaceDynamic->_mainCol);

			SetDynamicRatio();
		}

		void SetParameters(int colNum, int colMin, int colStepCheck, float limiarStepUp, float limiarStepDown, int stepMonitor)
		{
			searchSpaceDynamic->SetDynamicParameters(colNum, colMin, colStepCheck, limiarStepUp, limiarStepDown, stepMonitor);
		}

		void IndexSearchSpace(Dictionary<Int64, cli::array<float>^>^ spaceVectors, int bitsPerSample, double dynamicSpaceRatio)
		{
			double count = 0;

			std::map<__int64, float *> &spaceDynamic = *_spaceDynamic;
			std::map<__int64, float *> &spaceStatic = *_spaceStatic;

			for each(KeyValuePair<Int64, cli::array<float>^> pair in spaceVectors)
			{
				if (count > 1.0)
				{
					spaceDynamic[(__int64)pair.Key] = GetNative(pair.Value);
					count -= 1.0;
				}
				else
				{
					spaceStatic[(__int64)pair.Key] = GetNative(pair.Value);
				}

				count += dynamicSpaceRatio;
			}

			_dynamicSpaceRatio = dynamicSpaceRatio;

			if (FoundID != nullptr)
				delete[] FoundID;

			if (_dynamicSpaceRatio == 1)
			{
				searchSpaceDynamic->IndexSpace(spaceDynamic, bitsPerSample);
				FoundID = new __int64[spaceDynamic.size()];
			}
			else
			{
				searchSpaceStatic->IndexSpace(spaceStatic, bitsPerSample);
				searchSpaceDynamic->IndexSpace(spaceDynamic, bitsPerSample, searchSpaceStatic->_roundMem, searchSpaceStatic->_mainCol);

				FoundID = new __int64[spaceDynamic.size() + spaceStatic.size()];
			}

		}
	
		int Search(cli::array<Single>^ vect, int bruteTarget)
		{
			bool allZero = true;
			for (int i = 0; i < _spaceDim; i++)
			{
				_searchVect[i] = vect[i];
				if (_searchVect[i] != 0)
					allZero = false;
			}

			if (allZero)
				return 0;

			int numOfFound = 0;
			if (_dynamicSpaceRatio == 1)
			{
				numOfFound = searchSpaceDynamic->SearchDynamic(_searchVect, false, FoundID, bruteTarget);
			}
			else
			{				
				int bruteTargetDynamic = (int)((bruteTarget * _dynamicSpaceRatio) + 0.5);
				if (bruteTargetDynamic < 0)
					bruteTargetDynamic = 1;

				int nDynamic = searchSpaceDynamic->SearchDynamic(_searchVect, false, FoundID, bruteTargetDynamic); // round value when add 0.5
				int nStatic = searchSpaceStatic->SearchStaticW(searchSpaceDynamic->VectW, searchSpaceDynamic->OrderW, searchSpaceDynamic->ReptitionControl, searchSpaceDynamic->limiar, searchSpaceDynamic->colNum, FoundID + nDynamic);
				
				numOfFound = nDynamic + nStatic;
			}

			return numOfFound;
		}

		int SearchDebug(cli::array<Single>^ vect, EliminationLog^ elimInfo, int bruteTarget)
		{
			for (int i = 0; i < _spaceDim; i++)
				_searchVect[i] = vect[i];

			int nDynamic = searchSpaceDynamic->SearchDynamic(_searchVect, true, FoundID, bruteTarget);
			int nStatic = searchSpaceStatic->SearchStaticW(searchSpaceDynamic->VectW, searchSpaceDynamic->OrderW, searchSpaceDynamic->ReptitionControl, searchSpaceDynamic->limiar, searchSpaceDynamic->colNum, FoundID + nDynamic);

			int numOfFound = nDynamic + nStatic;

			elimInfo->ColElim->Clear();
			elimInfo->LimiarElim->Clear();
			elimInfo->CountElim->Clear();
			for (size_t i = 0; i < searchSpaceDynamic->keepColElim.size(); i++)
			{
				elimInfo->ColElim->Add(searchSpaceDynamic->keepColElim[i]);
				elimInfo->LimiarElim->Add(searchSpaceDynamic->keepLimiarElim[i]);
				elimInfo->CountElim->Add(searchSpaceDynamic->keepCountElim[i]);
			}

			elimInfo->RealColNum->Clear();
			elimInfo->RealColNum->Add(searchSpaceStatic->iTop + searchSpaceStatic->iBot + 1);

			return numOfFound;
		}

		__int64 * GetFoundID()
		{
			return FoundID;
		}

		cli::array<Int64>^ GetSortedSpaceID()
		{
			vector<__int64> idDynamic = searchSpaceDynamic->GetSortedSpaceID();
			vector<__int64> idStatic = searchSpaceStatic->GetSortedSpaceID();

			int len = (int)(idDynamic.size() + idStatic.size());

			cli::array<Int64>^ sortedId = gcnew cli::array<Int64>(len);

			int position = 0;
			for (int i = 0; i < idDynamic.size(); i++)
				sortedId[position++] = idDynamic[i];
			
			for (int i = 0; i < idStatic.size(); i++)
				sortedId[position++] = idStatic[i];

			return sortedId;
		}
	};
}
