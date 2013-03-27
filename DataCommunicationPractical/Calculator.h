#include <iostream>
#include <vector>
#include <utility>
#include <map>
#include <algorithm>
#include <fstream>
#include "FilePrinter.h"

using namespace std;

extern double Xsource, Ysource;

template<typename A, typename B>
class CompPairBySecondElement{
	bool operator() (const pair<A, B> &p1, const pair<A, B> &p2){
		return p1.second < p2.second;
	}
};

class Calculator{

private:
	int No_Residual_Surfaces;
	int No_Grps;
	int pairs;

	Calculator(){}

	template<typename A, typename B>
	static pair<B, A> flipPair(const pair<A, B> &p){
		return pair<B,A>(p.second, p.first);
	}

	template<typename A, typename B>
	static map<B, A> flipMap(const map<A,B> &src){
		map<B,A> dst;
		transform(src.begin(), src.end(), inserter(dst, dst.begin()), &Calculator::flipPair<A,B>);
		return dst;
	}

public:
	Calculator(int No_Residual_Surfaces, int No_Grps, int pairs){
		this->No_Residual_Surfaces = No_Residual_Surfaces;
		this->No_Grps = No_Grps;
		this->pairs = pairs;
	}

	bool operator() (pair<pair<int, int>, double>&p1, pair<pair<int, int>, double>&p2){
		return p1.second < p2.second;
	}

	void calculate(const vector< vector< vector<double> > > &Residual_Signal_Strength_Surfaces){
		cout << "Calculate" << endl;

		FilePrinter fp(this->pairs);

		vector< vector< pair<pair<int, int>, double> > > Residual_Signal_Strength_Surfaces_Vec;
		map<pair<int, int>, vector<pair<int, double> > > Residual_Signal_Strength_Surfaces_Map;
		map<pair<int, int>, pair<double, double> > Residual_Signal_Strength_Surfaces_Average_Variance_Map;

		map<pair<int, int>, double> rateAverageVariance;
		map<pair<int, int>, double> rateAverageSum;
		map<pair<int, int>, double> rateErrSum;
		map<pair<int, int>, double> rateIndexSum;
		map<pair<int, int>, double> rateWeightSum;

		pair< pair<int, int>, double > result;
		vector< pair<double, pair<int, int> > > results;

		getVectorOfResidualSignalStrengthSurfacesByXYPair(Residual_Signal_Strength_Surfaces, Residual_Signal_Strength_Surfaces_Vec);
		getMapOfResidualSignalStrengthSurfacesByXYPair(Residual_Signal_Strength_Surfaces_Vec, Residual_Signal_Strength_Surfaces_Map);
		getMapOfResidualSignalStrengthSurfacesVarianceByXYPair(Residual_Signal_Strength_Surfaces_Map, Residual_Signal_Strength_Surfaces_Average_Variance_Map);

		calAverageVarianceError(Residual_Signal_Strength_Surfaces_Average_Variance_Map, rateAverageVariance);
		calErrorSum(Residual_Signal_Strength_Surfaces_Vec, rateErrSum);
		calAverageErrorSum(rateAverageSum, rateErrSum);
		calIndexSum(Residual_Signal_Strength_Surfaces_Vec, rateIndexSum);
		calWeightedSum(Residual_Signal_Strength_Surfaces_Vec, rateWeightSum);

		findMin(rateAverageVariance, 5, results);
		fp.writeError("rateAverageVariance1.txt", rateAverageVariance);
		fp.writeResult("rateAverageVariance1Result.txt", results);
		fp.writeCorrectResult("rateAverateVariance1Result_Correct.txt", rateAverageVariance, Xsource, Ysource);

		findMin(rateAverageSum, 5, results);
		fp.writeError("rateAverageSum1.txt", rateAverageSum);
		fp.writeResult("rateAverageSum1Result.txt", results);
		fp.writeCorrectResult("rateAverageSum1Result_Correct.txt", rateAverageSum, Xsource, Ysource);

		findMin(rateErrSum, 5, results);
		fp.writeError("rateErrSum1.txt", rateErrSum);
		fp.writeResult("rateErrSum1Result.txt", results);
		fp.writeCorrectResult("rateErrSum1Result_Correct.txt", rateErrSum, Xsource, Ysource);


		findMin(rateIndexSum, 5, results);
		fp.writeError("rateIndexSum1.txt", rateIndexSum);
		fp.writeResult("rateIndexSum1Result.txt", results);
		fp.writeCorrectResult("rateIndexSum1Result_Correct.txt", rateIndexSum, Xsource, Ysource);

		findMin(rateWeightSum, 5, results);
		fp.writeError("rateWeightSum1.txt", rateWeightSum);
		fp.writeResult("rateWeightSum1Result.txt", results);
		fp.writeCorrectResult("rateWeightSum1Result_Correct.txt", rateWeightSum, Xsource, Ysource);

		// calErrorSum(Residual_Signal_Strength_Surfaces_Vec, rateErrSum);
		calIndexSum(Residual_Signal_Strength_Surfaces_Vec, 5, rateIndexSum);
		calWeightedSum(Residual_Signal_Strength_Surfaces_Vec, 5, rateWeightSum);

		findMin(rateIndexSum, 5, results);
		fp.writeError("rateIndexSum2.txt", rateIndexSum);
		fp.writeResult("rateIndexSum2Result.txt", results);
		fp.writeCorrectResult("rateIndexSum2Result_Correct.txt", rateIndexSum, Xsource, Ysource);

		findMin(rateWeightSum, 5, results);
		fp.writeError("rateWeightSum2.txt", rateWeightSum);
		fp.writeResult("rateWeightSum2Result.txt", results);
		fp.writeCorrectResult("rateWeightSum2Result_Correct.txt", rateWeightSum, Xsource, Ysource);

		// random print
	}

	void getVectorOfResidualSignalStrengthSurfacesByXYPair(
		const vector< vector< vector<double> > > &Residual_Signal_Strength_Surfaces,
		vector< vector< pair<pair<int, int>, double> > > & Residual_Signal_Strength_Surfaces_Vec){
			for(int tupleNo = 0; tupleNo<=No_Residual_Surfaces-1; tupleNo++){
				vector<pair<pair<int, int>, double> > vec;
				for(int xIdx = 0; xIdx <= No_Grps-1; xIdx++){
					for (int yIdx = 0; yIdx <= 49; yIdx++){
						vec.push_back(make_pair(make_pair(xIdx, yIdx), Residual_Signal_Strength_Surfaces[tupleNo][xIdx][yIdx]));
					}
				}
				Residual_Signal_Strength_Surfaces_Vec.push_back(vec);
			}

			for(auto i=0;i<Residual_Signal_Strength_Surfaces_Vec.size();++i){
				sort(Residual_Signal_Strength_Surfaces_Vec[i].begin(), Residual_Signal_Strength_Surfaces_Vec[i].end(), Calculator());
			}	
	}

	void getMapOfResidualSignalStrengthSurfacesByXYPair(
		vector< vector< pair<pair<int, int>, double> > > & Residual_Signal_Strength_Surfaces_Vec, 
		map<pair<int, int>, vector<pair<int, double> > > & Residual_Signal_Strength_Surfaces_Map){
			for(int x = 0; x<=No_Grps-1; ++x){
				for(int y = 0; y<=49;++y){
					Residual_Signal_Strength_Surfaces_Map.insert(make_pair(make_pair(x,y), vector<pair<int, double> >()));
				}
			}

			for(int tupleNo = 0; tupleNo<=No_Residual_Surfaces-1; tupleNo++){
				for(int i=0;i<Residual_Signal_Strength_Surfaces_Vec[tupleNo].size();++i){
					auto itr = Residual_Signal_Strength_Surfaces_Map.find(Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first);
					itr->second.push_back(make_pair(tupleNo, Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].second));
				}
			}
	}

	void getMapOfResidualSignalStrengthSurfacesVarianceByXYPair(
		map<pair<int, int>, vector<pair<int, double> > > &Residual_Signal_Strength_Surfaces_Map,
		map<pair<int, int>, pair<double, double> > &Residual_Signal_Strength_Surface_Variance_Map){
			for(auto itr = Residual_Signal_Strength_Surfaces_Map.begin(); itr != Residual_Signal_Strength_Surfaces_Map.end(); ++itr){
				double sum = 0;
				for(auto vitr = itr->second.begin(); vitr != itr->second.end(); ++vitr){
					sum += vitr->second;
				}
				double average = sum / itr->second.size();
				double variance = 0;
				for(auto vitr = itr->second.begin(); vitr != itr->second.end(); ++vitr){
					double diff = vitr->second - average;
					variance += (diff * diff) / itr->second.size();
				}
				Residual_Signal_Strength_Surface_Variance_Map.insert(make_pair(itr->first, make_pair(average, variance)));
			}
	}

	void calAverageVarianceError(map<pair<int, int>, pair<double, double> > & input,
		map<pair<int, int>, double> &rateAverageVarianceErr) {
			for(auto itr = input.begin(); itr != input.end(); ++itr){
				double result = itr->second.first * itr->second.second;
				rateAverageVarianceErr[itr->first] = itr->second.first + result;
			}
	}

	void calErrorSum(
		vector< vector< pair<pair<int, int>, double> > > & Residual_Signal_Strength_Surfaces_Vec, 
		map<pair<int, int>, double> &rateErrSum){
			rateErrSum.clear();
			for(int tupleNo = 0; tupleNo<=No_Residual_Surfaces-1; tupleNo++){
				for(auto i=0;i<Residual_Signal_Strength_Surfaces_Vec[tupleNo].size();++i){
					rateErrSum[Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first]+=Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].second;
				}
			}
	}

	void calAverageErrorSum(map<pair<int, int>, double > &rateAverageErrSum,
		map<pair<int, int>, double> &rateErrSum){
		rateAverageErrSum.clear();
		rateAverageErrSum[make_pair(0,0)] = rateErrSum[make_pair(0,0)];
			
		// X domains
		for(int yIdx = 0; yIdx <= 49; ++yIdx){
			for(int xIdx = 1; xIdx <= No_Grps-1; ++xIdx){
				rateAverageErrSum[make_pair(xIdx, yIdx)] = (rateErrSum[make_pair(xIdx, yIdx)] +
					rateErrSum[make_pair(xIdx-1, yIdx)])/2.0;
			}
		}
		// Y domains
		for (int xIdx=0; xIdx<=No_Grps-1; ++xIdx) {
			for (int yIdx=1; yIdx<=49; ++yIdx) {
				rateAverageErrSum[make_pair(xIdx,yIdx)] = (rateErrSum[make_pair(xIdx, yIdx)] +
					rateErrSum[make_pair(xIdx, yIdx-1)])/2.0;
			}
		}
	}

	void calIndexSum(
		vector< vector< pair<pair<int, int>, double> > > & Residual_Signal_Strength_Surfaces_Vec, 
		map<pair<int, int> , double> &rateIndexSum){
			rateIndexSum.clear();
			for(int tupleNo = 0; tupleNo <= No_Residual_Surfaces-1; ++tupleNo){
				for(int i = 0;i<Residual_Signal_Strength_Surfaces_Vec[tupleNo].size();++i){
					rateIndexSum[Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first]+=i;
				}
			}
	}

	void calIndexSum(
		vector< vector< pair<pair<int, int>, double> > > & Residual_Signal_Strength_Surfaces_Vec, 
		int k, map<pair<int, int> , double> &rateIndexSum){
			if(Residual_Signal_Strength_Surfaces_Vec.size()<0){
				return;
			}

			rateIndexSum.clear();
			map<pair<int, int>, int> indexCount;
			
			for(int tupleNo = 0; tupleNo <= No_Residual_Surfaces-1; ++tupleNo){
				for(int i=0;i<Residual_Signal_Strength_Surfaces_Vec[0].size();++i){
					indexCount[Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first]++;
					if(indexCount[Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first]>k){
						continue;
					}
					rateIndexSum[Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first] += i;
				}
			}
	}

	void calWeightedSum(
		vector< vector< pair<pair<int, int>, double> > > & Residual_Signal_Strength_Surfaces_Vec, 
		map<pair<int, int>, double> &rateWeightSum){
			rateWeightSum.clear();
			for(int tupleNo = 0; tupleNo<=No_Residual_Surfaces-1; tupleNo++){
				for(int i=0;i<Residual_Signal_Strength_Surfaces_Vec[tupleNo].size();++i){
					rateWeightSum[Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first]+= (i+1)*Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].second;
				}
			}
	}

	void calWeightedSum(
		vector< vector< pair<pair<int, int>, double> > > & Residual_Signal_Strength_Surfaces_Vec, 
		int k, map<pair<int, int>, double> &rateWeightSum){
			rateWeightSum.clear();
			if(Residual_Signal_Strength_Surfaces_Vec.size() == 0)
				return;
			map<pair<int, int>, int > indexCount;

			for(int tupleNo = 0; tupleNo <= No_Residual_Surfaces-1; ++tupleNo){
				for(int i=0;i<Residual_Signal_Strength_Surfaces_Vec[0].size();++i){
					indexCount[Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first]++;
					if(indexCount[Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first]>k){
						continue;
					}
					rateWeightSum[Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first]+= (i+1)*Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].second;
				}
			}
	}

	pair< pair<int, int>, double> findMin(map<pair<int, int> , double> & resultMap){
		pair<int, int> idx(0,0);
		double mini = 1e10;
		for(auto itr = resultMap.begin(); itr!=resultMap.end(); ++itr){
			if( mini > itr->second ){
				mini = itr->second;
				idx.first = itr->first.first;
				idx.second = itr->first.second;
			}
		}
		return make_pair(idx, mini);
	}

	void findMin(map<pair<int, int> , double> & resultMap, unsigned int k, vector<pair<double, pair<int,int> > > &arr){
		map<double, pair<int, int> > valueMap = flipMap(resultMap);
		arr.clear();	
		for(auto itr=valueMap.begin();itr!=valueMap.end();++itr){
			arr.push_back(*itr);
			if(arr.size()>=k)
				break;
		}
	}
};