#include <iostream>
#include <vector>
#include <utility>
#include <map>
#include <algorithm>
using namespace std;

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
	Calculator(int No_Residual_Surfaces, int No_Grps){
		this->No_Residual_Surfaces = No_Residual_Surfaces;
		this->No_Grps = No_Grps;
	}

	bool operator() (pair<pair<int, int>, double>&p1, pair<pair<int, int>, double>&p2){
		return p1.second < p2.second;
	}

	void calculate(const vector< vector< vector<double> > > &Residual_Signal_Strength_Surfaces){
		vector< vector< pair<pair<int, int>, double> > > Residual_Signal_Strength_Surfaces_Vec;
		map<pair<int, int>, vector<pair<int, double> > > Residual_Signal_Strength_Surfaces_Map;

		map<pair<int, int>, double> rateErrSum;
		map<pair<int, int>, double> rateIndexSum;
		map<pair<int, int>, double> rateWeightSum;

		pair< pair<int, int>, double > result;
		vector< pair<double, pair<int, int> > > results;

		getVectorOfResidualSignalStrengthSurfacesByXYPair(Residual_Signal_Strength_Surfaces, Residual_Signal_Strength_Surfaces_Vec);
		getMapOfResidualSignalStrengthSurfacesByXYPair(Residual_Signal_Strength_Surfaces_Vec, Residual_Signal_Strength_Surfaces_Map);

		calErrorSum(Residual_Signal_Strength_Surfaces_Vec, rateErrSum);
		calIndexSum(Residual_Signal_Strength_Surfaces_Vec, 5, rateIndexSum);
		calWeightedSum(Residual_Signal_Strength_Surfaces_Vec, rateWeightSum);

		findMin(rateErrSum, 5, results);
		printResults(results);

		findMin(rateIndexSum, 5, results);
		printResults(results);

		findMin(rateWeightSum, 5, results);
		printResults(results);

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

			for(int i=0;i<Residual_Signal_Strength_Surfaces_Vec.size();++i){
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

	void calErrorSum(
		vector< vector< pair<pair<int, int>, double> > > & Residual_Signal_Strength_Surfaces_Vec, 
		map<pair<int, int>, double> &rateErrSum){
			for(int tupleNo = 0; tupleNo<=No_Residual_Surfaces-1; tupleNo++){
				for(int i=0;i<Residual_Signal_Strength_Surfaces_Vec[tupleNo].size();++i){
					rateErrSum[Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first]+=Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].second;
				}
			}
	}

	void calIndexSum(
		vector< vector< pair<pair<int, int>, double> > > & Residual_Signal_Strength_Surfaces_Vec, 
		map<pair<int, int> , double> &rateIndexSum){
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
			for(int tupleNo = 0; tupleNo<=No_Residual_Surfaces-1; tupleNo++){
				for(int i=0;i<Residual_Signal_Strength_Surfaces_Vec[tupleNo].size();++i){
					rateWeightSum[Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].first]+= (i+1)*Residual_Signal_Strength_Surfaces_Vec[tupleNo][i].second;
				}
			}
	}

	void calWeightedSum(
		vector< vector< pair<pair<int, int>, double> > > & Residual_Signal_Strength_Surfaces_Vec, 
		int k, map<pair<int, int>, double> &rateWeightSum){
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

	void findMin(map<pair<int, int> , double> & resultMap, int k, vector<pair<double, pair<int,int> > > &arr){
		map<double, pair<int, int> > valueMap = flipMap(resultMap);
		arr.clear();	
		for(auto itr=valueMap.begin();itr!=valueMap.end();++itr){
			arr.push_back(*itr);
			if(arr.size()>=k)
				break;
		}
	}

	void printResults(vector<pair<double, pair<int, int> > > &vec){
		for (int i=0;i<vec.size();++i){
			printResult(vec[i]);
		}
		cout << endl;
	}

	void printResult(pair<double, pair<int, int> >&p){
		cout << "x index: " << p.second.first 
			<< ", y index: " << p.second.second
			<< ", value: " << p.first << endl;
	}

	void printResult(pair<pair<int, int>, double> &p){
		cout << "x index: " << p.first.first 
			<< ", y index: " << p.first.second
			<< ", value: " << p.second << endl;
	}
};