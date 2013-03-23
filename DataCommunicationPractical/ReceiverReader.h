#include <iostream>
#include <vector>
#include <string>
#include <fstream>
using namespace std;

class ReceiverReader{
	
	vector<vector<double> > Receiver_Tuple_Array;

public:
	ReceiverReader(int No_Residual_Surfaces, int num = 7){

		// number of receiver pairs * 7
		Receiver_Tuple_Array = vector<vector<double> >(No_Residual_Surfaces, vector<double>(num));

		readFromFile(No_Residual_Surfaces);
	}

	vector<vector<double> > & getReceiverTupleArray(){
		return Receiver_Tuple_Array;
	}

private:
	void readFromFile(const int num, string fileName = "out.txt"){
		ifstream fs(fileName);
		double x, y, m;
		int index = 0;
		if(fs.is_open()){
			while(fs.good()){
				fs >> x >> y >> m;
				Receiver_Tuple_Array[index][1] = x;
				Receiver_Tuple_Array[index][2] = y;
				Receiver_Tuple_Array[index][3] = m;

				fs >> x >> y >> m;
				Receiver_Tuple_Array[index][4] = x;
				Receiver_Tuple_Array[index][5] = y;
				Receiver_Tuple_Array[index][6] = m;
				index++;

				if(num == index){
					break;
				}
			}
			fs.close();
		}
	}
public: 
	void print(){
		if(Receiver_Tuple_Array.size()==0) return;
		for(int i=0;i<Receiver_Tuple_Array.size();++i){
			for(int j=1;j<Receiver_Tuple_Array[0].size();++j){
				cout << Receiver_Tuple_Array[i][j] << " ";
			}
			cout << endl;
		}
	}
};