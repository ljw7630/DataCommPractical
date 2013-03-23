#include <iostream>
#include <fstream>
#include <ctime>
#include <sys\stat.h>
#include <map>
#include <vector>
#include <direct.h>
#include <algorithm>
#include <string>
using namespace std;


#ifdef _WIN32
	#define makedir(path, mode) _mkdir(path)
#else
	#define makedir(path, mode) mkdir(path, mode)
#endif

class FilePrinter{
private:
	ofstream output;
	string folderPath;
public:
	FilePrinter(){
		folderPath = makeFolder();
		string command("copy out.txt \"");
		command.append(folderPath);
		command.append("out.txt\"");
		system(command.c_str());
	}

	void writeError(string fileName, map<pair<int, int>, double> & res){
		
		output.open((folderPath + fileName).c_str());
		printError(res, output);
		output.close();
	}

	void writeResult(string fileName, vector<pair<double, pair<int, int> > > &res){
		
		output.open((folderPath + fileName).c_str());
		printResults(res, output);
		output.close();
	}

	void writeCorrectResult(string fileName, map<pair<int, int>, double> &res,
		double Xsource, double Ysource){
			
			output.open((folderPath + fileName).c_str());
			printCorrectResult(res, Xsource, Ysource, output);
			output.close();
	}

private:
	/// Random folder with timeStamp
	string makeFolder(){
		time_t t = time(NULL);
		tm* lctime = localtime(&t);
		char * cht = asctime(lctime);
		string s(cht);
		replace(s.begin(), s.end(), ':', '-');
		s.pop_back();
		s.push_back('/');
		int returnValue = makedir(s.c_str(), 777);

		return s;
	}

	void printError(map<pair<int, int>, double> &res, ostream &out = cout){
		for(auto itr = res.begin(); itr != res.end(); ++itr){
			out << (*itr).first.first << " " << (*itr).first.second 
				<< " " << (*itr).second << endl;
		}
	}
	
	void printResults(vector<pair<double, pair<int, int> > > &vec, ostream &out = cout){
		for (int i=0;i<vec.size();++i){
			printResult(vec[i], out);
		}
		out << endl;
	}

	void printResult(pair<double, pair<int, int> >&p, ostream &out = cout){
		out << p.second.first << " "
			<< p.second.second << " "
			<< p.first << endl;
	}

	void printResult(pair<pair<int, int>, double> &p, ostream &out = cout){
		out << p.first.first << " "
			<< p.first.second << " "
			<< p.second << endl;
	}

	void printCorrectResult(map<pair<int, int>, double> &res, double Xsource, double Ysource, ostream &out = cout){
		out << (int)Xsource << " " << (int)Ysource
		<<  " " << res[make_pair((int)Xsource,(int)Ysource)] << endl;
	}
};