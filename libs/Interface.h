 #pragma once
#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#include <time.h>
#include <map>
#include <set>
#include <math.h>
#include <unordered_map>
#include <tuple>
#include <math.h>
#include <sstream>
#include <Graph.h>
#include <VarReader.h>
#include <Fasta.h>
#include <thread>
#include <future>

using namespace std;


struct  Results
{
	int total = 0;
	int aligned = 0;

	Results operator + (const Results &v1)
	{
		Results res;
		res.total = v1.total + this->total;
		res.aligned = v1.aligned + this->aligned;
		return res;
	}
};

class Worker
{
public:
	float coverage;
	future<Results> result;
	vector<string> vec_aligned_reads;
	bool enabled;
	int ID;
	string curr_read = "";
	Worker(int nID, float cov)
	{
		coverage = cov;
		ID = nID;
	}

	void Run(Graph *VG, ifstream *fin, mutex *mut, ofstream *fout, mutex *outmut)
	{
		enabled = false;
		result = async(launch::async, &Worker::ThreadTask, this, VG, fin, mut, fout, outmut);
	}

	Results ThreadTask(Graph *VG, ifstream *fin, mutex *mut, ofstream *fout, mutex *outmut)	
	{
		Results res;
		int aligned_reads=0;
		int total_reads=0;
		string read = "";
		string header = "";
		mut->lock();
		while (!fin->eof())
		{
			getline(*fin, header);
			getline(*fin, read);
			mut->unlock();
			total_reads++;
			string res_read = reverse(read);
			curr_read = read;
			int tmp_res = VG->SFinder(read, coverage);
			int tmp_res2 = VG->SFinder(res_read, coverage);
			if (tmp_res == 0)
			{
				outmut->lock();
				(*fout) << header << endl;
				(*fout) << read << endl;
				outmut->unlock();
				aligned_reads++;
			}
			if (tmp_res2 == 0)
			{
				outmut->lock();
				(*fout) << header<<"_rev"<< endl;
				(*fout) << read << endl;
				outmut->unlock();
				aligned_reads++;
			}
			mut->lock();
		}
		mut->unlock();
		res.aligned = aligned_reads;
		res.total = total_reads;
		return res;
	}

	Results Finalize()
	{
		future_status res_stat = result.wait_for(chrono::microseconds(100));
		if (res_stat == future_status::ready)
		{
			cout<<"ID = "<<ID<<endl;
			cout<<"status: ready"<<endl;
			return result.get();
		}
		if (res_stat == future_status::timeout)
		{
			cout<<"ID = "<<ID<<endl;
			cout<<"status: timeout"<<endl;
			cout<<"Read: "<<curr_read<<endl;
			result.wait();
			return result.get();
		}
		if (res_stat == future_status::deferred)
		{
			cout<<"ID = "<<ID<<endl;
			cout<<"status: deferred"<<endl;
			Results res;
			return res;
		}
	}

};


class ThreadPool
{
public:
	ifstream fin;
	ofstream fout;
	vector<Worker> threads;
	mutex mut;
	mutex outmut;

	ThreadPool(string path, string outpath, float cov, int numthreads)
	{
		fin = ifstream(path);
		fout = ofstream(outpath);
		for (int i = 0; i < numthreads; i++)
		{
			threads.push_back(Worker(i,cov));
		}
		cout<<"TPool construct: numthreads = "<<threads.size()<<endl;
	}

	void Run(Graph &VG)
	{
		for (int i = 0; i < threads.size(); i++)
		{
			threads[i].Run(&VG, &fin, &mut, &fout, &outmut);
		}
	}

	Results GetResults()
	{
		Results res;
		res.aligned = 0;
		res.total = 0;

		for (int i = 0; i < threads.size(); i++)
		{
			res = res + threads[i].Finalize();
		}
		return res;
	}

	void PrintTotalAlignedReads(string output)
	{
		ofstream f_out;
		f_out.open(output);
		int count = 0;
		for (int i = 0; i < threads.size(); i++)
		{
			for (int j = 0; j < threads[i].vec_aligned_reads.size(); j++)
			{
				//f_out << ">R" << count << endl;
				//f_out << threads[i].vec_aligned_reads[j]<<endl;
				count++;
			}
		}

		f_out<< "TotalResult: " << count << " reads aligned" << endl;
		f_out.close();
	}

};


int SRAalign(Graph &VG, string path, string outpulog, string outputreads, float cov, int numthreads)
{
	cout << "SRAAlign: numthreads " << numthreads << endl;
	ThreadPool TPool(path, outputreads, cov, numthreads);
	cout<<"SRAAlign: TPool created"<<endl;
	TPool.Run(VG);
	cout<<"SRAAlign: TPool runned"<<endl;
	Results res = TPool.GetResults();
	cout<<"SRAAlign: TPool results got"<<endl;

	ofstream f_out;

	f_out.open(outpulog);
	f_out << "Total:" << res.total << '\n';
	f_out << "Aligned:" << res.aligned << '\n';
	//f_out << "Nonaligned:" << TPool.non_aligned_reads << '\n';
	f_out.close();

	//TPool.PrintTotalAlignedReads(outputreads);
	cout<<"SRAAlign: file was writed"<<endl;
	

	return 0;
}


void RemoveDupsFromFile(string path, string output)
{
	map<string, bool> checked;

	ifstream fin(path);
	ofstream f_out(output);
	string buf;
	string title;
	string res_buf;

	while ((!fin.eof()))
	{
		getline(fin, title);
		getline(fin, buf);
		res_buf = reverse(buf);

		if ((checked.find(buf) == checked.end()) &&
			(checked.find(res_buf) == checked.end()))
		{
			checked[buf] = true;
			checked[res_buf] = true;

			f_out << title << endl;
			f_out << buf << endl;
		}

		
	}
	fin.close();
	f_out.close();

}

/*int SRAalign(Graph &VG, string path, string output)
{
	map<string, bool> checked;
	int count = 0;
	int total = 0;
	int failed_r = 0;
	int failed_s = 0;

	string buf;
	string res_buf;
	ifstream fin(path);
	while (!fin.eof())
	{

		getline(fin, buf);
		getline(fin, buf);
		res_buf = reverse(buf);



		if ((checked.find(buf) == checked.end()) &&
			(checked.find(res_buf) == checked.end()))
		{
			checked[buf] = true;
			checked[res_buf] = true;
			total++;
			int res = VG.SFinder(buf);
			if (res != 0)
			{

				res = VG.SFinder(res_buf);
				if (res != 0)
				{
					count++;
				}
			}
			else
			{
				count++;
			}
		}
	}
	fin.close();

	ofstream f_out;

	f_out.open(output);
	f_out << "Total:" << total << '\n';
	f_out << "Aligned:" << count << '\n';
	f_out << "Proportion:" << (float)(count / total) << '\n';
	f_out.close();

	return 0;

}*/