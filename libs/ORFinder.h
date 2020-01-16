#pragma once
#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>
#include <cstdio>
#include <time.h>
#include <math.h>
#include <tuple>
#include <math.h>
#include <sstream>

#include <Node.h>
#include <WArray.h>

class ORF
{
public:
	pair<int, int> begin;
	pair<int, int> end;
	bool closed;
	int length;

	ORF() {};

	ORF(int n, int p)
	{
		begin = pair<int, int>(n, p);
		closed = false;
		length = 3;
	}
	void Close(int n, int p, int offset)
	{
		if (!closed)
		{
			end = pair<int, int>(n, p);
			closed = true;
			if (n == begin.first)
			{
				length = p - begin.second;
			}
			else
			{
				length += p-offset;
			}
		}
		
	}

	void AddLength(int n, int p, int offset)
	{
		if (!closed)
		{
			if (n == begin.first)
			{
				length = p - begin.second;
			}
			else
			{
				length += p - offset;
			}
		}
	}

	bool operator == (ORF &b)
	{
		if (this->begin != b.begin) return false;
		if (this->end != b.end) return false;
		if (this->closed != b.closed) return false;
		if (this->length != b.length) return false;
		return true;
	}

	bool operator < (ORF &b)
	{
		return(this->begin < b.begin);
	}

	bool operator > (ORF &b)
	{
		return(this->begin > b.begin);
	}
};

class CodonList //qNode
{
public:
	int ID;
	vector<char> markers;
	//'b' - start, 'e' -end, '\0' - nothing

	map<int, vector<tuple<char, int, int>>> Links = map<int, vector<tuple<char, int, int>>>();

	//map<start point for jump, vector<SIGN, ID, POS>>
	CodonList()
	{

	}


	CodonList(Node &n)
	{
		int sz = 0;
		int tnmp = n.str.size() - 2;
		if (tnmp > 0)
		{
			sz = n.str.size() - 2;
		}
		markers = vector<char>(sz);
		ID = n.ID;

		vector<int> lastpos = vector<int>(3);
		for (int i = 0; i < markers.size(); i += 3)
		{
			for (int j = 0; j < 3; j++)
			{
				if ((i + j + 2 < n.str.size()))
				{
					lastpos[j] = i + j + 3;
					string codon = n.str.substr(i + j, 3);
					if (IsStartCodon(codon))
					{
						markers[i + j] = 'b';
					}
					if (IsStopCodon(codon))
					{
						markers[i + j] = 'e';
					}
				}
			}
		}

		for (int j = 0; j < 3; j++)
		{
			if (lastpos[j] != n.str.size())
			{
				vector<tuple<string, int, int>> tmp_res = vector<tuple<string, int, int>>();
				n.GetTrackR(lastpos[j], 3, tmp_res);
				for (int i = 0; i < tmp_res.size(); i++)
				{
					string s_tmp = get<0>(tmp_res[i]);
					if (IsStartCodon(s_tmp))
					{
						Links[lastpos[j]].push_back(tuple<char, int, int>('b', get<1>(tmp_res[i]), get<2>(tmp_res[i])));
					}
					else if (IsStopCodon(s_tmp))
					{
						Links[lastpos[j]].push_back(tuple<char, int, int>('e', get<1>(tmp_res[i]), get<2>(tmp_res[i])));
					}
					else
					{
						Links[lastpos[j]].push_back(tuple<char, int, int>('\0', get<1>(tmp_res[i]), get<2>(tmp_res[i])));
					}
				}
			}
			else
			{
				for (auto t = n.Next.begin(); t != n.Next.end(); t++)
				{
					Links[lastpos[j]].push_back(tuple<char, int, int>('\0', t->second->ID, 0));
				}
			}
		}
	}

	bool IsStartCodon(string codon)
	{
		transform(codon.begin(), codon.end(), codon.begin(), ::tolower);
		if (codon == "atg")
			return true;
		return false;
	}
	bool IsStopCodon(string codon)
	{
		transform(codon.begin(), codon.end(), codon.begin(), ::tolower);
		if ((codon == "taa") || (codon == "tag") || (codon == "tga"))
			return true;
		return false;
	}


	map<pair<int, int>, vector<ORF>> Run(int offset, vector<ORF> ORFs = vector<ORF>())
 	{
		map<pair<int, int>, vector<ORF>> tmpData;
		if (ID == -1)
		{
			return tmpData;
		}
		int i = offset;
		for (; i< markers.size(); i += 3)
		{
			if (markers[i] == 'b')
			{
				ORF tmp(ID, i);
				ORFs.push_back(tmp);
			}
			if (markers[i] == 'e')
			{
				for (int j = 0; j < ORFs.size(); j++)
				{
 					ORFs[j].Close(ID, i, offset);
				}
			}
		}
		for (int j = 0; j < ORFs.size(); j++)
		{
			ORFs[j].AddLength(ID, i, offset);
		}

		for (int j = 0; j < Links[i].size(); j++)
		{
			vector<ORF> tmpORFs = ORFs;
			if (get<0>(Links[i][j]) == 'e')
			{
				for (int j = 0; j < tmpORFs.size(); j++)
				{
					tmpORFs[j].Close(ID, i, offset);
				}

			}
			else if (get<0>(Links[i][j]) == 'b')
			{
				ORF tmp(ID, i);;
				tmpORFs.push_back(tmp);
			}
			sort(tmpORFs.begin(), tmpORFs.end());
			tmpORFs.resize(abs(distance(unique(tmpORFs.begin(), tmpORFs.end()), tmpORFs.begin())));


			if (tmpORFs.size() != 0)
			{
				pair<int, int> t = pair<int, int>(get<1>(Links[i][j]), get<2>(Links[i][j]));
				tmpData[t] = tmpORFs;
			}
		}

		if ((tmpData.size() == 0)&&(ORFs.size() != 0))
		{
			tmpData[pair<int, int>(-2, 0)] = ORFs;
		}
		return tmpData;
	}

};

class ORFBlock
{
public:
	map<pair<int, int>, vector<ORF>> Data;
	map<int,vector<ORF>> finalORF;
	ORFBlock()
	{

	}

	vector<ORF> GetByPosition(int n, int p)
	{
		pair<int, int> key = pair<int, int>(n, p);
		vector<ORF> ref = Data[key];
		Data.erase(key);
		return ref;
	}

	void Add(map<pair<int, int>, vector<ORF>> tmpData)
	{
		for (auto t2 = tmpData.begin(); t2 != tmpData.end(); t2++)
		{
			this->Add(t2->first, t2->second);
		}
	}

	void Add(pair<int, int> np, vector<ORF> ORFs)
	{
		int q = np.first;
		for (int i = 0; i < ORFs.size(); i++)
		{
			if ((!ORFs[i].closed)&&(q!=-2))
			{
				Data[np].push_back(ORFs[i]);
			}
			else
			{
				finalORF[ORFs[i].length].push_back(ORFs[i]);
			}
		}
	}

	void Extract()
	{
		vector<pair<int, int>> todelete;
		for (auto t = Data.begin(); t != Data.end(); t++)
		{
			for (int i = 0; i < t->second.size(); i++)
			{
				if (t->second[i].closed)
				{
					finalORF[t->second[i].length].push_back(t->second[i]);
					t->second.erase(t->second.begin() + i);
					i--;
				}
			}

			if (t->second.size() == 0)
			{
				todelete.push_back(t->first);
			}
		}

		for (int i = 0; i < todelete.size(); i++)
		{
			Data.erase(todelete[i]);
		}
	}


	void Print(string outputfile)
	{
		ofstream f_out;
		f_out.open(outputfile);
		for (auto t = finalORF.begin(); t != finalORF.end(); t++)
		{
			f_out << t->first << '\t' << t->second.size() << endl;
		}
		f_out.close();
	}

	void InitRun(CodonList startCodonList)
	{
		vector<ORF> tmp;
		this->Add(startCodonList.Run(0, tmp));
		tmp.clear();
		this->Add(startCodonList.Run(1, tmp));
		tmp.clear();
		this->Add(startCodonList.Run(2, tmp));
	}

	void Run(map<int, CodonList> &ORFBody, ORFBlock &tmpBlock)
	{
		for (auto t = Data.begin(); t != Data.end(); t++)
		{
			map<pair<int, int>, vector<ORF>> tmp = ORFBody[t->first.first].Run(t->first.second, t->second);
			for (auto t2 = tmp.begin(); t2 != tmp.end(); t2++)
			{
				tmpBlock.Add(t2->first, t2->second);
			}
			
		}
	}
};

class ORFAnnotator
{
public:

	int startID;
	int endID;
	map<int, CodonList> CodonsBody;
	ORFBlock orfblock;

	ORFAnnotator(map<int, Node> &Body)
	{
		int sz = Body.size();
		for (auto t = Body.begin(); t != Body.end(); t++)
		{
			CodonsBody[t->first] = CodonList(t->second);
			if (t->second.End) endID = t->second.ID;
			if (t->second.ID == 0) startID = 0;
			cout << (t->first) << '\t' << sz << '\r';
		}
	}

	void Run(map<int, Node> &Body)
	{
		int sz = Body.size();
		vector<int> Nexts;
		vector<int> tmpNexts;
		Nexts.push_back(0);
		do
		{
			for (int i = 0; i < Nexts.size(); i++)
			{
				cout << Nexts[i] << '\t' << sz << '\r';
				for (int j = 0; j < 3; j++)
				{
					orfblock.Add(CodonsBody[Nexts[i]].Run(j, orfblock.GetByPosition(Nexts[i], j)));
				}
				for (auto t1 = Body[Nexts[i]].Next.begin(); t1 != Body[Nexts[i]].Next.end(); t1++)
				{
					tmpNexts.push_back(t1->second->ID);
				}
			}
			Nexts = tmpNexts;
			sort(Nexts.begin(), Nexts.end());
			Nexts.resize(abs(distance(unique(Nexts.begin(), Nexts.end()), Nexts.begin())));
			tmpNexts.clear();

		} while (Nexts.size() != 0);

		orfblock.Extract();
		orfblock.Print("outGraph_annot3.txt");
	}



};







