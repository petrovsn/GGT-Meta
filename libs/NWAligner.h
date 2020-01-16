#pragma once
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
#include <algorithm>

#include <GenerHash.h>
#include <Node.h>
#include <BubbleIndex.h>




void printNWmart(string ref, string read, vector<vector<int>> NWMatrix1)
{
	cout << ref.size() << '/' << NWMatrix1.size() << endl;
	cout << read.size() << '/' << NWMatrix1[0].size() << endl;
	cout << '\t' << '\t';
	for (int i = 0; i < ref.size(); i++)
	{
		cout << ref[i] << '\t';
	}
	cout << endl;
	cout << '\t';
	for (int i = 0; i < NWMatrix1.size(); i++)
	{
		cout << NWMatrix1[i][0] << '\t';
	}
	cout << endl;
	for (int j = 1; j < read.size(); j++)
	{
		cout << read[j - 1] << '\t';
		for (int i = 0; i < NWMatrix1.size(); i++)
		{
			cout << NWMatrix1[i][j] << '\t';
		}
		cout << endl;
	}
}

class Weights
{
public:
	map<char, map<char, int>> CharMatrix;
	int gap;

	Weights()
	{
		CharMatrix = InitCharMatrix();
		gap = -5;
	}

	map<char, map<char, int>> InitCharMatrix()
	{
		map<char, map<char, int>> CharMatrix;
		for (char c1 : { 'a', 'c', 'g', 't', 'A', 'C', 'G', 'T' })
		{
			for (char c2 : { 'a', 'c', 'g', 't', 'A', 'C', 'G', 'T' })
			{
				if (tolower(c1) == tolower(c2))
				{
					CharMatrix[c1][c2] = 2;
				}
				else
				{
					CharMatrix[c1][c2] = -2;
				}
			}
		}
		return CharMatrix;
	}
};

Weights WGlob = Weights();

class GVariation
{
public:

	int ID1;
	int pos1;
	int ID2;
	int pos2;

	string alt;
	bool active = false;

	GVariation()
	{

	}
	
	GVariation(int eID, int epos)
	{
		ID2 = eID;
		pos2 = epos;
		active = true;
	}

	void Extend(char c)
	{
		alt = c + alt;
	}

	void Extend()
	{
	}

	void Finish(int sID, int spos)
	{
		ID1 = sID;
		pos1 = spos;
		active = false;
	}
};


struct NWpart
{
	vector<GVariation> vars;
	GVariation lastvar;

	string alnread = "";
	string alnref = "";

	int ID1 = -1;
	int pos1 = -1;
	int ID2 = -1;
	int pos2 = -1;

	int lastloc = -1;
	int score = 0;
	int ID = -1;

	int lastalignedID = -1;
	int lastalignedpos = -1;

	bool finished = false;
};




class NWNode
{
public:
	
	string ref;
	int ID;
	int startpos = 0;

	vector<int> NWNext;
	
	bool status = false;

	NWNode()
	{
		ID = -1;
		ref = "empty";
	}
	 
	NWNode(string &reference, int newID)
	{
		ID = newID;
		ref = reference;
	}

	NWNode(Node &N, vector<int>& links, map<int, Node> &Body, int N1, int pos1, int N2, int pos2, int lenread, bool first = false)
	{
		if (first)
		{
			startpos = pos1;
		}

		ID = N.ID;

		if ((ID == N1) and (ID == N2))//in one node
		{
			if (lenread*3 < (pos2 - pos1 + 1))
			{
				return;
			}
			ref = N.str.substr(pos1, pos2 - pos1+1);
			status = true;
		}
		else if ((ID == N1) and (ID != N2))//first node of subgraph
		{
			int tail = N.str.size() - pos1;
			ref = N.str.substr(tail);
			status = true;
		}
		else if ((ID != N1) and (ID == N2))//last node of subraph
		{
			if (lenread*3 < pos2)
			{
				return;
			}
			ref = N.str.substr(0,pos2);
			status = true;
		}
		else
		{
			ref = N.str;
			status = true;
		}

		for (auto p1 : N.Next)
		{
			bool internode = (*(p1.second)).IsInside(N1, N2, Body);
			bool finalnode = ((*(p1.second)).ID == N2);

			if (internode or finalnode)
			{
				links.push_back((*(p1.second)).ID);
				NWNext.push_back((*(p1.second)).ID);
			}
		}
	}

	void Link()
	{

	}

	vector<int> FirstIds;
	vector<vector<int>> NWMatrix;
	

	void InitMatrix(string &read, map<int, vector<int>> prevcolmns = map<int, vector<int>>())
	{
		FirstIds = vector<int>();
		
		int Ilen = ref.size() + 1;
		int Jlen = read.size() + 1;

		NWMatrix = vector<vector<int>>(Ilen);

		for (int i = 0; i < NWMatrix.size(); i++)
		{
			NWMatrix[i] = vector<int>(Jlen);
			for (int j = 0; j < NWMatrix[i].size(); j++)
			{
				NWMatrix[i][j] = 0;
			}
		}

		//fill first column
		if (prevcolmns.size() > 0) 
		{
			vector<int> previds;

			for (map<int, vector<int>>::iterator it = prevcolmns.begin(); it != prevcolmns.end(); ++it)
			{
				previds.push_back(it->first);
			}

			for (int j = 0; j < read.size() + 1; j++)
			{
				int curr_id = previds[0];
				int curr_max = prevcolmns[curr_id][j];

				for (int k = 0; k < previds.size(); k++)
				{
					if (prevcolmns[previds[k]][j] > curr_max)
					{
						curr_id = previds[k];
						curr_max = prevcolmns[previds[k]][j];
					}
				}

				FirstIds.push_back(curr_id);
				NWMatrix[0][j] = curr_max;
			}
		}  //fill first
		else
		{
			NWMatrix[0][0] = 0;
			for (int j = 1; j < NWMatrix[0].size(); j++)
			{
				NWMatrix[0][j] = NWMatrix[0][j - 1] + WGlob.gap;
			}
		}

		//fill first raw
		for (int i = 1; i < NWMatrix.size(); i++)
		{
			NWMatrix[i][0]=NWMatrix[i-1][0]+ WGlob.gap;
		}

		//fill other matrix
		for (int i = 1; i < NWMatrix.size(); i++)
		{
			for (int j = 1; j < NWMatrix[i].size(); j++)
			{
				vector<int> vals = { NWMatrix[i - 1][j - 1] + WGlob.CharMatrix[ref[i]][read[j]],
					NWMatrix[i - 1][j] + WGlob.gap ,
					NWMatrix[i][j - 1] + WGlob.gap };

				vector<int>::iterator res = max_element(vals.begin(), vals.end());
				int tmpidx = distance(vals.begin(), res);
				NWMatrix[i][j] = vals[tmpidx];
			}
		}
	}

	void AddVar(char readsimb, NWpart &head)
	{

	}

	NWpart NWTrace(string &read, NWpart head)
	{
		int i_pos = NWMatrix.size() - 1;
		int j_pos = NWMatrix.back().size() - 1;
		int score = head.score;

		if (head.lastloc!=-1)
		{
			j_pos = head.lastloc;
		}

		string alignedref = head.alnref;
		string alignedread = head.alnread;

		do
		{
			vector<int> vals = { NWMatrix[i_pos - 1][j_pos - 1], //match
								 NWMatrix[i_pos - 1][j_pos], //left
								 NWMatrix[i_pos][j_pos - 1] }; //up

			vector<int>::iterator res = max_element(vals.begin(), vals.end());
			int tmpidx = distance(vals.begin(), res);

			if (tmpidx == 0)
			{
				alignedref = ref[i_pos-1] + alignedref;
				alignedread = read[j_pos-1] + alignedread;

				if (tolower(ref[i_pos - 1]) == tolower(read[j_pos - 1]))
				{
					score = score + 1;
					head.lastalignedID = ID;
					head.lastalignedpos = i_pos - 1 + startpos;

					if (head.lastvar.active) //if variation is over
					{
						head.lastvar.Finish(ID, i_pos - 1 + startpos);
						head.vars.push_back(head.lastvar);
					}
				}
				else
				{
					if (head.lastvar.active) //if variation active extende, else - start new and extend. 
					{
						head.lastvar.Extend(read[j_pos - 1]);
					}
					else
					{
						head.lastvar = GVariation(head.lastalignedID, head.lastalignedpos);
						head.lastvar.Extend(read[j_pos - 1]);
					}
				}
				i_pos = i_pos - 1;
				j_pos = j_pos - 1;
				
				

			}
			else if (tmpidx == 1)
			{
				if (head.lastvar.active) //if variation active extende, else - start new and extend. 
				{
					head.lastvar.Extend();
				}
				else
				{
					head.lastvar = GVariation(head.lastalignedID, head.lastalignedpos);
					head.lastvar.Extend();
				}

				alignedref = ref[i_pos-1] + alignedref;
				alignedread = '-' + alignedread;

				i_pos = i_pos - 1;
			}
			else if (tmpidx == 2)
			{
				if (head.lastvar.active) //if variation active extende, else - start new and extend. 
				{
					head.lastvar.Extend(read[j_pos - 1]);
				}
				else
				{
					head.lastvar = GVariation(head.lastalignedID, head.lastalignedpos);
					head.lastvar.Extend(read[j_pos - 1]);
				}



				alignedread = read[j_pos-1] + alignedread;
				alignedref = '-' + alignedref;
				j_pos = j_pos - 1;
			}


		} while ((i_pos != 0) and (j_pos != 0));

		NWpart res;
		
		if (j_pos == 0)
		{
			head.finished = true;
		}
		head.score = score + head.score;
		head.alnread = alignedread;
		head.alnref = alignedref;
		head.lastloc = j_pos;
		

		if (FirstIds.size() != 0)
		{
			head.ID = FirstIds[j_pos];
		}
		return head;
	}
	 

	
};

class NWAligner
{
public:
	map<int, NWNode> NWBody;
	int lastID;
	bool loaded = false;

	int nwID1 = -1;
	int nwpos1 = -1;
	int nwID2 = -1;
	int nwpos2 = -1;

	NWAligner()
	{

	}
	void Init(int ID1, int pos1, int ID2, int pos2, map<int, Node> &Body, string &read)
	{
		loaded = false;
		NWBody.clear();
		lastID = ID2;
		vector<int> links;
		NWNode startNode = NWNode(Body[ID1], links, Body, ID1, pos1, ID2, pos2, read.size(), true);
		
		if (!startNode.status) return;

		map<int, vector<int>> prevIds;
		for (int j = 0; j < startNode.NWNext.size(); j++)
		{
			int t = startNode.NWNext[j];
			prevIds[t].push_back(startNode.ID);
		}

		NWBody[ID1] = startNode;


		

		while (links.size() != 0)
		{
			vector<int> new_links;
			for (int i = 0; i < links.size(); i++)
			{
				NWNode tmpNode = NWNode(Body[links[i]], new_links, Body, ID1, pos1, ID2, pos2, read.size());
				if (!tmpNode.status) continue;
				for (int j = 0; j < tmpNode.NWNext.size(); j++)
				{
					int t = tmpNode.NWNext[j];
					prevIds[t].push_back(tmpNode.ID);
				}

				NWBody[tmpNode.ID] = tmpNode;
			}
			sort(new_links.begin(), new_links.end());
			new_links.erase(unique(new_links.begin(), new_links.end()), new_links.end());
			links = new_links;
		}

		int altsumlen = 0;
		for (auto p1 : NWBody)
		{
			altsumlen = altsumlen + p1.second.ref.size();
		}

		if (altsumlen > read.size() * 3)
		{
			return;
		}

		links = { ID1 };
		startNode.InitMatrix(read);
		while (links.size() != 0)
		{
			vector<int> new_links;
			for (int i = 0; i < links.size(); i++)
			{
				map<int, vector<int>> prevcolmns;
				for (int j = 0; j < prevIds[links[i]].size(); j++)
				{

					prevcolmns[prevIds[links[i]][j]] = NWBody[prevIds[links[i]][j]].NWMatrix.back();
				}

				NWBody[links[i]].InitMatrix(read, prevcolmns);
			} 
			sort(new_links.begin(), new_links.end());
			new_links.erase(unique(new_links.begin(), new_links.end()), new_links.end());
			links = new_links;
		}

		nwID1 = ID1;
		nwID2 = ID2;
		nwpos1 = pos1;
		nwpos2 = pos2;
		loaded = true;
	}

	NWpart NWTrace(string &read)
	{
		NWpart res = NWBody[lastID].NWTrace(read, NWpart());
		while (res.finished != true)
		{
			res = NWBody[res.ID].NWTrace(read, res);
		}
		res.ID1 = nwID1;
		res.ID2 = nwID2;
		res.pos1 = nwpos1;
		res.pos2 = nwpos2;
		return res;
	}

};