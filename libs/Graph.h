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

#include <GenerHash.h>
#include <Node.h>
#include <WArray.h>
#include <TAligner.h>
#include <SAligner.h>

using namespace std;

class Graph
{
public:
	map<int, Node> Body;
	unordered_map<unsigned long long, vector<WArray>> hashtable;
	bool IdxBuilded;

	void Link(Node* N1, Node* N2)
	{
		(N1->Next).insert(pair<int, Node*>(N2->ID, N2));
		(N2->Prev).insert(pair<int, Node*>(N1->ID, N1));
	}
	pair<int, int> defNodePos(int ThreadId, int pos);
	pair<int, int> defNodePos(int ThreadId, int pos, Node curr_Node1, int curr_pos1);
	void LoadNodep1(ifstream &f_in);
	void LoadNodep2(ifstream &f_in);


	GenerHash GH;
	void CallNode(Node n1, GenerHash GH2);
	void GenPrimeTable(Node n1, int base, int len);
	void AddHash(int ID, int pos, int fID, int fpos, unsigned long long hash);

	int ASTRun(string str, int pos_end, vector<vector<WArray>> hit_place, int err_count = 0);
	vector<WArray> NextPos(vector<WArray>& VP1_last, vector<WArray> hits);
	vector<WArray> NextPos_back(vector<WArray> VP1_last, vector<WArray> hits);
	void NodeJump(vector<WArray> VP1, int len, vector<WArray> &res, int c_pos = 0);
	void NodeJump(Node n, int pos, int s_n, int s_pos, int len, vector<WArray> &res, int c_pos = 0);
	vector<WArray> TailCheck(vector<WArray> VP1, vector<WArray> end, int len);
	void TailCheck(Node n, int pos, int s_n, int s_pos, int len, vector<WArray> end, vector<WArray> &res, int c_pos = 0);

	void BubbleIndexBuild();
	void BubbleIndexBuildTest();
	int GetRelation(pair<int, int> p1, pair<int, int> p2);




	Graph() {};
	void LoadReference(string ref);

	
	void LoadVariance(int ThreadID1, int ThreadID2, vector<tuple<int, int, string>> VCF);
	
	void BuildIndex(int base, int len);

	void Save(string filename);
	void Load(string filename);

	void SaveIndex(string filename);
	void LoadIndex(string filename);

	void Clear();

	int NotSimpleFinderSNAP(string str, int err_count = 0);

	int TFinder(string str, int err_count = 0);


	int SFinder(string str, float cov);

	void printThread(int ID, int ThreadID, string output, string ref);
	

};