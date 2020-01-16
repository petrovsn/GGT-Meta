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
#include <thread>
#include <future>


#include <ctime>
#include <cstdio>

#include <GenerHash.h>
#include <Fasta.h>
#include <Node.h>
#include <WArray.h>
#include<NWAligner.h>
//#include <NWAligner.h>

using namespace std;

class FPoint
{
public:
	int score;
	int coverage;
	vector<int> seedIDs;
	vector<NWpart> NWlinks;

	FPoint()
	{
		score = 0;
		coverage = 0;
	}
};



class  FAlignment
{
public:
	vector<TArray> seeds;
	vector<NWpart> NWlinks;
	vector<GVariation> vars;
	int score;
	int coverage;

	//NWAligner nwaln;

	FAlignment(FPoint &fmap, vector<TArray>& StartTracks, string read, map<int, Node> &Body)
	{
		for (int i = 0; i < fmap.seedIDs.size(); i++)
		{
			seeds.push_back(StartTracks[fmap.seedIDs[i]]);
		}

		for (int i = 0; i < fmap.NWlinks.size(); i++)
		{
			for (int j = 0; j < fmap.NWlinks[i].vars.size(); j++)
			{
				vars.push_back(fmap.NWlinks[i].vars[j]);
			}
		}

		NWlinks = fmap.NWlinks;
		score = fmap.score;
		coverage = fmap.coverage;
		
	}


};

class Graph
{
public:
	map<int, Node> Body;
	map<unsigned long long, vector<WArray>> Hashtable;
	
	int hashBase = -1;
	int hashLen = -1;
	
	
	Graph(){}

	//node operations

	void Link(Node* N1, Node* N2)
	{
		(N1->Next).insert(pair<int, Node*>(N2->ID, N2));
		(N2->Prev).insert(pair<int, Node*>(N1->ID, N1));
	}

	void LoadReference(string ref)
	{
		Node startNode(ref, 0);
		startNode.B.insert(pair<int, int>(0, -1));
		Body.insert(pair<int, Node>(0, startNode));

		Node endNode("x", -1);
		endNode.End = true;
		Body.insert(pair<int, Node>(-1, endNode));

		Link(&Body[0], &Body[-1]);

	}

	

	//all for hash generator

	void BuildIndex(int base, int len)
	{
		Hashtable.clear();
		hashBase = base;
		hashLen = len;
		for (auto t1 : Body)
		{
			if ((t1.second.ID != (-1)) && (!t1.second.NaNbp))
			{
				GenHashIndexForNode(t1.second);
			}
			else
			{
				cout << "try to hash -1 or NaNbp Node\n";
			}
		}
		BubbleIndexBuild();
	}

	void GenHashIndexForNode(Node n1)
	{

		GenerHash GH = GenerHash(hashBase, hashLen);
		GH.ReInit(n1.ID, n1.str);

		int res = 0;
		
		WArray tmpWarray = WArray(n1.ID, 0);

		for (int i = 0; i < n1.str.length(); i++)
		{
			if (GH.Next(n1.str[i]) == 0)
			{	
				tmpWarray.finish(i);
				AddHash(tmpWarray, GH.currhash);


				tmpWarray.move_right();
			}
			else
			{
				//	cout<<"oooops"<<endl;
			}
		}

		for (auto node : n1.Next)
		{
			CallNode(*(node.second), GH, tmpWarray);
		}
	}

	void CallNode(Node n1, GenerHash GH2, WArray wray)
	{
		if (n1.ID == -1) return;
		if (n1.NaNbp) return;
		unsigned long long res = 0;

		wray.extend(n1.ID);
		for (int i = 0; i < n1.str.length(); i++)
		{
			res = GH2.Next(n1.str[i]);
			if (res == -3)//pointer to position at the end of initial node. return;
			{
				return;
			}
			else
			{
				wray.finish(i);
				AddHash(wray, GH2.currhash);

				wray.move_right();
				
			}
		}

		for (auto p1 : n1.Next)
		{
			CallNode(*(p1.second), GH2, wray);
		}
	}

	void AddHash(WArray warray, unsigned long long hash)
	{
		Hashtable[hash].push_back(warray);
	}

	

	//B-index build

	void BubbleIndexBuild()
	{
		vector<int> markedNode;

		Body[0].bidx = BubbleIndex(true);
		markedNode.push_back(0);

		for (int i = 0; i < markedNode.size(); i++)
		{

			bool bubble = (Body[markedNode[i]].Next.size() > 1);
			int c = 0;
			for (auto n : Body[markedNode[i]].Next)
			{
				Body[n.first].GenBIdx(Body[markedNode[i]].bidx, c, bubble);
				markedNode.push_back(n.first);
				c++;
			}
		}
	}

	void BubbleIndexBuildTest()
	{
		int Bs = Body.size();
		Body[0].bidx = BubbleIndex(true);
		int c = 0;
		for (auto t : Body[0].Next)
		{
			t.second->GenBIdx(Body[0].bidx, c, true);
			c++;
		}

		for (int i = 2; i < Bs; i++)
		{
			bool bubble = (Body[i].Next.size() > 1);
			c = 0;
			for (auto n : Body[i].Next)
			{
				n.second->GenBIdx(Body[i].bidx, c, bubble);
				c++;
			}
		}
	}

	//alignment&SNPcalling


	vector<TArray> GenerateSeeds(string& read)
	{
		GenerHash GH2 = GenerHash(hashBase, hashLen);
		GH2.ReInit(-2, read);

		vector<unsigned long long> num_read;
		for (int i = 0; i < read.size(); i++)
		{
			int res_hash = GH2.Next(read[i]);
			if (res_hash == 0)
			{
				num_read.push_back(GH2.currhash);
			}
		}

		vector<vector<WArray>> marks_read;
		for (int i = 0; i < num_read.size(); i++)
		{
			marks_read.push_back(Hashtable[num_read[i]]);
		}

		vector<TArray> StartTracks;
		for (int i = 0; i < marks_read[0].size(); i++)
		{
			StartTracks.push_back(TArray(marks_read[0][i], 0, hashLen - 1));
		}

		for (int i = 1; i < marks_read.size(); i += 1)
		{
			for (int j = 0; j < StartTracks.size(); j++)
			{
				for (int k = 0; k < marks_read[i].size(); k++)
				{
					int res = StartTracks[j].TryExtRight(marks_read[i][k], Body, i, i + hashLen - 1);
					if (res == 0)
					{
						break;
					}
				}
			}

			for (int k = 0; k < marks_read[i].size(); k++)
			{
				if (marks_read[i][k].assembled == false)
				{
					StartTracks.push_back(TArray(marks_read[i][k], i, i + hashLen));
				}
			}
		}

		return StartTracks;
	}

	void TrimSeeds(vector<TArray>& StartTracks)
	{
		for (int i = 0; i < StartTracks.size(); i++)
		{
			for (int j = 0; j < StartTracks.size(); j++)
			{
				if (i != j)
				{
					bool comp = StartTracks[i].IntersecSeed(StartTracks[j], Body);
					if (comp)
					{
						//cout << "GetIntersec\n";
						StartTracks[j].CutFromLeft(StartTracks[i], Body);
					}
				}
			}
		}
	}

	NWAligner nwa;

	FPoint FillLine(vector<vector<int>> &WeightMatrix, vector<TArray> &StartTracks, int startline, string &read)
	{
		vector<int> starts = { startline };
		vector<FPoint> node_weights = vector<FPoint>(WeightMatrix.size());


		node_weights[startline].seedIDs.push_back(startline);
		node_weights[startline].score = StartTracks[startline].read_endPos - StartTracks[startline].read_startPos;
		node_weights[startline].coverage = node_weights[startline].score;


		while (starts.size() != 0)
		{
			vector<int> new_starts;
			for (int idx = 0; idx < starts.size(); idx++)//for each of the edges
			{
				int i = starts[idx];
				for (int j = 0; j < WeightMatrix[i].size(); j++)
				{
					
					if (WeightMatrix[i][j] > 0)
					{
						int dist = StartTracks[j].read_startPos - StartTracks[i].read_endPos;

						int ID1 = StartTracks[i].NodeIDs.back();
						int pos1 = StartTracks[i].endPos;
						int ID2 = StartTracks[j].NodeIDs[0];
						int pos2 = StartTracks[j].startPos;
						int rcoord = StartTracks[i].read_endPos;
						int rlen = StartTracks[j].read_startPos - StartTracks[i].read_endPos + 1;
						string nwread = read.substr(rcoord, rlen);

					
						nwa.Init(ID1, pos1, ID2, pos2, Body, nwread);

						int score = 0;
						int coverage = 0;

						NWpart nw_res;

						if (nwa.loaded)
						{
							nw_res = nwa.NWTrace(nwread);
							score = score + nw_res.score-2;
							coverage = coverage + nw_res.score-2;

						}
						else
						{
							score = score - 1000;
							coverage = coverage - 1000;
						}

						int RelScore = WeightMatrix[i][j] + coverage + node_weights[i].score;

						if (RelScore >= node_weights[j].score)
						{
							node_weights[j].seedIDs = node_weights[i].seedIDs;
							node_weights[j].seedIDs.push_back(j);
							node_weights[j].score = RelScore;
							node_weights[j].coverage = node_weights[i].coverage + WeightMatrix[i][j];
							node_weights[j].NWlinks = node_weights[i].NWlinks;
							node_weights[j].NWlinks.push_back(nw_res);
							new_starts.push_back(j);
						}
					}
				}
			}
			starts = new_starts;
		}
		
		int fcov = 0; 
		FPoint res; 
		for (int i = 0; i < node_weights.size(); i++)
		{
			if (node_weights[i].coverage > fcov)
			{
				fcov = node_weights[i].coverage;
				res = node_weights[i];
			}
		}
		return res;
	}

	vector<FAlignment> AlignHashSmWtmn(string read, float coverage_treshhold = 0.9)
	{
		
		vector<TArray> StartTracks = GenerateSeeds(read); //call uninterrupted seeds

		TrimSeeds(StartTracks); //trim seeds for indels with repeat;
		

		//plot adjance matrix
		vector<vector<int>> WeightMatrix = vector<vector<int>>(StartTracks.size());

		for (int i = 0; i < WeightMatrix.size(); i++)
		{
			WeightMatrix[i] = vector<int>(StartTracks.size());
			for (int j = 0; j < WeightMatrix[i].size(); j++)
			{
				WeightMatrix[i][j] = 0;
			}
		}

		for (int i = 0; i < StartTracks.size(); i++)
		{
			for (int j = 0; j < StartTracks.size(); j++)
			{
				if (i != j)
				{
					if (StartTracks[i].read_endPos < StartTracks[j].read_startPos)
					{
						int ID1 = StartTracks[i].NodeIDs.back();
						int pos1 = StartTracks[i].endPos;
						int ID2 = StartTracks[j].NodeIDs[0];
						int pos2 = StartTracks[j].startPos;

						bool comp = Node1LessNode2(ID1, pos1, ID2, pos2, Body);
						
						if (comp)
						{
							int w_plus = StartTracks[j].read_endPos - StartTracks[j].read_startPos +1;
							WeightMatrix[i][j] = w_plus;
						}
					}
				}
			}
		}

		vector<int> idxs_start; //all possible start points of D-Run;
		
		for (int j = 0; j < WeightMatrix.size(); j++)
		{
			int r = 0;
			for (int i = 0; i < WeightMatrix.size(); i++)
			{
				r = r + WeightMatrix[i][j];
			}
			if (r == 0)
			{
				idxs_start.push_back(j);
			}
		}

	
		vector<FPoint> maps; 
		for (int i = 0; i < idxs_start.size(); i++)
		{
			FPoint res = FillLine(WeightMatrix, StartTracks, idxs_start[i], read);
			maps.push_back(res);
		}

		vector<FAlignment> total_aln;
		vector<FAlignment> filtered_aln;

		int maxscore = -1;
		int maxmapid = -1;
		for (int i = 0; i < maps.size(); i++)
		{
			if (maps[i].score > maxscore)
			{
				maxscore = maps[i].score;
				maxmapid = i;
			}
		} 
		
		FAlignment tmp = FAlignment(maps[maxmapid], StartTracks, read, Body);

		float rescov = tmp.score * 1.0 / read.size();
		if (rescov >= coverage_treshhold)
		{
			filtered_aln.push_back(tmp);
		}
						
		return filtered_aln;
	}
};




int main(int argc, const char* argv[])
{
	string ref = loadfasta("data/Ecoli_O157.fasta");
		Graph GMAIN = Graph();

	GMAIN.LoadReference(ref);

	GMAIN.BuildIndex(5, 13);
	
	map<string, vector<FAlignment>> results;
	
	//vector<pair<string, string>> reads = loadmultifasta("data/SRR65_reads_ssample.fa");
	//vector<pair<string, string>> reads = loadmultifasta("data/test.fa");
	vector<pair<string, string>> reads = loadmultifasta("data/SRR65_reads.fa"); 
	ofstream fout = ofstream("data/misaligned_total.fa");

 
	clock_t res = 0;

	for (int i = 0; i < reads.size(); i++)
	{
		clock_t time_a = clock();
		vector<FAlignment> m1 = GMAIN.AlignHashSmWtmn(reads[i].second, 0.8);
		if (m1.size() == 0)
		{
			string revread = reverse(reads[i].second);
			m1 = GMAIN.AlignHashSmWtmn(revread, 0.8);
		}
		clock_t time_b = clock() - time_a;
		//cout << time_b << '\r';
		res = res + time_b;
		cout << res *1.0 / (i+1) << '\t' << i<< '\r';
	}
	fout.close();


	
	/*string s1 = "abcdefgjklmn"+ref;
	int len = s1.size();
	//char* s2 = new char(len);;
	char* s2 = (char*)(ref.c_str());

	clock_t time_a = clock();
	string s1_n = s1.substr(1000, 50000);
	clock_t time_b = clock();

	cout  << '\t'<< (time_b-time_a)<< endl;
	time_a = clock();
	char* s2_n = s2 + 4;
	time_b = clock();
	cout  << '\t' << (time_b - time_a) << endl;*/
}


//int main2(int argc, const char* argv[])
int main2()
{
	return 0;
	//string ref = "atcaattccggaaatttcccgggaa";
	string read = "atcttAcgggaaattcccggaAa";

	//nwpart res1 = NWAlign(ref, read);

	//cout << res1.alnref << endl;
	//cout << res1.alnread << endl;
	//cout << res1.lastpos << endl;

	Graph GMAIN = Graph();

	string ref1 = "atcaattccgg";
	
	string ref2 = "aaatttcccgggaa";

	Node n1 = Node("act", 0);
	Node n2 = Node("atcaattccgg", 2);
	Node n3 = Node("atcaattccggg", 3);
	Node n4 = Node("aaatttcccgggaa", 4);
	Node n_end = Node("", -1);


	GMAIN.Body[n1.ID] = n1;
	GMAIN.Body[n2.ID] = n2;
	GMAIN.Body[n3.ID] = n3;
	GMAIN.Body[n4.ID] = n4;
	GMAIN.Body[n_end.ID] = n_end;

	GMAIN.Link(&GMAIN.Body[n1.ID], &GMAIN.Body[n2.ID]);
	GMAIN.Link(&GMAIN.Body[n1.ID], &GMAIN.Body[n3.ID]);

	GMAIN.Link(&GMAIN.Body[n2.ID], &GMAIN.Body[n4.ID]);
	GMAIN.Link(&GMAIN.Body[n3.ID], &GMAIN.Body[n4.ID]);

	GMAIN.Link(&GMAIN.Body[n4.ID], &GMAIN.Body[n_end.ID]);

	GMAIN.BuildIndex(5, 7);
	

	NWAligner nwal = NWAligner();
	nwal.Init(0, 1, 4, 13, GMAIN.Body, read);
	NWpart res = nwal.NWTrace(read);
	cout << res.alnref << endl;
	cout << res.alnread << endl; 

	/*
	r1.InitMatrix(read);
	r1alt.InitMatrix(read);

	map<int, vector<int>> prevcolmns;
	prevcolmns[1] = r1.NWMatrix.back();
	prevcolmns[3] = r1alt.NWMatrix.back();

	r2.InitMatrix(read, prevcolmns);

	printNWmart(r1.ref, read, r1.NWMatrix);
	cout << endl;
	printNWmart(r2.ref, read, r2.NWMatrix);


	map<int, NWNode> nwbody;

	nwbody[1] = r1;
	nwbody[2] = r2;
	nwbody[3] = r1alt;

	NWpart res;

	res = r2.NWTrace(read, res);
	cout << res.alnref << endl;
	cout << res.alnread << endl;

	res = nwbody[res.ID].NWTrace(read, res);
	cout << res.alnref << endl;
	cout << res.alnread << endl;*/

}