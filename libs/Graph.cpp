#include <Graph.h>

#include <chrono>
#include <ctime>

using namespace std;



void Graph::LoadReference(string ref)
{
	Node startNode(ref, 0);
	startNode.B.insert(pair<int, int>(0, -1));
	Body.insert(pair<int, Node>(0, startNode));
	Node endNode("x", -1);
	endNode.End = true;

	Body.insert(pair<int, Node>(-1, endNode));

	Body[-1].End = true;

	Link(&Body[0], &Body[-1]);
}

pair<int, int> Graph::defNodePos(int ThreadId, int pos)
{

	Node curr_Node = Body[0];
	int curr_len = curr_Node.str.length();
	int curr_pos = 0;

	int f_ID = 0;

	while (pos>curr_len)
	{
	//	curr_Node.Next[ThreadId]->print();
		//	system("pause");
		int tmp = curr_Node.Next[ThreadId]->ID;
		curr_Node = Body[tmp];
		curr_len += curr_Node.str.length();
	}

	return pair<int, int>(curr_Node.ID, curr_Node.str.length() - 1 - (curr_len - pos));
}

pair<int, int> Graph::defNodePos(int ThreadId, int pos, Node curr_Node1, int curr_pos1)
{
	int curr_ID = curr_Node1.ID;
	int curr_len = Body[curr_ID].str.length();
	int curr_pos = curr_pos1;

	do
	{
		if ((curr_pos<pos) && (pos <= curr_pos + curr_len))
		{
			return pair<int, int>(curr_ID, (pos - curr_pos));
		}
		else
		{
			curr_ID = Body[curr_ID].B[ThreadId];
			curr_pos += curr_len;
			curr_len = Body[curr_ID].str.length();
		}

	} while (!Body[curr_ID].End);

	return pair<int, int>(-1, -1);
}

void Graph::LoadVariance(int ThreadID1, int ThreadID2, vector<tuple<int, int, string>> VCF)
{
	cout << "Start variance loading" << endl;
	//walk down the thread1, incepting variance node from position list.
	//ints in VCF should be sorted

	int curr_ID = 0;
	int curr_len = Body[curr_ID].str.length();
	int curr_pos = 0;

	int vcf_id = 0;

	bool check = false;
	int ID1 = 0;
	int ID2 = 0;
	int max = VCF.size();
	do
	{
		int pos1 = get<0>(VCF[vcf_id]);
		int pos2 = get<1>(VCF[vcf_id]);
		string str = get<2>(VCF[vcf_id]);
		cout << vcf_id << '/' << max << '\n';
		
		//chrono::time_point<//chrono::system_clock> start1, end1;
		//start1 = //chrono::system_clock::now();
		if ((curr_pos<pos1) && (pos1 <= curr_pos + curr_len))
		{
			//chrono::time_point<//chrono::system_clock> start2, end2;
		//	start2 = //chrono::system_clock::now();
			if (pos1 == curr_pos + curr_len) //начало вариации со стыка нодов
			{
				//chrono::time_point<//chrono::system_clock> start3, end3;
		//		start3 = //chrono::system_clock::now();
				ID1 = Body[curr_ID].ID;
				//end3 = //chrono::system_clock::now();
		//		int elsec3 = //chrono::duration_cast<std:://chrono::milliseconds>(end3 - start3).count();
			//	cout << "Intro1a: " << elsec3 << '\t';
			}
			else
			{
				
				ID1 = Body[curr_ID].ID;
				Node tmpNode1(Body.size());
				Body.insert(pair<int, Node>(tmpNode1.ID, tmpNode1));
				//chrono::time_point<//chrono::system_clock> start3, end3;
			//	start3 = //chrono::system_clock::now();
				Body[curr_ID].Split(pos1 - curr_pos, &Body[tmpNode1.ID]);
		//		end3 = //chrono::system_clock::now();
			//	int elsec3 = //chrono::duration_cast<std:://chrono::milliseconds>(end3 - start3).count();
			//	cout << "Intro1a: " << elsec3 << '\t';
			}
			//end2 = //chrono::system_clock::now();
		//	int elsec = //chrono::duration_cast<std:://chrono::milliseconds>(end2 - start2).count();
		//	cout << "Intro1: " << elsec << '\t';

		//	start2 = //chrono::system_clock::now();

			pair<int, int> tmp_pair1 = defNodePos(ThreadID1, pos2, Body[curr_ID], curr_pos);

		//	end2 = //chrono::system_clock::now();
		//	elsec = //chrono::duration_cast<std:://chrono::milliseconds>(end2 - start2).count();
		//	cout << "Intro2: " << elsec << '\t';

		//	start2 = //chrono::system_clock::now();
			if (Body[tmp_pair1.first].str.length() == tmp_pair1.second) //конец вариации на стык нодов
			{
				Node tmp = Body[tmp_pair1.first];
				ID2 = Body[tmp_pair1.first].B[ThreadID1];
			}
			else
			{
				Node tmpNode1(Body.size());
				ID2 = tmpNode1.ID;
				Body.insert(pair<int, Node>(tmpNode1.ID, tmpNode1));
				Body[tmp_pair1.first].Split(tmp_pair1.second, &Body[tmpNode1.ID]);
			}

		//	end2 = //chrono::system_clock::now();
	//		elsec = //chrono::duration_cast<std:://chrono::milliseconds>(end2 - start2).count();
	//		cout << "Intro3: " << elsec << '\t';

	//		start2 = //chrono::system_clock::now();
			//	RecurcivePrint();
			Node incNode(str, Body.size());
			int incID = incNode.ID;
			
			
			Body.insert(pair<int, Node>(incID, incNode));
			Link(&Body[ID1], &Body[incID]);
			Link(&Body[incID], &Body[ID2]);
			Body[incID].B.insert(pair<int, int>(ThreadID2, ID2));
			Body[ID1].B.insert(pair<int, int>(ThreadID2, incID));

			check = false;

			if (vcf_id<VCF.size() - 1)
			{
				vcf_id++;
				
			}
		//	end2 = //chrono::system_clock::now();
		//	elsec = //chrono::duration_cast<std:://chrono::milliseconds>(end2 - start2).count();
		//	cout << "Intro4: " << elsec << '\r';

		}
		else
		{
			if (check)
			{
				Body[curr_ID].B.insert(pair<int, int>(ThreadID2, Body[curr_ID].B[ThreadID1]));
			}
		}
	//	end1 = //chrono::system_clock::now();
	//	int elapsed_seconds = //chrono::duration_cast<std:://chrono::milliseconds>(end1 - start1).count();
	//	cout << "Total execution: " << elapsed_seconds << '\r';
			
		curr_pos += Body[curr_ID].str.length();
		curr_ID = Body[curr_ID].B[ThreadID1];
		curr_len = Body[curr_ID].str.length();
		check = (curr_ID == ID2);
	} while (Body[curr_ID].End != true);
}



void Graph::Save(string filename)
{
	ofstream f_out;
	f_out.open(filename+".gngg");
	for (auto n : Body)
	{

		f_out << n.second.ID << '.' << n.second.str << '.' << n.second.NaNbp << '.' << n.second.End << '.' << n.second.Next.size();
		for (auto n1 : n.second.Next)
		{
			f_out << '.' << n1.second->ID;
		}

		f_out << '.' << n.second.Prev.size();
		for (auto n1 : n.second.Prev)
		{
			f_out << '.' << n1.second->ID;
		}
		f_out << '.' << n.second.B.size();
		for (auto n1 : n.second.B)
		{
			f_out << '.'<< n1.first << '/' << n1.second;
		}
		f_out << '.' << endl;
	}
	f_out.close();
}


void Graph::Load(string filename)
{
	ifstream f_in;
	f_in.open(filename);
	do
	{
		LoadNodep1(f_in);
	}
	while (!f_in.eof());
	f_in.close();
	f_in.open(filename);
	do
	{
		LoadNodep2(f_in);
	} while (!f_in.eof());
	f_in.close();
}

void Graph::LoadNodep1(ifstream &f_in)
{
	string buf;
	getline(f_in, buf, '.');
	if (buf == "") return;
	int ID = stoi(buf);
	getline(f_in, buf, '.');
	string str = buf;
	getline(f_in, buf, '.');
	bool NaNbp = (buf == "1");
	getline(f_in, buf, '.');
	bool End = (buf == "1");
	
	Node nn(str, ID);
	nn.End = End;
	
	Body[nn.ID] = nn;

	getline(f_in, buf);
	return;
}

void Graph::LoadNodep2(ifstream &f_in)
{
	string buf;
	getline(f_in, buf, '.');
	if (buf == "") return;
	int ID = stoi(buf);
	getline(f_in, buf, '.');
	string str = buf;
	getline(f_in, buf, '.');
	bool NaNbp = (buf == "1");
	getline(f_in, buf, '.');
	bool End = (buf == "1");

	getline(f_in, buf, '.');
	int sz = stoi(buf);  //NEXT
	for (int i = 0; i < sz; i++)
	{
		getline(f_in, buf, '.');
		int lnk = stoi(buf);
		Body[ID].Next[lnk] = &(Body[lnk]);
	}
	getline(f_in, buf, '.');
	sz = stoi(buf); //PREV
	for (int i = 0; i < sz; i++)
	{
		getline(f_in, buf, '.');
		int lnk = stoi(buf);
		Body[ID].Prev[lnk] = &(Body[lnk]);
	}
	getline(f_in, buf, '.');
	sz = stoi(buf); //B
	for (int i = 0; i < sz; i++)
	{
		getline(f_in, buf, '/');
		int t1 = stoi(buf);
		getline(f_in, buf, '.');
		int n1 = stoi(buf);
		Body[ID].B.insert(pair<int, int>(t1, n1));
	}
	getline(f_in, buf);
}

void Graph::BuildIndex(int base, int len)
{
	GH.Init(base, len);
	for (auto t1 : Body)
	{
		if ((t1.second.ID != (-1)) && (!t1.second.NaNbp))
		{
			GenPrimeTable(t1.second, base, len);
		}
		else
		{
			cout << "try to hash -1 or NaNbp Node\n";
		}
	}
}


void Graph::GenPrimeTable(Node n1, int base, int len)
{

	GH.ReInit(n1.ID, n1.str);
	int pos = 0;
	unsigned long long res = 0;

	for (int i = 0; i<n1.str.length(); i++)
	{
		if (GH.Next(n1.str[i]) == 0)
		{
			AddHash(n1.ID, GH.pos, n1.ID, i, GH.currhash);
		}
		else
		{
			//	cout<<"oooops"<<endl;
		}
	}


	for (auto node : n1.Next)
	{
		CallNode(*(node.second), GH);
	}
}

void Graph::CallNode(Node n1, GenerHash GH2)
{
	if (n1.ID == -1) return;
	if (n1.NaNbp) return;
	unsigned long long res = 0;
	for (int i = 0; i<n1.str.length(); i++)
	{
		res = GH2.Next(n1.str[i]);
		if (res == -2)//pointer to position at the end of node. return;
		{
			return;
		}
		else
		{
			AddHash(GH2.id, GH2.pos, n1.ID, i, GH2.currhash);
		}
	}

	for (auto p1 : n1.Next)
	{
		CallNode(*(p1.second), GH2);
	}
}

void Graph::AddHash(int ID, int pos, int fID, int fpos, unsigned long long hash)
{
	WArray tmp;
	tmp.value[0] = ID;
	tmp.value[1] = pos;
	tmp.value[2] = fID;
	tmp.value[3] = fpos;

	try
	{
		hashtable[hash].push_back(tmp);
	}
	catch (...)
	{
		cout << "oooooops, exeption" << endl;
	}
}


void Graph::SaveIndex(string filename)
{
	ofstream f_out;
	f_out.open(filename);
	f_out << GH.base << '.' << GH.len << endl;
	for (auto t : hashtable)
	{
		f_out << t.first << '.' << t.second.size();
		for (int i = 0; i < t.second.size(); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				f_out << '.' << t.second[i].value[j];
			}
		}
		f_out << '.' << endl;
	}
	f_out.close();
}


void Graph::LoadIndex(string filename)
{
	hashtable.clear();
	ifstream f_in;
	f_in.open(filename);

	string buf;
	getline(f_in, buf, '.');
	int ghbase = stoi(buf);
	getline(f_in, buf);
	int ghlen = stoi(buf);
	GH.Init(ghbase, ghlen);
	do
	{
		getline(f_in, buf, '.');
		if (buf == "") return;
		unsigned long long h1 = stoull(buf);
		getline(f_in, buf, '.');
		int sz = stoi(buf);
		WArray tmp1;
		vector<WArray> wr_tmp;
		hashtable[h1] = wr_tmp;
		for (int i = 0; i < sz; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				getline(f_in, buf, '.');
				tmp1.value[j] = stoi(buf);
			}
			hashtable[h1].push_back(tmp1);
		}
		getline(f_in, buf);
	} while (!f_in.eof());
}


void Graph::Clear()
{
	hashtable.clear();
	Body.clear();
}
//==========================BubbleIndexBuild=======================================================
void Graph::BubbleIndexBuild()
{
	vector<int> markedNode;

	Body[0].bidx = BubbleIndex(true);
	markedNode.push_back(0);
	
	for (int i = 0; i<markedNode.size(); i++)
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

void Graph::BubbleIndexBuildTest()
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


//================================================================================================

int Graph::NotSimpleFinderSNAP(string str, int err_count)
{
	cout << "SNAPrun->" << endl;// << str << endl;
	if (str.length()<GH.len) return -1;
	int res = -1;
	vector<unsigned long long> hash_seeds;
	GH.ReInit(-2, str);
	vector<int> positions;

	for (int i = 0; i<str.length(); i++)
	{
		int res = GH.Next(str[i]);
		if (res == 0)
		{
			hash_seeds.push_back(GH.currhash);
		}
	}



	vector<vector<WArray>> hit_place; //(hash_seeds.size());

	for (int i = 0; i<hash_seeds.size(); i++)
	{
		vector<WArray> tmp;
		hit_place.push_back(tmp);
		for (int j = 0; j<hashtable[hash_seeds[i]].size(); j++)
		{
			hit_place[i].push_back(hashtable[hash_seeds[i]][j]);
		}
	}


	res = ASTRun(str, str.length(), hit_place, err_count);
	return res;
}

int Graph::ASTRun(string str, int pos_end, vector<vector<WArray>> hit_place, int err_count)
{
	cout << "Astrun->" << endl;// str << endl;

	bool mishit;

	int z_count = 0;
	int pos = 0;
	int pos_check = GH.len - 1;
	int step = GH.len - 1;

	int htsz = hit_place.size();
	cout << "htsz: " << htsz << endl;
	if (htsz <= step)
	{
		cout << "ASTRun. Short" << endl;
		return 0;
	}

	vector<vector<WArray>> VP1;
	vector<WArray> VP2;
	vector<int> positions;

	VP1.push_back(hit_place[0]);

	cout << "forward_run" << endl;
	do
	{
		pos += step;
		positions.push_back(pos);
		//cout << "pos: " << pos << endl;
		//cout << "Running dangerous function" << endl;
		VP2 = NextPos(VP1.back(), hit_place[pos]);
		//cout << "We are here! We are returned" << endl;
		if (VP2.size() == 0)
		{
			if (z_count >= err_count)
				return -1;
			pos -= step;

			//cout << "err_block, pos: " << pos << endl;
			step = GH.len * 2 + 1;

			int tmp_hps = hit_place.size();

			if (pos + step >= tmp_hps)
			{
				//cout << "err_block, tail of read" << endl;
				step = tmp_hps - pos - 1;
			}

			//cout << "err_block, NJ" << endl;
			NodeJump(VP1.back(), step, VP2); //обычный джамп

			z_count++;

			vector<WArray> tmpVP3;

			//cout << "err_block, NP: " << VP2.size() << '/' << hit_place[pos + step].size() << "\t" << pos + step << "\t" << hit_place.size() << endl;
			//cout << "second entry into dangerous function" << endl;
			tmpVP3 = NextPos(VP2, hit_place[pos + step]); //check mis hit (sorry, but I read "check my shit" all times)

														  //cout << "err_block, NP_end:"; //<< tmpVP3.size();
			if (tmpVP3.size() == 0)
			{
				//cout << "err_block, NW+" << endl;
				step += 1;
				if (pos + step<tmp_hps)
					tmpVP3 = NextPos(VP2, hit_place[pos + step]); //checkdeletion
			}
			if (tmpVP3.size() == 0)
			{
				////cout << "err_block, NW-" << endl;
				step -= 2;
				if (pos + step<tmp_hps)
					tmpVP3 = NextPos(VP2, hit_place[pos + step]); //checkinseption
			}

			if (tmpVP3.size() == 0) return -1;

			VP2 = NextPos_back(VP2, tmpVP3);

			if (VP2.size() == 0) return -1;

			//VP2 = NW_check(VP2, str.substr(pos, step));

			if (VP2.size() == 0) return -1;

			VP1.push_back(VP2);
			VP1.push_back(tmpVP3);
			pos += step;
			pos_check += step;
		}
		else
		{

			VP1.push_back(VP2);

		}
		pos_check += step;
		step = GH.len - 1;

	} while (pos_check <= pos_end - GH.len);


	if (pos != hit_place.size() - 1)
	{
		//VP2 = TailCheck(VP1.back(), hit_place.back(), hit_place.size() - 1 - pos);
	}
	VP1.pop_back();
	VP1.push_back(VP2);


	//backward run
	for (int i = VP1.size() - 1; i >= 0; i--)
	{

	}





	cout << "ASTRun completed" << endl;

	for (int i = 0; i < VP1.size(); i++)
	{
		for (int j = 0; j < VP1[i].size(); j++)
		{
			cout << VP1[i][j].value[0] << '.' << VP1[i][j].value[1] << '\t';
		}
		//cout << '\n';
	}


	return 0;
}

vector<WArray> Graph::NextPos(vector<WArray>& VP1_last, vector<WArray> hits)
{
	vector<WArray> res;

	for (auto v : VP1_last)
	{
		int node = v.value[2];
		int pos = v.value[3];
		for (int i = 0; i<hits.size(); i++)
		{
			//cout << "NP: " << node << '.' << pos << '.' << hits[i].value[0] << '.' << hits[i].value[1] << endl;
			if ((node == hits[i].value[0]) && (pos == hits[i].value[1]))
			{
				//cout << "NP: checked" << endl;
				WArray *tmp = new WArray();
				tmp->value[0] = hits[i].value[0];
				tmp->value[1] = hits[i].value[1];
				tmp->value[2] = hits[i].value[2];
				tmp->value[3] = hits[i].value[3];
				//cout << "NP: tmp created" << endl;
				res.push_back(*tmp);
				//cout << "pushed" << endl;
			}
		}
	}
	//cout << "Prepared for sort" << endl;
	sort(res.begin(), res.end());
	//cout << "Sorted" << endl;
	res.resize(abs(distance(unique(res.begin(), res.end()), res.begin())));
	//cout << "Resized" << endl;
	return res;


}

vector<WArray> Graph::NextPos_back(vector<WArray> VP1_last, vector<WArray> hits)
{
	vector<WArray> res;
	for (auto v : VP1_last)
	{
		int node = v.value[2];
		int pos = v.value[3];
		for (int i = 0; i<hits.size(); i++)
		{
			if ((node == hits[i].value[0]) && (pos == hits[i].value[1]))
			{
				res.push_back(v);
			}
		}
	}

	sort(res.begin(), res.end());
	res.resize(unique(res.begin(), res.end()) - res.begin());
	return res;


}

void Graph::NodeJump(vector<WArray> VP1, int len, vector<WArray> &res, int c_pos)
{
	for (auto v : VP1)
	{
		int node = v.value[0];
		int pos = v.value[1];
		NodeJump(Body[node], pos, node, pos, len, res, c_pos);
	}

	sort(res.begin(), res.end());
	res.resize(unique(res.begin(), res.end()) - res.begin());

}

void Graph::NodeJump(Node n, int pos, int s_n, int s_pos, int len, vector<WArray> &res, int c_pos)
{
	if ((n.str.length() - pos)>(len - c_pos))
	{
		//cout << "NJ_push" << s_n << '.' << s_pos << '.' << n.ID << '.' << pos + len - c_pos << endl;

		WArray tmp;
		tmp.value[0] = s_n;
		tmp.value[1] = s_pos;
		tmp.value[2] = n.ID;
		tmp.value[3] = pos + len - c_pos;
		res.push_back(tmp);
		//cout << "return here" << endl;
		return;
	}
	//cout << "final node jump first block" << endl;
	for (auto n1 : n.Next)
	{
		NodeJump(*n1.second, 0, s_n, s_pos, len, res, c_pos + n.str.length() - pos);
	}
}

vector<WArray> Graph::TailCheck(vector<WArray> VP1, vector<WArray> end, int len)
{
	vector<WArray> res;
	for (auto v : VP1)
	{
		int node = v.value[2];
		int pos = v.value[3];
		int c_pos = 0;
		TailCheck(Body[node], pos, node, pos, len, end, res);
	}
	return res;
}

void Graph::TailCheck(Node n, int pos, int s_n, int s_pos, int len, vector<WArray> end, vector<WArray> &res, int c_pos)
{
	for (int i = 0; i < end.size(); i++)
	{
		if ((n.ID == end[i].value[0]) && (len - c_pos))
		{
			WArray tmp;
			tmp.value[0] = s_n;
			tmp.value[1] = s_pos;
			tmp.value[2] = n.ID;
			tmp.value[3] = pos + len - c_pos;
			res.push_back(tmp);
			return;
		}
	}
	for (auto n1 : n.Next)
	{
		TailCheck(*n1.second, 0, s_n, s_pos, len, end, res, c_pos + n.str.length() - pos);
	}
}



//============SFinder=========================================================================================
int Graph::SFinder(string str, float cov)
 {
	chrono::time_point<chrono::system_clock> start1, end1;
	int elsec;
	start1 = chrono::system_clock::now();
	if (str.length()<GH.len) return -1;
	int res = -1;

	int countN = 0;
	for (int i = 0; i < str.size(); i++)
	{
		if ((str[i] == 'N') || (str[i] == 'n'))
		{
			countN++;
		}
	}
	if (countN > 5) return -1;
	
	vector<vector<unsigned long long>> hash_seeds;

	StringVectorHash SVH = StringVectorHash(GH.base, GH.len);
	hash_seeds = SVH.GetHash(str);
	
	end1 = chrono::system_clock::now();
	elsec = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
	cout << "Hash generating: " << elsec << '\n';
	
	start1 = chrono::system_clock::now();
	vector<vector<WArray>> hit_place;
	for (int i = 0; i<hash_seeds.size(); i++)
	{
		vector<WArray> tmp;
		hit_place.push_back(tmp);
		for (int k = 0; k < hash_seeds[i].size(); k++)
		{
			for (int j = 0; j < hashtable[hash_seeds[i][k]].size(); j++)
			{
				hit_place[i].push_back(hashtable[hash_seeds[i][k]][j]);
			}
		}
	}
	end1 = chrono::system_clock::now();
	elsec = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
	cout << "Hit placing: " << elsec << '\n';
	
	start1 = chrono::system_clock::now();


	SAligner saligner(GH.len, str, hit_place, cov);

	end1 = chrono::system_clock::now();
	elsec = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
	cout << "SAligner build: " << elsec << '\n';

	start1 = chrono::system_clock::now();

	int i = saligner.Run(Body);

	end1 = chrono::system_clock::now();
	elsec = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
	cout << "SAligner run: " << elsec << '\n';

	start1 = chrono::system_clock::now();

     res = saligner.Ligation(str.length(), Body);


	 end1 = chrono::system_clock::now();
	 elsec = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
	 cout << "SAligner ligation: " << elsec << '\n';
	return res;
}

//===========NWPreTracker======================================================================================

int Graph::GetRelation(pair<int,int> p1, pair<int, int> p2)
{
	 
	int q = Body[p1.first].GetRelation(Body[p2.first]); //-1 = prev, 0 = parallel, 1 = next;
	if (q == 3)
	{
		if (p1.second < p2.second)
		{
			q = 1;
		}
		else if (p1.second > p2.second)
		{
			q = -1;
		}
		else
		{
			q = 0;
		}
	}

	return q;
}

void Graph::printThread(int ID, int ThreadID, string output, string ref)
{
	string res = "";
	ofstream f_out;
	f_out.open(output);
	int curr_pos = 0;
	f_out << "seq_ref" + to_string(ThreadID) << '\n';
	int currID = ID;
	do
	{

		Node tmp = Body[currID];
		res = res + Body[currID].str;
		curr_pos = curr_pos + tmp.str.length();
		string t_ref = ref.substr(1, curr_pos);
		if (res != t_ref)
		{
			cout << "trackmatching error." << endl;
			cout << "ref.length: " << t_ref.length() << endl;
			cout << "track.length: " << res.length() << endl;
			return;
		}
		currID = Body[currID].B[ThreadID];
	} while (!Body[currID].End);

	f_out << '\n';
	f_out.close();


}

