#pragma once
#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <cstdio>
#include <map>
#include <math.h>
#include <WArray.h>
#include <Node.h>
#include <NWAligner.h>
#include <chrono>
#include <ctime>

class CoordMatcher
{
public:
	map<pair<int, int>, vector<pair<int, int>>> m;

	CoordMatcher(vector<WArray> warray)
	{
		for (int i = 0; i < warray.size(); i++)
		{
			AddTrack(warray[i]);
		}
	}

	void AddTrack(WArray w1)
	{
		pair<int, int> p1 = pair<int, int>(w1.value[0], w1.value[1]);
		pair<int, int> p2 = pair<int, int>(w1.value[2], w1.value[3]);
		m[p1].push_back(p2);
	}

	vector<pair<int, int>> Next(pair<int, int> c_begin)
	{
		if (m.find(c_begin) == m.end()) {
			return vector<pair<int, int>>();
		}
		else
		{ 
			// found
			return m[c_begin];
		}
		return vector<pair<int, int>>();
	}
};

class STrack
{

public:
	int loc_first, loc_last;
	pair<int, int> pos_first, pos_last;
	bool failed_elong;
	int length;
	int hash_length;

	STrack() {};
	STrack(int hash, int loc, pair<int, int> p_first, pair<int, int> p_last)
	{
		loc_first = loc;
		loc_last = loc + hash;
		pos_last = p_last;
		pos_first = p_first;
		failed_elong = false;
		hash_length = hash;
		length = hash;
	}

	void Increase(int loc, pair<int, int> p1, pair<int, int>p2)
	{
		loc_last = loc;
		pos_last = p2;
		length += hash_length-1;
		failed_elong = false;
	}
	void IncreaseTail(int loc, pair<int, int> p)
	{
		length += (loc - loc_last);
		loc_last = loc;
		pos_last = p;
		failed_elong = false;
	}

    void Expand(map<int, Node> &Body, string read)
	{
		int n_for = pos_last.first;
		int p_for = pos_last.second;
		chrono::time_point<chrono::system_clock> start1, end1;
		int elsec;
		start1 = chrono::system_clock::now();

		tuple<int, int, int> tmp_for = Body[n_for].checkThreadForward(p_for, read.substr(loc_last));
		
		end1 = chrono::system_clock::now();
		elsec = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
		cout << "Expand check Forward time" << elsec << '\n';


		int add_len_for = get<0>(tmp_for);
		n_for = get<1>(tmp_for);
		p_for = get<2>(tmp_for);

		int n_back = pos_first.first;
		int p_back = pos_first.second;
		tuple<int, int, int>tmp_back = Body[n_back].checkThreadBackward(p_back, read.substr(0, loc_first+1));

		int add_len_back = get<0>(tmp_back);
		n_back = get<1>(tmp_back);
		p_back = get<2>(tmp_back);

		
		length = length + add_len_back + add_len_for;
		loc_last = loc_last + add_len_for - 1;
		loc_first = loc_first - add_len_back + 1;

	}
};

class SBundle
{
public:
	vector<STrack> Tracks;
	bool failed_elong;
	bool cont;

	SBundle() {};
	SBundle(STrack strack)
	{
		Tracks.push_back(strack);
		failed_elong = false;
		cont = true;
	};
	
	vector<SBundle> GetIncreasedBundles(int loc, pair<int,int> p1, vector<pair<int, int>> conts)
	{
		vector<SBundle> res;
		for (int i = 0; i < conts.size(); i++)
		{
			SBundle tmp = *this;
			tmp.Increase(loc, p1, conts[i]);
			res.push_back(tmp);
		}
		return res;
	}

	void Increase(int loc, pair<int, int> p1, pair<int, int>p2)
	{
		for (int i = 0; i < Tracks.size(); i++)
		{
			Tracks[i].Increase(loc, p1, p2);
		}
		failed_elong = false;
	}

	void IncreaseTail(int loc, pair<int, int> p1)
	{
		for (int i = 0; i < Tracks.size(); i++)
		{
			Tracks[i].IncreaseTail(loc, p1);
		}
	}

	void AddTrack(STrack track)
	{
		Tracks.push_back(track);
	}

	void Merge(SBundle sb)
	{
		Tracks.insert(Tracks.end(), sb.Tracks.begin(), sb.Tracks.end());
		failed_elong = false;
	}

	vector<STrack> Clear(int hash_length)
	{
		vector<STrack> res;
		if (failed_elong == false) return res;
		for(int i =0; i<Tracks.size(); i++)
		{
				if (Tracks[i].length > hash_length * 2)
				{
					res.push_back(Tracks[i]);
				}
				Tracks.erase(Tracks.begin() + i);
				i--;
		}
		return res;
	}


	
};

class FinTrackPool
{
public:
	int seedkey;
	map<int, pair<bool, int>> values;
	FinTrackPool()
	{
	};

	void AddPoint(pair<int, int> p, bool seed = false) // ID/value
	{
		if (seed) seedkey = p.first;
		if (values[p.first].second < p.second)
		{
			values[p.first] = pair<bool, int>(true, p.second);
		}
	}

	void AddVectorPoint(vector<pair<int, int>> p_vec)
	{
		ReFlagVisited();
		for (int i = 0; i < p_vec.size(); i++)
		{
			AddPoint(p_vec[i]);
		}
	}

	void ReFlagVisited()
	{
		for (auto k : values)
		{
			values[k.first].first = false;
		}
	}

	bool NeedToUpdate()
	{
		
		for (auto k : values)
		{
			if (values[k.first].first)
				return true;
		}
		return false;
	}

	int GetMaxWeight()
	{
		int maxlen = values[seedkey].second;
		for (auto k : values)
		{
			if (k.second.second > maxlen) maxlen = k.second.second;
		}
		return maxlen;
	}

	pair<int, int>  GetLastPosition()
	{
		int maxlen = values[seedkey].second;
		pair<int, int> res = values[seedkey]; 
		for (auto k : values)
		{
			if (k.second.second > maxlen) maxlen = k.second.second;
		}
		
	}
};



class STrackList
{
public:
	int hash_length;
	vector<STrack> stored_path;

	vector<STrack> linked_blocks;


	//vector<STrack> current_path;
	map<pair<int, int>, SBundle> curr_path;
	map<pair<int, int>, SBundle> elonged_path;

	vector<vector<int>> sorted_by_end;
	vector<vector<int>> sorted_by_begin;

	vector<vector<int>> relations;

	FinTrackPool FinTracks;
	FinTrackPool FinBackTracks;

	STrackList() {};
	STrackList(int hash, int loc, CoordMatcher cm)
	{
		hash_length = hash;
		for (auto t : cm.m)
		{
			for (int i = 0; i < t.second.size(); i++)
			{
				STrack tmp(hash, loc, t.first, t.second[i]);
				curr_path[t.second[i]].AddTrack(tmp);
			}
		}
	}

	void AddPoint(int loc, pair<pair<int, int>, pair<int, int>> point)
	{
		STrack tmp(hash_length, loc, point.first, point.second);
		
		if (curr_path.find(point.second) == curr_path.end())
		{
			curr_path[point.second] = tmp;
		}
		else
		{
			curr_path[point.second].AddTrack(tmp);
		}
	}

	void AddPoint(int loc, pair<int, int> p1, pair<int, int> p2)
	{
		STrack tmp(hash_length, loc, p1, p2);
		SBundle sb(tmp);
		if (elonged_path.find(p2) == elonged_path.end())
		{
			elonged_path[p2] = sb;
		}
		else
		{
			elonged_path[p2].Tracks.push_back(tmp);
		}
	}

	void Reflag()
	{
		for (auto t : curr_path)
		{
			curr_path[t.first].failed_elong = true;
			curr_path[t.first].cont = false;
		}
	}

	void Clear()
	{
		vector<pair<int, int>> todelete;
		for (auto t : curr_path)
		{
			vector<STrack> res = curr_path[t.first].Clear(hash_length);
			stored_path.insert(stored_path.end(), res.begin(), res.end());
			if (curr_path[t.first].Tracks.size() == 0)
			{
				todelete.push_back(t.first);
			}
		}
		for (int i = 0; i < todelete.size(); i++)
		{
			curr_path.erase(todelete[i]);
		}
	}

	void Step(int loc, CoordMatcher cm)
	{
		for (auto p : cm.m)
		{
			if (curr_path.find(p.first) == curr_path.end())
			{
				for (int i = 0; i < p.second.size(); i++)
				{
					AddPoint(loc, p.first, p.second[i]);
				}
			}
			else
			{
				curr_path[p.first].cont = true;
				for (int i = 0; i < p.second.size(); i++)
				{
					SBundle tmp = curr_path[p.first];
					tmp.Increase(loc, p.first, p.second[i]);
					elonged_path[p.second[i]].Merge(tmp);
				}
				
			}
		}
		SaveBreakedBundles();
		curr_path = elonged_path;
		elonged_path.clear();
		Clear();
		Reflag();

	}

	void FinalStep(int loc, CoordMatcher cm, map<int, Node> &Body)
	{
		for (auto p : cm.m)
		{
			for (int i = 0; i < p.second.size(); i++)
			{
				for (auto sb : curr_path)
				{
					if (CanBeElonged(sb.second,
									p.first,
									p.second[i],
									Body))
					{
						SBundle tmp = sb.second;
						tmp.IncreaseTail(loc, p.second[i]);
						elonged_path[p.second[i]].Merge(tmp);
					}
				}
			}
		}
		SaveBreakedBundles();
		curr_path = elonged_path;
		elonged_path.clear();
		Clear();
		Reflag();
	}


	bool CanBeElonged(SBundle sb, pair<int, int> p1, pair<int, int> p2, map<int,Node> &Body)
	{
		pair<int, int> fpoint = sb.Tracks[0].pos_last;
		int rel_prev = GetNodeRelation(fpoint, p1, Body);
		int rel_next = GetNodeRelation(fpoint, p2, Body);

		if ((rel_prev == -1) && (rel_next == 1))
		{
			return true;
		}
		return false;
	}


	void SaveBreakedBundles()
	{
		for (auto t : curr_path)
		{
			if (curr_path[t.first].cont == false)
			{
				vector<STrack> res = curr_path[t.first].Clear(hash_length);
				stored_path.insert(stored_path.end(), res.begin(), res.end());
			}
		}
	}

	void InitRelations(map<int, Node> &Body)
	{
		int l = stored_path.size();
		relations = vector<vector<int>>(l);
		for (int i = 0; i < l; i++)
		{
			relations[i] = vector<int>(l);
		}


		for (int i = 0; i < l; i++)
		{
			for (int j = 0; j < l; j++)
			{
				relations[i][j] = GetSTracksRelation(i, j, Body);
			}
		}

	}

	int GetNodeRelation(pair<int, int> p1, pair<int, int> p2, map<int, Node> &Body)
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

	int GetSTracksRelation(int sT1ID, int sT2ID, map<int, Node> &Body)
	{
		
		int c1 = GetNodeRelation(stored_path[sT1ID].pos_last, stored_path[sT2ID].pos_first, Body);
		int c2 = GetNodeRelation(stored_path[sT1ID].pos_first, stored_path[sT2ID].pos_last, Body);
		if (c1 == 1) return 1;
		if (c2 == -1) return -1;
		return 0;
	}

	int GetGapEstimate(int sT1ID, int sT2ID)
	{
		return stored_path[sT1ID].loc_last - stored_path[sT2ID].loc_first;
	}

	int GetInitialSeed()
	{
		int SPl = stored_path.size();
		if (SPl == 0) return -1;
		int len = stored_path[0].length;
		int i_res = 0;
		for (int i = 0; i < SPl; i++)
		{
			if (stored_path[i].length > len)
			{
				i_res = i;
				len = stored_path[i].length;
			}
		}
		return i_res;
	}

	void SortByBeginEnds(int len)
	{
		sorted_by_begin = vector<vector<int>>(len);
		sorted_by_end = vector<vector<int>>(len);
		for (int i = 0; i < stored_path.size(); i++)
		{
			sorted_by_begin[stored_path[i].loc_first].push_back(i);
			sorted_by_end[stored_path[i].loc_last].push_back(i);
		}
	}

	int GetNWAligneWeight(map<int, Node> &Body, pair<int, int>n1, pair<int, int>n2, string read)
	{
		NWAligner NWAl;
		int err = NWAl.ExtractSubgraph(Body, n1, n2, read.size());
		if (err != 0) return 0;
		NWAl.MatrixInitNR(read);
		int res = NWAl.StartAlignmentNR(read);
		return res;
	}

	vector<pair<int,int>> GetListofPrecs(string read, map<int, Node> &Body, int seedId, int weight)
	{
		vector<pair<int, int>> res;
		for (int i = 0; i < stored_path.size(); i++)
		{
			int loc_dist = (stored_path[i].loc_first - stored_path[seedId].loc_last);
			if ((relations[seedId][i] == 1) &&
				(loc_dist > 0))
			{
				pair<int, int> n1 = stored_path[seedId].pos_last;
				pair<int, int> n2 = stored_path[i].pos_first;

				int frag_pos = stored_path[seedId].loc_last;
				int len = stored_path[seedId].loc_first - stored_path[i].loc_last;
				string read_fragment = read.substr(frag_pos, len);

				int NWweight = GetNWAligneWeight(Body, n1, n2, read_fragment);

				res.push_back(pair<int,int>(i, weight + stored_path[i].length - loc_dist +1+ NWweight));
			}
		}
		return res;
	}


	vector<pair<int, int>> GetListofAns(string read, map<int, Node> &Body, int seedId, int weight)
	{
		vector<pair<int, int>> res;
		for (int i = 0; i < stored_path.size(); i++)
		{
			int loc_dist = (stored_path[seedId].loc_first - stored_path[i].loc_last);
			if ((relations[seedId][i] == -1) &&
				(loc_dist > 0))
			{

				pair<int, int> n1 = stored_path[seedId].pos_first;
				pair<int, int> n2 = stored_path[i].pos_last;

				int frag_pos = stored_path[i].loc_last;
				int len = stored_path[i].loc_first - stored_path[seedId].loc_last;
				string read_fragment = read.substr(frag_pos, len);

				int NWweight = GetNWAligneWeight(Body, n2, n1, read_fragment);

				res.push_back(pair<int, int>(i, weight + stored_path[i].length -loc_dist + NWweight + 1));
			}
		}
		return res;
	}

	float LigationBurn(string read, map<int, Node> &Body)
	{
		int readlen = read.size();

		int seedID = GetInitialSeed();
		FinTracks.AddPoint(pair<int, int>(seedID, stored_path[seedID].length), true);
		bool needtoupdate = false;
		vector<pair<int, int>> nexts;
		do
		{
			for (auto k : FinTracks.values)
			{
				if (k.second.first)
				{
					vector<pair<int, int>> tmp = GetListofPrecs(read, Body, k.first, k.second.second);
					nexts.insert(nexts.end(), tmp.begin(), tmp.end());
				}
			}
			FinTracks.AddVectorPoint(nexts);
			nexts.clear();
			needtoupdate = FinTracks.NeedToUpdate();
		} 
		while (needtoupdate);

		
		FinBackTracks.AddPoint(pair<int, int>(seedID, stored_path[seedID].length), true);
		needtoupdate = false;
		do
		{
			for (auto k : FinBackTracks.values)
			{
				if (k.second.first)
				{
					vector<pair<int, int>> tmp = GetListofAns(read, Body, k.first, k.second.second);
					nexts.insert(nexts.end(), tmp.begin(), tmp.end());
				}
			}
			FinBackTracks.AddVectorPoint(nexts);
			nexts.clear();
			needtoupdate = FinBackTracks.NeedToUpdate();
		} while (needtoupdate);
		




		int TotalWeight = FinBackTracks.GetMaxWeight() + FinTracks.GetMaxWeight() - stored_path[seedID].length;

		float coverage = (float)TotalWeight / (float)readlen;

		return coverage;

	}


	void ExpandAll(map<int,Node> &Body, string read)
	{
		cout << "Expand size:" << stored_path.size() << endl;
		for (auto t = stored_path.begin(); t != stored_path.end(); t++)
		{
			t->Expand(Body, read);
		}
	}

};


class SAligner
{
public:
	bool failed = false;

	int hash_length;
	int start_loc;
	string read;
	STrackList STracks;
	vector<CoordMatcher> r_matches;
	float coverage;

	SAligner(int hash, string input_read, vector<vector<WArray>> warray, float cov)
	{
		coverage = cov;
		read = input_read;
		hash_length = hash;
		for (int i = 0; i < warray.size(); i++)
		{
			CoordMatcher tmp(warray[i]);
			r_matches.push_back(tmp);
		}

		int i = 0; 
		while (r_matches[i].m.size() == 0)
		{
			i++;
			if (i == r_matches.size())
			{
				failed = true;
				return;
			}
		}
		start_loc = i;
		STracks = STrackList(hash_length, start_loc, r_matches[i]);
	}

	int Run(map<int, Node> &Body)
	{

		int curr_loc= start_loc+hash_length - 1;
		int max_pos = r_matches.size()-1;
		if (curr_loc > max_pos)
		{
			STracks.Clear();
			return 0;
		}
		do
		{
			STracks.Step(curr_loc, r_matches[curr_loc]);
			curr_loc += hash_length-1;
		} while (curr_loc <= max_pos);

		if (curr_loc != max_pos)
		{
			STracks.FinalStep(max_pos, r_matches[max_pos], Body);
		}

		STracks.Clear();
		return 0;
	}

	int Ligation(int readlen, map<int, Node> &Body)
	{
		chrono::time_point<chrono::system_clock> start1, end1;
		int elsec;
		start1 = chrono::system_clock::now();
		//cout<<"SAligner: Ligation start"<<endl;
		if (STracks.stored_path.size() == 0) return -1;

		STracks.InitRelations(Body);

		end1 = chrono::system_clock::now();
		elsec = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
		cout << "Ligation: Initrelations " << elsec << '\n';

		start1 = chrono::system_clock::now();
		STracks.ExpandAll(Body, read);

		end1 = chrono::system_clock::now();
		elsec = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
		cout << "Ligation: Expand " << elsec << '\n';
		start1 = chrono::system_clock::now();

		float tmp_coverage = STracks.LigationBurn(read, Body);

		end1 = chrono::system_clock::now();
		elsec = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
		cout << "Ligation: Burn " << elsec << '\n';

		if (tmp_coverage >= coverage) return 0;

		return -1;
	}


};