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
#include <algorithm>
#include <GenerHash.h>
#include <BubbleIndex.h>

class Node
{
public:
	int ID;    //id == -1 - final node of thread
	string  str;
	bool NaNbp;
	bool End;
	map<int, Node*> Prev;
	map<int, Node*> Next;
	unordered_map<int, int> B; //threads 

	BubbleIndex bidx;

	Node();

	Node(int id);

	Node(string str, int id);

	const Node operator = (Node& n);

	Node Split(int pos, int id);

	void Split(int pos, Node* n_pr);

	void Link(Node *node);

	void GenBIdx(BubbleIndex bidx_last, int c, bool p);

	int GetRelation(Node &n);

	bool IsInside(Node &n1, Node &n2);

	bool IsInside(int N1, int N2, map<int, Node> &Body);

	void GetTrackR(int i, int len, vector<tuple<string, int, int>> &res, string last_track = "")
	{
		if (End)
		{
			return;
		}
		if (len - last_track.size() < str.size()-i)
		{
			res.push_back(tuple<string, int, int>(last_track, ID, i + (len - last_track.size())));
			last_track += str.substr(i, len - last_track.size());
		}
		else
		{

			last_track += str.substr(i);
			for (auto t = Next.begin(); t != Next.end(); t++)
			{
				t->second->GetTrackR(0, len, res, last_track);
			}
		}
		
	}

/*	void GetBackTrackR(int i, int len, vector<tuple<string,int,int>> &res, string last_track = "")
	{
		int offset = len - last_track.size();
		if (i - len >= 0)
		{
			int p1 = (i - offset);
			last_track = str.substr(p1,offset) + last_track;
			res.push_back(tuple<string,int,int>(last_track, ID, p1));
		}
		else
		{
			last_track = str.substr(0, i+1) + last_track;
			for (auto t = Prev.begin(); t != Prev.end(); t++)
			{
				t->second->GetBackTrackR(t->second->str.size()-1, len, res, last_track);
			}
		}
	}*/

	tuple<int,int,int> checkThreadForward(int p, string read)
	{
		tuple<int, int, int> res;
		int final_count = 0;
		int pos = 0;
		if (read.size() <= str.substr(p).size())
		{
			while (pos < read.size())
			{
				if (read[pos] == str[p + pos])
				{
					final_count++;
					pos++;
				}
				else
				{
					res =  tuple<int,int,int>(final_count-1,ID,p+pos-1);
					return res;
				}
			}
			res = tuple<int, int, int>(final_count - 1, ID, p + pos - 1);
			return res;
			
		}
		else if (read.size() > str.substr(p).size())
		{
			while (pos < str.size()-p)
			{
				if (read[pos] == str[p + pos])
				{
					final_count++;
					pos++;
				}
				else
				{
					res =  tuple<int, int, int>(final_count, ID, p + pos-1);
					return res;
				}
			}

			int count = 0;
			for (auto t : Next)
			{
				tuple<int,int,int> tmp_res = t.second->checkThreadForward(0, read.substr(pos));
				int tmp = get<0>(tmp_res);
				if (tmp > count)
				{
					count = tmp;
					res = tmp_res;
				}
			}
			
			final_count += count;
			int fin_ID = get<1>(res);
			int fin_pos = get<2>(res);
			res = tuple<int, int, int>(final_count, fin_ID, fin_pos);
		}
		return res;
	}

	tuple<int, int, int> checkThreadBackward(int p, string read)
	{
		tuple<int, int, int> res;
		int final_count = 0;
		int rlen = read.size();
		int pos = read.size()-1;
		if (read.size() <= p+1)
		{
			while (pos >= 0)
			{
				if (read[pos] == str[(p + 1) - rlen+pos])
				{
					final_count++;
					pos--;
				}
				else
				{
					res = tuple<int, int, int>(final_count, ID, (p + 1) - rlen + pos + 1);
					return res;
				}
			}
			res = tuple<int, int, int>(final_count + 1, ID, (p + 1) - rlen + pos + 1);
			return res;
		}
		else 
		{
			while (pos>rlen-p-2)
			{
				if (read[pos] == str[p -rlen+pos+1])
				{
					final_count++;
					pos--;
				}
				else
				{
					res = tuple<int, int, int>(final_count, ID, (p + 1) - rlen + pos + 1);
					return res;
				}
			}

			if (pos < 0)
			{
				res = tuple<int, int, int>(final_count, ID, (p + 1) - rlen + pos + 1);
				return res;
			}


			int count = 0;
			for (auto t : Prev)
			{
				tuple<int, int, int> tmp_res = t.second->checkThreadBackward(t.second->str.size()-1, read.substr(0, pos+1));
				int tmp = get<0>(tmp_res);
				if (tmp > count)
				{
					count = tmp;
					res = tmp_res;
				}
			}

			final_count += count;
			int fin_ID = get<1>(res);
			int fin_pos = get<2>(res);
			res = tuple<int, int, int>(final_count, fin_ID, fin_pos);
		}
		return res;
	}
}; 

bool Node1LessNode2(int ID1, int pos1, int ID2, int pos2, map<int, Node>& Body);

bool Node1LessorEqualNode2(int ID1, int pos1, int ID2, int pos2, map<int, Node>& Body);