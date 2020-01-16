#pragma once

#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <cstdio>
#include <map>
#include <math.h>
#include <Node.h>


class TBlock
{
public:
	int first;
	int last;
	int loc1, loc2; //in node
	bool Inced;
	TBlock(int p, int l)
	{
		loc1 = l;
		loc2 = l;
		first = p;
		last = p;
		Inced = true;
	};
	void Inc()
	{
		last++;
		loc2++;
		Inced = true;
	};

	
};


class Tracker
{
public:
	vector<int> ch_nodes;
	int first;
	int last;

	Tracker(int n, int p)
	{
		ch_nodes.push_back(n);
		first = p;
		last = p;
	}



};


class NodeRaw
{
public:
	int max_size;
	vector<vector<int>> dots;
	vector<TBlock> tblocks;
	vector<int> nexts;
	
	NodeRaw()
	{

	}

	NodeRaw(Node &node)
	{
		max_size = node.str.length();
		dots = vector<vector<int>>(max_size, vector<int>(0));
		
		for (auto b : node.Next)
		{
			nexts.push_back(b.first);
		}

	}
	
	void AddDot(int pos, int loc)
	{
		dots[pos].push_back(loc);
	}

	void Run()
	{
		for (int i = 0; i < dots.size(); i++)
		{
			ExpTblockVector(dots[i],i);
			ClearTBlocks();
		}
	}

	void ExpTblockVector(vector<int> &d, int loc)
	{
		for (int i = 0; i < d.size(); i++)
		{
			bool added = false;
			for (int j = 0; j < tblocks.size(); j++)
			{
				if ((d[i] == (tblocks[j].last + 1))&&(!tblocks[j].Inced))
				{
					tblocks[j].Inc();
					added = true;
					continue;
				}
			}
			if (!added)
			{
				tblocks.push_back(TBlock(d[i],loc));
			}
			
		}
	}

	void ClearTBlocks()
	{
		for (int i = 0; i < tblocks.size(); i++)
		{
			if ((!tblocks[i].Inced)&&((tblocks[i].last - tblocks[i].first)<2))
			{
				tblocks.erase(tblocks.begin()+i);
				i--;
			}
		}
		for (int i = 0; i < tblocks.size(); i++)
		{
			tblocks[i].Inced = false;
		}
	}

	int CheckTblocks(int pos, bool ended)
	{
		for (int i = 0; i < tblocks.size(); i++)
		{
			
			if (tblocks[i].first == pos)
			{
				if (ended)
				{
					if (tblocks[i].loc2 == max_size - 1)
					{
						return tblocks[i].first;
					}
				}
				else
				{
					return tblocks[i].first;
				}
			}
			
		}

		return -1;
	}
};



class NodeList
{
public:
	map<int, NodeRaw> plot;

	vector<Tracker> tracks;


	void AddPoint(Node &node, int pos, int loc)//pos - position in NodeRaw, loc - location of substring in read;
	{
		if (plot.find(node.ID) == plot.end())//not exists;
		{
			NodeRaw tmpNodeRaw(node);
			plot[node.ID] = tmpNodeRaw;
		}
		plot[node.ID].AddDot(pos, loc);
	};

	void TRunner()
	{
		for (auto n : plot)
		{
			plot[n.first].Run();
		}
	}

	void TAligner()
	{
		for (auto n : plot)
		{
			int q = plot[n.first].CheckTblocks(0, true);
			if (q != -1)
			{
				tracks.push_back(Tracker(n.first, q));
			}
		}
	}
};