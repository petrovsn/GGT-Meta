#pragma once
#include <vector>

using namespace std;
class BubbleIndex
{
public:

	vector<pair<int,int>> idx;
	bool loaded;

	BubbleIndex(bool first = false)
	{
		loaded = false;
		if (first)
		{
			idx.push_back(pair<int, int>(0, 0));
			loaded = true;
		}
	};

	BubbleIndex(BubbleIndex bidx_last, int c, bool p)
	{
		idx = bidx_last.idx;
		if (p)//Bubble branch
		{
			idx.back().first++;
			idx.back().second = c;
			idx.push_back(pair<int, int>(0, 0));

		}
		else //collapse point
		{
			idx.pop_back();
			idx.back().first++;
			idx.back().second = 0;
		}
		loaded = true;
	}
	
	int Compare(const BubbleIndex &b) //  -1 -prev; 0 - parralel, 1 - is the next, 3 - is the same, -2 - error;
	{
		if ((loaded == false) or (b.loaded == false))
		{
			return -2;
		}
		for (int i = 0; i < min(idx.size(), b.idx.size()); i++)
		{
			if (idx[i].first < b.idx[i].first) return 1;
			else if (idx[i].first > b.idx[i].first) return -1;
			else if (idx[i].second != b.idx[i].second) return 0;
		}

		return 3;
	}

};