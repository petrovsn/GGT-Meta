#pragma once
#include <tuple>;
#include <Node.h>
using namespace std;

class WArray
{
public:
	int startPos;
	int endPos;
	vector<int> NodeIDs;

	bool assembled;

	WArray(int ID, int pos)
	{
		startPos = pos;
		NodeIDs.push_back(ID);
	}

	void extend(int ID)
	{
		NodeIDs.push_back(ID);
	}

	void move_right()
	{
		startPos++;
	}

	void finish(int pos)
	{
		endPos = pos;
	}
}; 

class TArray
{
public:
	int startPos;
	int endPos;
	vector<int> NodeIDs;

	bool assembled;
	int read_startPos;
	int read_endPos;

	TArray(WArray add_wray, int read_sPos, int read_ePos)
	{
		assembled = false;
		startPos = add_wray.startPos;
		endPos = add_wray.endPos;
		NodeIDs = add_wray.NodeIDs;
		read_startPos = read_sPos;
		read_endPos = read_ePos;
	}

	void extend(int ID)
	{
		NodeIDs.push_back(ID);
	}

	void move_right()
	{
		startPos++;
	}

	void finish(int pos)
	{
		endPos = pos;
	}

	int TryAdd(WArray& add_wray)
	{
		if ((endPos == add_wray.startPos) and (NodeIDs.back() == add_wray.NodeIDs[0]))
		{
			if (add_wray.NodeIDs.size() > 1)
			{
				for (int i = 1; i < add_wray.NodeIDs.size(); i++)
				{
					NodeIDs.push_back(add_wray.NodeIDs[i]);
				}
			}
			add_wray.assembled = true;
			endPos = add_wray.endPos;
			return 0;
		}
		return -1;

	}

	int TryExtRight(WArray& add_wray, map<int, Node>& Body, int read_pos, int read_ePos)
	{
		//case one, same node, right position
		int lastID = NodeIDs.back();
		int lastID_new = add_wray.NodeIDs.back();
		if ((lastID == lastID_new) and (endPos == add_wray.endPos - 1))
		{
			endPos = endPos + 1;
			add_wray.assembled = true;
			read_endPos = read_ePos;
			return 0;
		}


		//case two, difference node
		int curnodelen = Body[lastID].str.size();
		bool checkside = false;

		if (curnodelen - 1 == endPos)
		{
			for (auto node : Body[lastID].Next)
			{
				if (lastID_new == node.first)
				{
					checkside = true;
				}
			}
		}

		if (checkside)
		{
			endPos = 0;
			NodeIDs.push_back(lastID_new);
			add_wray.assembled = true;
			read_endPos = read_ePos;
			return 0;
		}


		return -1;

	}

	int TryExtLeft(WArray& add_wray, map<int, Node>& Body)
	{
		//case one, same node, left position
		int lastID = NodeIDs[0];
		int lastID_new = add_wray.NodeIDs[0];
		if ((lastID == lastID_new) and (startPos == add_wray.startPos + 1))
		{
			startPos = startPos - 1;
			add_wray.assembled = true;
			return 0;
		}

		//case two, difference node
		bool checkside = false;

		if (0 == startPos)
		{
			for (auto node : Body[lastID].Prev)
			{
				if ((lastID_new == node.first) and (add_wray.startPos == node.second->str.size() - 1))
				{
					checkside = true;
				}
			}
		}

		if (checkside)
		{
			startPos = add_wray.startPos;
			vector<int> IDs_tmp;

			IDs_tmp.push_back(lastID_new);
			for (int i = 0; i < NodeIDs.size(); i++)
			{
				IDs_tmp.push_back(NodeIDs[i]);
			}
			NodeIDs = IDs_tmp;
			add_wray.assembled = true;
			return 0;
		}


		return -1;
	}

	bool IntersecSeed(TArray &tr, map<int, Node>& Body)
	{
		if (read_startPos < tr.read_startPos)
		{
			
			bool comp1 = Node1LessNode2(NodeIDs[0], startPos, tr.NodeIDs[0], tr.startPos, Body);
			bool comp2 = Node1LessorEqualNode2(tr.NodeIDs[0], tr.startPos, NodeIDs.back(), endPos, Body);
			if (comp1 and comp2)
			{
				return true;
			}
		}
		return false;
	}

	void CutFromLeft(TArray& tr_left, map<int, Node>& Body)
	{
		int F_ID = tr_left.NodeIDs.back();
		int F_pos = tr_left.endPos;

		if (F_ID == NodeIDs[0])
		{
			int dif = F_pos - startPos+1;
			read_startPos = read_startPos + dif;
			startPos = startPos + dif;
		}
		else
		{
			int pnt = 0;
			int currID = NodeIDs[pnt];
			int dist = Body[currID].str.length() - startPos;
			pnt = 1;
			currID = NodeIDs[pnt];

			while (currID != F_ID)
			{
				dist = dist + Body[currID].str.length();
				pnt = pnt + 1;
				currID = NodeIDs[pnt];
			}

			dist = dist + F_pos + 1;
			read_startPos = read_startPos + dist;
			startPos = F_pos + 1;
			vector<int> NewIds;
			for (int i = pnt; i < NodeIDs.size(); i++)
			{
				NewIds.push_back(NodeIDs[i]);
			}
			NodeIDs = NewIds;
		}
	}
};

