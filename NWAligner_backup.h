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

class NWNode
{
public:
	string ref;

	int ID;    //id == -1 - final node of thread
	bool NaNbp;

	bool End;
	bool First;

	bool proseeded;

	map<int, NWNode*> Prev;
	map<int, NWNode*> Next;

	map<char, map<char, int>>* CharMatrix;

	NWNode() {};
	NWNode(Node N1, map<char, map<char, int>>* CharMatrixGlob)
	{
		ref = N1.str;
		ID = N1.ID;
		NaNbp = N1.NaNbp;
		First = false;
		End = false;
		proseeded = false;
		matrixloaded = false;
		CharMatrix = CharMatrixGlob;

	}

	NWNode(Node N1, bool first, int pos, map<char, map<char, int>>* CharMatrixGlob)
	{
		ID = N1.ID;
		NaNbp = N1.NaNbp;

		if (first)
		{
			ref.append(&N1.str[pos]);
			End = false;
			First = true;
		}
		else
		{
			ref = N1.str.substr(0, pos);
			First = false;
			End = true;
		}
		proseeded = false;
		matrixloaded = false;
		CharMatrix = CharMatrixGlob;
	}


	vector<vector<int>> NWmatrix;
	vector<int> backIDs;

	bool matrixloaded;

	void InitMatrix(string read)
	{
		if (First)
		{
			InitMatrix_first(read);
		}
		else
		{
			vector<vector<int>> lastcols;
			vector<int> lastIDs;
			for (auto t1 = Prev.begin(); t1 != Prev.end(); t1++)
			{
				if (!t1->second->matrixloaded)
					t1->second->InitMatrix(read);
				lastcols.push_back(t1->second->NWmatrix.back());
				lastIDs.push_back(t1->second->ID);
			}

			InitMatrix(read, lastcols, lastIDs);
		}

		matrixloaded = true;
	}

	void InitMatrix_first(string read)
	{
		//zeroRaw and zeroColomn - for Init
		NWmatrix = vector<vector<int>>(ref.size() + 1);
		for (int i = 0; i < NWmatrix.size(); i++)
		{
			NWmatrix[i] = vector<int>(read.size() + 1);
		}
		NWmatrix[0][0] = 0;
		for (int i = 1; i < NWmatrix.size(); i++)
		{
			NWmatrix[i][0] = -i;
		}
		for (int i = 1; i < NWmatrix[0].size(); i++)
		{
			NWmatrix[0][i] = -i;
		}

		for (int i = 1; i < NWmatrix.size(); i++)
		{
			for (int j = 1; j < NWmatrix[0].size(); j++)
			{
				NWmatrix[i][j] = GetScore(NWmatrix[i - 1][j], NWmatrix[i][j - 1], NWmatrix[i - 1][j - 1], ref[i - 1], read[j - 1]);
			}
		}


	}

	void InitMatrix(string read, vector<vector<int>> lastcolomns, vector<int> colIDs)
	{
		NWmatrix = vector<vector<int>>(ref.size() + 1);
		backIDs = vector<int>(read.size() + 1);

		for (int i = 0; i < NWmatrix.size(); i++)
		{
			NWmatrix[i] = vector<int>(read.size() + 1);
		}

		for (int i = 0; i < NWmatrix[0].size(); i++)
		{
			int tmp = lastcolomns[0][i];
			int tmp_id = colIDs[0];
			for (int j = 0; j < lastcolomns.size(); j++)
			{
				if (tmp < lastcolomns[j][i])
				{
					tmp = lastcolomns[j][i];
					tmp_id = colIDs[j];
				}
			}

			NWmatrix[0][i] = tmp;
			backIDs[i] = tmp_id;
		}

		for (int i = 1; i < NWmatrix.size(); i++)
		{
			NWmatrix[i][0] = NWmatrix[i - 1][0] - i;
		}

		for (int i = 1; i < NWmatrix.size(); i++)
		{
			for (int j = 1; j < NWmatrix[0].size(); j++)
			{
				NWmatrix[i][j] = GetScore(NWmatrix[i - 1][j], NWmatrix[i][j - 1], NWmatrix[i - 1][j - 1], ref[i - 1], read[j - 1]);
			}
		}
	}

	int GetScore(int Up, int Left, int Diag, char c1, char c2)
	{
		int gap = -1;
		vector<int> vars = { Up + gap,Left + gap, (*CharMatrix)[c1][c2] };
		int res = *max_element(vars.begin(), vars.end());
		return res;
	}





	int StartAlignment(string read)
	{
		int weight = 0;
		vector<char> totalres;
		GetTraceBackRecursive(read, read.size() + 1, totalres);
		reverse(totalres.begin(), totalres.end());
		int j = 0;
		for (int i = 0; i < totalres.size(); i++)
		{
			switch (totalres[i])
			{
			case 'u':
			{
				//cout << read[j];
				j++;
				break;
			}
			case 'd':
			{
				//cout << read[j];
				weight++;
				j++;
				break;
			}
			case 'l':
			{
				//cout << '-';
				break;
			}
			}
		}
		//cout << endl;
		//for (char c : totalres)
			//cout << c;
		return weight;

	}

	pair<int, int> GetTraceBack(string read, pair<int, int> start_pos, vector<char>& veclastNode)
	{
		int i = 0;
		int j = 0;

		if (!First)
		{
			i = NWmatrix.size() - 1;
			j = read.size();
		}
		else
		{
			i = NWmatrix.size() - 1;
			j = start_pos.second;
		}

		pair<int, int> res = pair<int, int>(0, 0);


		while ((i != 0) && (j != 0))
		{
			char dir = GetDirection(i, j);

			veclastNode.push_back(dir);
			switch (dir)
			{
			case 'd':
			{
				i--;
				j--;
				break;
			}
			case 'u':
			{
				j--;
				break;
			}
			case 'l':
			{
				i--;
				break;
			}
			}
		}

		if ((i != 0) && (j == 0))
		{
			for (int k = i; k <= 0; k--)
			{
				veclastNode.push_back('l');
			}
			i = 0;
		}
		else if ((i == 0) && (j != 0))
		{
			if (First)
			{
				for (int k = j; j <= 0; j--)
				{
					veclastNode.push_back('u');
				}
				res = pair<int, int>(0, 0);
			}
			else
			{
				res = pair<int, int>(backIDs[j], j);
			}
		}
		else if ((i == 0) && (j == 0))
		{
			if (First)
			{
				res = pair<int, int>(0, 0);
			}
			else
			{
				res = pair<int, int>(backIDs[j], j);
			}
		}
		return res;
	}

	void GetTraceBackRecursive(string read, int start_pos, vector<char>& veclastNode)
	{
		int i = NWmatrix.size() - 1;
		int j = start_pos - 1;

		while ((i != 0) && (j != 0))
		{
			char dir = GetDirection(i, j);

			veclastNode.push_back(dir);
			switch (dir)
			{
			case 'd':
			{
				i--;
				j--;
				break;
			}
			case 'u':
			{
				j--;
				break;
			}
			case 'l':
			{
				i--;
				break;
			}
			}
		}

		if ((i != 0) && (j == 0))
		{
			for (int k = i; k <= 0; k--)
			{
				veclastNode.push_back('l');
			}
			i = 0;
		}
		else if ((i == 0) && (j != 0))
		{
			if (First)
			{
				for (int k = j; j <= 0; j--)
				{
					veclastNode.push_back('u');
				}
				return;
			}
			else
			{
				Prev[backIDs[j]]->GetTraceBackRecursive(read, j, veclastNode);
			}
		}
		else if ((i == 0) && (j == 0))
		{
			if (First)
			{
				return;
			}
			else
			{
				Prev[backIDs[j]]->GetTraceBackRecursive(read, j, veclastNode);
			}
		}
	}

	char GetDirection(int i, int j)
	{
		vector<int> tiles = { NWmatrix[i - 1][j - 1], NWmatrix[i - 1][j], NWmatrix[i][j - 1] };
		auto it = max_element(tiles.begin(), tiles.end());
		int res = distance(tiles.begin(), it);
		switch (res)
		{
		case 0: return 'd';
			break;
		case 1: return 'l';
			break;
		case 2: return 'u';
			break;
		}
	}
};

class NWAligner
{
public:

	map<char, map<char, int>> CharMatrix;

	void InitCharMatrix()
	{
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
					CharMatrix[c1][c2] = 0;
				}
			}
		}
	}

	NWAligner()
	{
		InitCharMatrix();
	};

	float ApprxLen = 0.0;

	map<int, NWNode> NWBody;
	NWNode* endPtr;
	NWNode* beginPtr;

	int ExtractSubgraph(map<int, Node>& Body, pair<int, int> np1, pair<int, int> np2, int readlen)
	{

		if (np1.first == np2.first)
		{
			NWBody[np1.first] = NWNode(Body[np1.first], true, np1.second, &CharMatrix);
			beginPtr = &NWBody[np1.first];
			endPtr = &NWBody[np2.first];
		}
		else
		{
			NWBody[np1.first] = NWNode(Body[np1.first], true, np1.second, &CharMatrix);
			beginPtr = &NWBody[np1.first];
			NWBody[np2.first] = NWNode(Body[np2.first], false, np2.second, &CharMatrix);
			endPtr = &NWBody[np2.first];
		}

		ApprxLen = beginPtr->ref.size();
		vector<int> addedIds = { beginPtr->ID, endPtr->ID };

		int tid = 0;
		while (tid != addedIds.size())
		{
			NWNode* t1 = &(NWBody[addedIds[tid]]);
			if (!t1->proseeded)
			{

				Node tmp = Body[t1->ID];
				float tmp_len = 0.0f;
				int tmp_cnt = 0;
				for (auto n_tmp : Body[t1->ID].Next)
				{
					if (n_tmp.second->IsInside(Body[np1.first], Body[np2.first]))
					{

						if (NWBody.find(n_tmp.second->ID) == NWBody.end())
						{
							NWBody[n_tmp.first] = NWNode(Body[n_tmp.first], &CharMatrix);
						}
						t1->Next[n_tmp.second->ID] = &NWBody[n_tmp.first];
						NWBody[n_tmp.first].Prev[t1->ID] = t1;

						tmp_len += NWBody[n_tmp.first].ref.size();
						tmp_cnt++;

						addedIds.push_back(n_tmp.first);
					}
				}
				if (tmp_cnt != 0)
				{
					tmp_len = tmp_len / tmp_cnt;
					ApprxLen += tmp_len;
				}
				if (ApprxLen >= 4 * readlen)
				{
					return -1;
				}
				t1->proseeded = true;
			}
			tid++;;
		}
		return 0;//GraphExtractSuccesfull
	}

	void MatrixInit(string read)
	{
		endPtr->InitMatrix(read);
	}

	void MatrixInitNR(string read)
	{
		vector<int> matrixed;

		beginPtr->InitMatrix_first(read);
		matrixed.push_back(beginPtr->ID);
		int currID = beginPtr->ID;

		while (currID != endPtr->ID)
		{
			if (NWBody[currID].Next.size() != 0)
			{
				int currID_tmp = currID;
				for (auto t1 = NWBody[currID].Next.begin(); t1 != NWBody[currID].Next.end(); t1++)
				{
					t1->second->InitMatrix(read);
					matrixed.push_back(t1->second->ID);
					currID_tmp = t1->second->ID;
				}
				currID = currID_tmp;
			}
		}

		endPtr->InitMatrix(read);
	}

	int  StartAlignment(string read)
	{
		if (read.size() * 4 < ApprxLen) return 0;
		int res = endPtr->StartAlignment(read);

		return res;
	}

	int StartAlignmentNR(string read)
	{
		if (read.size() * 4 < ApprxLen) return 0;
		vector<char> totalres;
		pair<int, int> next_pos;

		int currID = endPtr->ID;
		while (currID != beginPtr->ID)
		{
			next_pos = NWBody[currID].GetTraceBack(read, next_pos, totalres);
			currID = next_pos.first;
		}
		next_pos = NWBody[currID].GetTraceBack(read, next_pos, totalres);

		int weight = 0;
		reverse(totalres.begin(), totalres.end());
		int j = 0;
		for (int i = 0; i < totalres.size(); i++)
		{
			switch (totalres[i])
			{
			case 'u':
			{
				//cout << read[j];
				j++;
				break;
			}
			case 'd':
			{
				//cout << read[j];
				weight++;
				j++;
				break;
			}
			case 'l':
			{
				//cout << '-';
				break;
			}
			}
		}
		//	cout << endl;
		//	for (char c : totalres)
		//		cout << c;
		return weight;
	}

	int AlignRead(string read, int ID1, int pos1, int ID2, int pos2, map<int, Node>& Body)
	{

	}
}; 
