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

using namespace std;
class GenerHash
{
public:
	vector<unsigned long long> p;
	int base;
	int len;

	unsigned long long currhash;
	int curr_len;

	string str;
	int id;
	int pos;
	int m_pos;

	bool formed;

	GenerHash(int base, int len)
	{
		this->base = base;
		this->len = len;
		for (int i = len - 1; i >= 0; i--)
		{
			p.push_back(pow(base, i));
		}
	}

	int ord(char c)
	{
		switch (c)
		{
		case 'a': return 1;
			break;
		case 'c': return 2;
			break;
		case 'g': return 3;
			break;
		case 't': return 4;
			break;
		case 'A': return 1;
			break;
		case 'C': return 2;
			break;
		case 'G': return 3;
			break;
		case 'T': return 4;
			break;
		default:
			cout << "cashe error. unknown simbol: " << c << endl;
			break;
		}

	}

	void Init(int base, int len)
	{
		this->str = "";
		this->curr_len = 0;
		this->currhash = 0;

		p.clear();

		this->base = base;
		this->len = len;
		for (int i = len - 1; i >= 0; i--)
		{
			p.push_back(pow(base, i));
		}
	}

	void ReInit(int Id, string str)
	{
		this->str = str;
		this->m_pos = str.length() - 1;

		this->curr_len = 0;
		this->currhash = 0;

		this->pos = 0;
		this->id = Id;
		this->formed = false;
	}

	int Inc(char c)
	{
		currhash += ord(c)*p[curr_len];
		curr_len++;
		if (curr_len == len)
		{
			formed = true;
			return 0;
		}
		return -1;
	}

	int Next(char c1)
	{
		int res;

		if (!formed)
		{
			res = Inc(c1);
		}
		else
		{
			if (pos<m_pos)
			{
				currhash = (currhash - ord(str[pos])*p[0])*base + ord(c1);
				pos++;
				res = 0;
			}
			else
			{
				res = -3;
			}
		}
		return res;
	}

	unsigned long long Hash(string str)
	{
		unsigned long long res = 0;
		for (int i = 0; i<len; i++)
		{
			res += ord(str[i])*p[i];
		}
		return res;
	}

	int left()
	{
		return len - curr_len;
	}
};
