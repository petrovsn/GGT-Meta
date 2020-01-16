#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#include <sstream>

using namespace std;

string loadfasta(string path)
{
	string buf;
	ifstream fin(path);

	string head;
	char c;
	getline(fin, head);
	while (!fin.eof())
	{
		c = tolower(fin.get());
		if ((c == 'a') || (c == 'c') || (c == 'g') || (c == 't')) buf += c;
	}
	//cout << head << '\t' << buf.size() << '\t';
	return buf;
}


string reverse(string read)
{
	string rev_read = "";
	for (int i = read.length() - 1; i >= 0; i--)
	{
		switch (read[i])
		{
		case 'a':
		case 'A':
		{
			rev_read += 't';
			break;
		}
		case 't':
		case 'T':
		{
			rev_read += 'a';
			break;
		}
		case 'c':
		case 'C':
		{
			rev_read += 'g';
			break;
		}
		case 'g':
		case 'G':
		{
			rev_read += 'c';
			break;
		}
		}
	}
	return rev_read;
}

vector<string> split(string str, char delimiter)
{
	vector<string> internal;
	stringstream ss(str); // Turn the string into a stream.
	string tok;
	while (getline(ss, tok, delimiter))
	{
		internal.push_back(tok);
	}

	return internal;
}

vector<pair<string, string>> loadmultifasta(string path)
{
	vector<pair<string, string>> res;

	ifstream fin(path);
	string head;
	string read;

	getline(fin, head);
	getline(fin, read);
	while (!fin.eof())
	{
		
		//head.erase(0);
		auto tmp = pair<string, string>(head, read);
		res.push_back(tmp);
		getline(fin, head);
		getline(fin, read);
	}

	return res;
}

vector<pair<string, string>> loadmultifastQ(string path)
{
	vector<pair<string, string>> res;


	return res;
}