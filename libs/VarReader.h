#pragma once
#include <string>
#include <fstream>
#include <sstream>



struct Block
{
	int coord1, coord2;
	char dir;
	void reset()
	{
		coord1 = -1;
		coord2 = -1;
		dir = ' ';
	}
};
class SnpBlock
{
public:
	int lastpos;
	vector<tuple<int, int, string>> data;
	SnpBlock()
	{
		lastpos = 0;
	};
	void Add(tuple<int, int, string> t1)
	{
		data.push_back(t1);
		lastpos = get<1>(t1);
	}

	tuple<int, int, string> getSnpBlock()
	{
		string str = "";
		for (int i = 0; i<data.size(); i++)
		{
			str = str + get<2>(data[i]);
		}
		int pos1 = get<0>(data[0]);
		int pos2 = lastpos;

		return tuple<int, int, string>(pos1, pos2, str);
	}

	void clear()
	{
		data.clear();
	}
};

class MaffReader
{
public:
	int state;
	string var;
	char last_c;
	int pos1, pos2, c_pos;
	vector<tuple<int, int, string>> *res;

	MaffReader(char c, vector<tuple<int, int, string>> *res_init, int offset)
	{
		res = res_init;
		last_c = c;
		pos1 = 0;
		pos2 = 0;
		c_pos = offset - 1;
		var = "";
		state = 0;
	};

	void ReInit(int offset)
	{
		pos1 = 0;
		pos2 = 0;
		c_pos = offset - 1;
		var = "";
		state = 0;
	}

	void Next(char r1, char r2, char c)
	{
		if (r1 != '-')
		{
			c_pos++;
			if (r2 != '-')
			{
				last_c = r2;
			}
		}

		switch (c)
		{
		case '*':
		{
			if (var.length() != 0)
			{
				pos2 = c_pos;
				res->push_back(tuple<int, int, string>(pos1, pos2, var));
				var = "";
			}
			break;
		}
		case '.':
		case ' ':
		{
			if (r1 == '-')
			{
				if (var.length() == 0)
				{
					pos1 = c_pos;
					var += last_c;
					var += r2;
				}
				else
				{
					var += r2;
				}
			}
			else if (r2 == '-')
			{
				if (var.length() == 0)
				{
					pos1 = c_pos - 1;
					var += last_c;
				}
				else {};
			}
			else if (r1 != r2)
			{
				if (var.length() == 0)
				{
					pos1 = c_pos;
				}
				var += r2;
			}
			break;
		}
		}
	}


};
class VariationData
{
public:


	vector<pair<Block, Block>> data;
	vector<tuple<int, int, string>> res;
	string ref1, ref2;
	VariationData()
	{
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

	bool my_sort_function(tuple<int, int, string> &t1, tuple<int, int, string> &t2)
	{
		return get<0>(t1)<get<0>(t2);
	}

	void load(string path)
	{

		Block bl1;
		Block bl2;
		bl1.reset();
		bl2.reset();


		ifstream fin(path);
		string buf;

		int state = 0;

		do
		{
			getline(fin, buf);
			switch (buf[0])
			{
			case '>':
			{
				vector<string> tmp1 = split(buf, ' ');
				vector<string> tmp2 = split(tmp1[1], ':');
				vector<string> tmp3 = split(tmp2[1], '-');

				if (state == 0)
				{
					bl1.dir = tmp1[2][0];
					bl1.coord1 = stoi(tmp3[0]) - 1;
					bl1.coord2 = stoi(tmp3[1]) - 1;
					state = 1;
				}
				else
				{
					bl2.dir = tmp1[2][0];
					bl2.coord1 = stoi(tmp3[0]) - 1;
					bl2.coord2 = stoi(tmp3[1]) - 1;
				}
				break;
			}
			case '=':
			{
				if (bl1.dir == bl2.dir)
				{
					data.push_back(pair<Block, Block>(bl1, bl2));
				}
				bl1.reset();
				bl2.reset();
				state = 0;
				break;
			}
			case '#':
				break;
			default:
				break;
			}
		} while (!fin.eof());


	}

	int len(string str)
	{
		int res = str.length();
		for (int i = 0; i<str.length(); i++)
		{
			if (str[i] == '.') res--;
		}
		return res;
	}


	void load_block_snp(string path, int offset1, bool indel)
	{


		string buf = "_";
		string buf2 = "_";
		ifstream fin(path);

		int off = 1;
		if (indel)
		{
			off = 1;
		}

		do
		{
			getline(fin, buf, '\t');

			if (buf.length() == 0)
			{
				continue;
			}

			getline(fin, buf, '\t');
			int pos1 = stoi(buf) + offset1 - off;
			getline(fin, buf, '\t');
			int pos2 = pos1 + buf.length();
			getline(fin, buf, '\t');
			transform(buf.begin(), buf.end(), buf.begin(), (int(*)(int))tolower);
			getline(fin, buf2);

			if (buf2 == "DELETION")
			{
				pos1 = pos1 - 1;
				pos2 = pos2 - 1;
			}
			res.push_back(tuple<int, int, string>(pos1, pos2, buf));




		} while (!fin.eof());

	}

	void load_block_snp(MaffReader &mr1, string path, int offset)
	{
		mr1.ReInit(offset);
		string r1 = "_";
		string r2 = "_";
		string s = "";
		ifstream fin(path);

		getline(fin, r1);
		getline(fin, r1);
		getline(fin, r1);
		while (!fin.eof())
		{
			getline(fin, r1);
			getline(fin, r2);
			getline(fin, s);
			r1 = split(r1, ' ').back();
			r2 = r2.substr(r2.length() - r1.length(), r1.length());
			s = s.substr(s.length() - r1.length(), r1.length());
			for (int i = 0; i<s.length(); i++)
			{
				mr1.Next(r1[i], r2[i], s[i]);
			}
			getline(fin, s);
		}
	}

	vector<tuple<int, int, string>> getVCF(string ref1, string ref2, string path)
	{

		MaffReader mr1('_', &res, -1);

		if (data[0].first.coord1 = 0)
		{
			//первый блок одновременно с началом первого
		}

		for (int i = 0; i<data.size() - 1; i++)//переход по блокам от первого до предпоследнего
		{
			vector<tuple<int, int, string>> tmp;
			//	load_block_snp(path+"block"+to_string(i)+"_snp.vcf", data[i].first.coord1,false);
			//	load_block_snp(path+"block"+to_string(i)+"_indel.vcf", data[i].first.coord1,true);

			load_block_snp(mr1, path + "block" + to_string(i) + ".mft", data[i].first.coord1);
			int t1 = data[i].first.coord2 + 1;
			int t2 = data[i + 1].first.coord1;
			string t3 = ref2.substr(data[i].second.coord2 + 1, data[i + 1].second.coord1 - data[i].second.coord2 - 1);
			res.push_back(tuple<int, int, string>(t1, t2, t3));
		}

		/*	load_block_snp(path+"block"+to_string(data.size()-1)+"_snp.vcf",data[data.size()-1].first.coord1,false);
		load_block_snp(path+"block"+to_string(data.size()-1)+"_indel.vcf",data[data.size()-1].first.coord1,true);*/
		load_block_snp(mr1, path + "block" + to_string(data.size() - 1) + ".mft", data[data.size() - 1].first.coord1);


		if (data[data.size() - 1].first.coord2 != ref1.size()) //конец второго блока
		{

		}



		sort(res.begin(), res.end());

		res = prepost_vcf(res);

		for (int i = res.size() - 1; i>res.size() - 30; i--)
		{
			//cout << get<0>(res[i]) << '\t' << get<1>(res[i]) << '\t' << get<2>(res[i]) << endl;
		}

		return res;
	}


	vector<tuple<int, int, string>> prepost_vcf(vector<tuple<int, int, string>> &vcf)
	{
		vector<tuple<int, int, string>> tmp_res;

		int pos1 = get<0>(vcf[0]);
		int pos2 = get<1>(vcf[0]);
		string str = get<2>(vcf[0]);

		for (int i = 1; i<vcf.size(); i++)
		{
			if (get<0>(vcf[i]) == pos2)
			{
				pos2 = get<1>(vcf[i]);
				str = str + get<2>(vcf[i]);
			}
			else
			{
				tmp_res.push_back(tuple<int, int, string>(pos1, pos2, str));
				pos1 = get<0>(vcf[i]);
				pos2 = get<1>(vcf[i]);
				str = get<2>(vcf[i]);
			}
		}



		return tmp_res;
	}


}; 
