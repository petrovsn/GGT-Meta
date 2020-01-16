#include <Node.h>


Node::Node()
{
	this->ID = -1;
	this->str = "";
	this->End = false;
	this->NaNbp = false;
}

Node::Node(int id)
{
	this->ID = id;
	this->str = "";
	this->End = false;
	this->NaNbp = false;
}

Node::Node(string str, int id)
{
	this->str = str;
	this->ID = id;
	this->End = false;
	if (str[0] == 'n') this->NaNbp = true;
	else this->NaNbp = false;
}

const Node Node::operator = (Node& n)
{
	this->ID = n.ID;
	this->str = n.str;

	this->B = n.B;
	this->Next = n.Next;
	this->Prev = n.Prev;

	return *this;
}

Node Node::Split(int pos, int id)
{
	string tmp = this->str.substr(pos, str.length() - pos);
	this->str = this->str.substr(0, pos);

	Node n_pr(tmp, id);

	(n_pr.Prev).insert(pair<int, Node*>(this->ID, this));



	for (auto t1 : this->Next)
	{
		t1.second->Prev.erase(this->ID);
		t1.second->Prev.insert(pair<int, Node*>(id, &n_pr));
		n_pr.Next.insert(t1);
	}



	Next.clear();
	Next.insert(pair<int, Node*>(id, &n_pr));

	return n_pr;
}

void Node::Split(int pos, Node* n_pr)
{
	string tmp = this->str.substr(pos, str.length() - pos);
	this->str = this->str.substr(0, pos);

	n_pr->str = tmp;
	n_pr->NaNbp = this->NaNbp;

	(n_pr->Prev).insert(pair<int, Node*>(this->ID, this));

	for (auto t1 : this->Next)
	{
		t1.second->Prev.erase(this->ID);
		t1.second->Prev.insert(pair<int, Node*>(n_pr->ID, n_pr));
		n_pr->Next.insert(t1);
	}

	Next.clear();
	Next.insert(pair<int, Node*>(n_pr->ID, n_pr));



	for (auto t2 : this->B)
	{
		n_pr->B.insert(t2);
		B[t2.first] = n_pr->ID;
	}

}

void Node::Link(Node *N2)
{
	(this->Next).insert(pair<int, Node*>(N2->ID, N2));
	(N2->Prev).insert(pair<int, Node*>(this->ID, this));
}


void Node::GenBIdx(BubbleIndex bidx_last, int c, bool p)
{
	if ((bidx.loaded) || (ID == -1)) return;
	bidx = BubbleIndex(bidx_last, c, p);
}

int Node::GetRelation(Node &n)
{
	return bidx.Compare(n.bidx); //  -1 -prev; 0 - parralel, 1 - is the next, 3 - is the same
}

bool  Node::IsInside(Node &n1, Node &n2)
{
	int n1rel = this->GetRelation(n1);
	int n2rel = this->GetRelation(n2);

	if (((n1rel == -1)||(n1rel == 3))
		&& ((n2rel == 1) || (n2rel == 3)))
		return true;
	return false;
}

bool Node::IsInside(int N1, int N2, map<int, Node> &Body)
{
	int n1rel = this->GetRelation(Body[N1]);
	int n2rel = this->GetRelation(Body[N2]);

	if (((n1rel == -1) || (n1rel == 3))
		&& ((n2rel == 1) || (n2rel == 3)))
		return true;
	return false;
}


bool Node1LessNode2(int ID1, int pos1, int ID2, int pos2, map<int, Node>& Body)
{
	if (ID1 == ID2)
	{
		return pos1 < pos2;
	}
	else
	{
		int rel = Body[ID1].GetRelation(Body[ID2]);
		if (rel == 1)
		{
			return true;
		}
	}
	return false;
}

bool Node1LessorEqualNode2(int ID1, int pos1, int ID2, int pos2, map<int, Node>& Body)
{
	if (ID1 == ID2)
	{
		return pos1 <= pos2;
	}
	else
	{
		int rel = Body[ID1].GetRelation(Body[ID2]);
		if (rel == 1)
		{
			return true;
		}
	}
	return false;
}
