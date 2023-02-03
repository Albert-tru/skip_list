#pragma once

#ifndef EDGE_H
#define EDGE_H

#include"Point.h"

class Edge
{
public:
	int p1;		//�ߵĶ˵������
	int p2;
	bool isBad = false;
	bool isBound = true;	//�������Ƿ�Ϊ�߽�
	Edge() {

	}
	Edge(const int p1, const int p2)
		:p1(p1), p2(p2)
	{
	}

	bool operator<(const Edge& other)const   //����<�����
	{
		if (this->p1 < other.p1)
			return true;
		else if (this->p1 == other.p1)
		{
			if (this->p2 < other.p2)
				return true;
		}
		return false;
	}

	~Edge(){
		//cout << this->p1 << " " << this->p2 << "destroyed ! " << endl;
	}
};


#endif // !EDGE_H

