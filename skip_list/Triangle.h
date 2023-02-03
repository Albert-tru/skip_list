#pragma once
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include"Edge.h"
#include<vector>
#include<iostream>

using namespace std;

struct Triangle {

	int a;		//三角形端点的索引,顺时针方向
	int b;
	int c;
	bool isBad = false;
	//Edge edges[3];

	Triangle() = default;
	Triangle(const Triangle&) = default;
	Triangle(Triangle&&) = default;
	Triangle(const int a, const int b, const int c)
		:a(a), b(b), c(c)
	{
		this->a = a;
		this->b = b;
		this->c = c;

		//Edge e1(a, b);
		//Edge e2(b, c);
		//Edge e3(c, a);
		//cout << e1.p1 << endl;
		//edges[0].p1 = a, edges[0].p2 = b;
		//edges[1].p1 = b, edges[1].p2 = c;
		//edges[2].p1 = c, edges[2].p2 = a;

		//this->edges.push_back(&e1);
		//this->edges.push_back(&e2);
		//this->edges.push_back(&e3);
	}

	Triangle& operator=(const Triangle&) = default;
	Triangle& operator=(Triangle&&) = default;
	bool operator==(Triangle t) {
		return t.a == a && t.b == b && t.c == c;
	}

	bool operator > (const Triangle & t) const{

		return this->a >= t.a;
	}
	bool operator < (const Triangle & t) const{

		return this->a < t.a;
	}

	~Triangle()
	{
		//cout << this->a << " " << this->b << " " << this->c << "destroyed! " << endl;
	}
};

#endif // !TRIANGLE_H


