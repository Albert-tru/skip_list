#pragma once
#include <vld.h>
#include"Point.h"
#include"Edge.h"
#include"Triangle.h"
#include"Point.h"
#include"skiplist.h"

#include<vector>
#include<algorithm>
#include<iostream>
#include<map>
#include<set>
#include<cmath>
#include<set>
#include<stack>

using namespace std;

typedef pair<Edge, Triangle> PLT;



class Delaunay_skiplist {

	map<int,Skiplist> ht;				//条带跳表
	double ht_w;						//条带宽度
	int ht_m;							//条带最大下标，下表为0~ht_m
	
	double min_y, min_x, max_x, max_y;	//点集中的边界坐标	
	const int len = 600;

	vector<Triangle> _triangles;		//三角形面片集合
	set<Edge> _edges;					//边集合
	vector<Point> _points;

	//映射关系建立
	map<Edge, Triangle> LT;				//边->边左边的三角形


public:

	double DistTwoPoints(int qid, int pid) {
		if (pid == -1 || qid == -1)
			return 1e11;
		return sqrt(pow(_points[pid].x - _points[qid].x, 2) + pow(_points[pid].y - _points[qid].y, 2));
	}

	double JudgePointLine_Cross(Point p1, Point p2, Point p) {
		double temp = (p2.x - p1.x) * (p.y - p1.y) - (p.x - p1.x) * (p2.y - p1.y);		//使用向量叉乘来判别方向

		return temp;
	}

	//三阶行列式
	double Det3(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3) {
		return a1 * (b2 * c3 - c2 * b3) + a2 * (b3 * c1 - b1 * c3) + a3 * (b1 * c2 - c1 * b2);
	}

	//P是否在三角形外接圆内的判别准则,true表示圆内
	bool JudgeCricumcircle_Det3(Triangle Tri, Point P) {
		double ax = _points[Tri.a].x;
		double ay = _points[Tri.a].y;
		double aa = pow(ax, 2) + pow(ay, 2);
		double bx = _points[Tri.b].x;
		double by = _points[Tri.b].y;
		double bb = pow(bx, 2) + pow(by, 2);
		double cx = _points[Tri.c].x;
		double cy = _points[Tri.c].y;
		double cc = pow(cx, 2) + pow(cy, 2);
		double px = P.x;
		double py = P.y;
		double pp = pow(px, 2) + pow(py, 2);
		//double temp = ax * Det3(by, bb, 1, cy, cc, 1, px, pp, 1) - ay * Det3(bx, bb, 1, cx, cc, 1, px, pp, 1)
		//	+ aa * Det3(bx, by, 1, cx, cy, 1, px, py, 1) - Det3(bx, by, bb, cx, cy, cc, px, py, pp);

		//cout << ax << " " << ay << " " << aa << " " << bx << " " << by << " "
		//	<< bb << " " << cx << " " << cy << " " << cc << " " << px << " " << py << " " << pp << " " << temp << endl;

		double rx = Det3(aa, ay, 1, bb, by, 1, cc, cy, 1) / (2.0 * Det3(ax, ay, 1, bx, by, 1, cx, cy, 1));
		double ry = Det3(ax, aa, 1, bx, bb, 1, cx, cc, 1) / (2.0 * Det3(ax, ay, 1, bx, by, 1, cx, cy, 1));

		const Point circumCenter(rx, ry, 0);		//圆心
		const double circumRadius = _points[Tri.a].Distance2(circumCenter);		//a点到圆心的距离
		const double dist = P.Distance2(circumCenter);	//p点到圆心的距离

		//cout << "dist:  " << sqrt(dist) << "   R:  " << sqrt(circumRadius) << "   圆心： " << circumCenter.x << "  " << circumCenter.y << "  Point: " << P.Num << "  Tri: " << Tri.a << " " << Tri.b << " " << Tri.c << endl;

		return dist < circumRadius;

	}


	bool PointInTri_Cross(Point p, Triangle tri) {
		Point p1 = _points[tri.a];
		Point p2 = _points[tri.b];
		Point p3 = _points[tri.c];

		int flag = 0;

		if (JudgePointLine_Cross(p1, p2, p) <= 0) {	//位于这条边的右侧
			flag++;
		}
		if (JudgePointLine_Cross(p2, p3, p) <= 0) {
			flag++;
		}
		if (JudgePointLine_Cross(p3, p1, p) <= 0) {
			flag++;
		}
		if (flag == 3)
			return true;
		else
			return false;
	}

	//一、在ht、_points中插入初始矩形的四个顶点，并对映射关系初始化
	void Init(const vector<Point>& points) {

		 max_x = 0, max_y = 0;
		 min_y = 1e10, min_x = 1e10;

		//确定最值
		for (auto t : points) {
			//_points.push_back(t);
			if (t.x > max_x)
				max_x = t.x;
			if (t.y > max_y)
				max_y = t.y;
			if (t.y < min_y) {
				min_y = t.y;
			}
			if (t.x < min_x) {
				min_x = t.x;
			}
		}

		int n = points.size();
		double Threshold = 5.5 * sqrt(n);	//启发式算法确定的阈值
		ht_m = n / Threshold;				//条带数量
		ht_w = (max_y-min_y) / ht_m;		//条带宽度

		//建立四个新顶点及两个三角形，并push进去
		Point *p1 = new Point(0.9 * min_x, min_y, 0, 0);
		Point p2(0.8 * min_x, max_y, 0, 1);
		Point p3(1.1 * max_x, max_y, 0, 2);
		Point p4(1.2 * max_x, min_y, 0, 3);
		Triangle tri1(0, 1, 2);
		Triangle tri2(2, 3, 0);
		_triangles.push_back(tri1);
		_triangles.push_back(tri2);
		_points.push_back(*p1);
		_points.push_back(p2);
		_points.push_back(p3);
		_points.push_back(p4);

		vector<int>v = { 3,1,2,4 };
		//四个顶点push到ht中
		for (auto p : _points){
			int k = (p.y - min_y) / ht_w;
			cout << "k: " << k << endl;
			ht[k].add(p.x, p.Num);
		//	this->Print();
		}
		

		//建立LT映射关系，三边与init_tri																		
		LT.insert(PLT(Edge(0, 2), tri1));
		LT.insert(PLT(Edge(2, 0), tri2));

		//建立ST映射关系
		_points[0].VT.push_back(tri1), _points[0].VT.push_back(tri2);
		_points[1].VT.push_back(tri1);
		_points[2].VT.push_back(tri1), _points[2].VT.push_back(tri2);
		_points[3].VT.push_back(tri2);

		//建立SP映射关系
		_points[0].VP.push_back(_points[1]);
		_points[0].VP.push_back(_points[2]);
		_points[0].VP.push_back(_points[3]);
		_points[1].VP.push_back(_points[0]);
		_points[1].VP.push_back(_points[2]);
		_points[2].VP.push_back(_points[1]);
		_points[2].VP.push_back(_points[3]);
		_points[3].VP.push_back(_points[0]);
		_points[3].VP.push_back(_points[2]);
		_points[2].VP.push_back(_points[0]);
		
		delete p1;
	}

	//找到在本条带ht[k]中距离p最近的点，并返回该点的索引
	int NearestPointStrip(int pid, int k) {	
		int cnt = _points[pid].Num;
		//int k = (_points[pid].y - min_y) / ht_w;
		ht[k].add(_points[pid].x, cnt);
		double dmin = 1e9;
		int id_pre = -1, id_cur = -1;
		SkiplistNode* tmp = ht[k].search(_points[pid].x);	//找到pid所在的位置
		SkiplistNode* pre = tmp->_frontV[0];
		SkiplistNode* cur = tmp->_nextV[0];
		int id = -1;
		double dx_left = 1e10, dx_right = 1e10;
		double d_left = 1e10, d_right = 1e10;
		if (pre && pre->val != -1)
			dx_left = tmp->val - pre->val, d_left = DistTwoPoints(pid, pre->id);
		if (cur && cur->val != -1)
			dx_right = cur->val - tmp->val, d_right = DistTwoPoints(pid, cur->id);
		if (d_right < dmin)
			dmin = d_right, id = cur->id;
		if (d_left < dmin)
			dmin = d_left, id = pre->id;
		while ((dmin > dx_left || dmin > dx_right ) && cur && pre) {
			if (pre->id != -1 && pre != nullptr) {
				d_left = DistTwoPoints(pid, pre->id);
				if (d_left < dmin)
					dmin = d_left, id = pre->id;
			}
			if (cur->id != -1 && cur != nullptr) {
				d_right = DistTwoPoints(pid, cur->id);
				if (d_right < dmin)
					dmin = d_right, id = cur->id;
			}
			if (cur && cur->_nextV[0] != nullptr) {
				cur = cur->_nextV[0];
				dx_right = cur->val - tmp->val;
			}
			else
				break;
			if (pre && pre->_frontV[0] != nullptr) {
				pre = pre->_frontV[0];
				dx_left = tmp->val - pre->val;
			}
			else
				break;
		}
		cout << "***********************************************************************************************" << endl;
		cout << "before: " << endl;
		this->Print();
		cout << "-------------------------------------------------------------" << endl;
		ht[k].erase(_points[pid].x, pid);
		cout << "erase: " << _points[pid].x << endl;
		this->Print();
		cout << "***********************************************************************************************" << endl;
		return id;
	}

	
	//全域范围内寻找距离pid最近的点,而非pid点
	int FindNearestPoint(int pid) {
		int k = (_points[pid].y - min_y) / ht_w;
		int kid = NearestPointStrip(pid, k);
		double dmin = 1e10;
		if(kid != -1)
			dmin = DistTwoPoints(pid, kid);

		bool moveUp = true, moveDown = true;
		int upper = k - 1, lower = k + 1;
		while (moveUp || moveDown) {
			double d_upper = _points[pid].y - k * ht_w;
			double d_lower = (k + 1) * ht_w - _points[pid].y;
			int tid = -1;
			if (moveUp) {
				if (d_upper < dmin && upper >= 0) {
					tid = NearestPointStrip(pid, upper);
					if (tid != -1 && DistTwoPoints(tid, pid) < DistTwoPoints(kid, pid)) {
						kid = tid;
					}
					upper -= 1;
				}
				else {
					moveUp = false;
				}
			}
			if (moveDown) {
				if (d_lower < dmin && lower <= ht_m) {
					tid = NearestPointStrip(pid, lower);
					if (tid != -1 && DistTwoPoints(tid, pid) < DistTwoPoints(kid, pid)) {
						kid = tid;
					}
					lower += 1;
				}
				else {
					moveDown = false;
				}
			}
		}
		cout << "p: " << _points[pid].x << " " << _points[pid].y << endl;
		cout << "q: " << _points[kid].x << " " << _points[kid].y << endl;
		return kid;
	}


	int FindOtherPoint(int pid, int kid) {
		int qid = -1;
		double ans = -1e10;
		for (auto p : _points[kid].VP) {
			double tmp = (_points[pid].x - _points[kid].x) * (p.x - _points[kid].x) +
				(_points[pid].y - _points[kid].y) * (p.y - _points[kid].y);
			if (tmp > ans) {
				ans = tmp;
				qid = p.Num;
			}
		}
		
		return qid;
	}

	
	Triangle PointInTri(int pid, int kid, int qid) {
		vector<Triangle> result;			//存放两点共同存在的两个三角形
		for (auto tri : _points[kid].VT) {
			vector<Triangle>::iterator iter = find(_points[qid].VT.begin(), _points[qid].VT.end(), tri);
			if (iter != _points[qid].VT.end()) {
				if (PointInTri_Cross(_points[pid], tri) == true) {
					return tri;
				}
			}
		}


		
		/*
		for (auto tri : result) {
			if (PointInTri_Cross(_points[pid], tri) == true) {
				return tri;
			}
		}
		*/
	}

	
	
	//连接p与三角形顶点，生成三个新的三角形，是否有特殊情况？ eg：点在边上
	void BuildNewTri(int pid, Triangle tri) {

		//更新LT，连接生成新的三角形
		Triangle tri1(_points[pid].Num, tri.a, tri.b);
		Triangle tri2(_points[pid].Num, tri.b, tri.c);
		Triangle tri3(_points[pid].Num, tri.c, tri.a);
		_triangles.push_back(tri1);
		_triangles.push_back(tri2);
		_triangles.push_back(tri3);

		//删除旧三角形并删除与之相关的LT映射
		vector<Triangle>::iterator iter = find(_triangles.begin(), _triangles.end(), tri);
		_triangles.erase(iter);
		LT.erase(Edge(tri.b, tri.a));
		LT.erase(Edge(tri.c, tri.b));
		LT.erase(Edge(tri.a, tri.c));

		//更新LT映射关系,外圈三对，侧边六对		
		LT[Edge(tri.b, tri.a)] = tri1;
		LT[Edge(tri.c, tri.b)] = tri2;
		LT[Edge(tri.a, tri.c)] = tri3;
		LT[Edge(pid, tri.a)] = tri3, LT[Edge(tri.a, pid)] = tri1;
		LT[Edge(tri.b, pid)] = tri2, LT[Edge(pid, tri.b)] = tri1;
		LT[Edge(tri.c, pid)] = tri3, LT[Edge(pid, tri.c)] = tri2;

		

		//删除ST，旧的顶点-三角形映射关系
		vector<Triangle>::iterator iter1 = find(_points[tri.a].VT.begin(), _points[tri.a].VT.end(), tri);
		_points[tri.a].VT.erase(iter1);
		vector<Triangle>::iterator iter2 = find(_points[tri.b].VT.begin(), _points[tri.b].VT.end(), tri);
		_points[tri.b].VT.erase(iter2);
		vector<Triangle>::iterator iter3 = find(_points[tri.c].VT.begin(), _points[tri.c].VT.end(), tri);
		_points[tri.c].VT.erase(iter3);

		
		//更新ST，顶点-三角形映射关系
		_points[pid].VT.push_back(tri1);				//以p为顶点
		_points[pid].VT.push_back(tri2);
		_points[pid].VT.push_back(tri3);
		_points[tri.a].VT.push_back(tri1);				//以a、b、c为顶点
		_points[tri.a].VT.push_back(tri3);
		_points[tri.b].VT.push_back(tri1);
		_points[tri.b].VT.push_back(tri2);
		_points[tri.c].VT.push_back(tri2);
		_points[tri.c].VT.push_back(tri3);

		

		//更新SP，点的连接关系
		_points[pid].VP.push_back(_points[tri.a]);		//以p为顶点
		_points[pid].VP.push_back(_points[tri.b]);
		_points[pid].VP.push_back(_points[tri.c]);
		_points[tri.a].VP.push_back(_points[pid]);		//以a、b、c为顶点
		_points[tri.b].VP.push_back(_points[pid]);
		_points[tri.c].VP.push_back(_points[pid]);

		
	}

	void CreatePoly(Point p, vector<Edge>& Poly) {

		for (int i = 0; i < Poly.size(); i++) {
			if (LT.find(Poly[i]) != LT.end()) {			//如果这条边有邻接三角形，才考虑是否在其外接圆内

				Triangle tri = LT.at(Poly[i]);		//tri表示这条边外部的邻接三角形

				/***************************************************************************************在这里修改空外接圆检测判断***************/
				if (JudgeCricumcircle_Det3(tri, p) == false) {	//p位于圆外,需要覆盖指向这条边内的那对键值

					if (LT.find(Poly[i]) != LT.end()) {		//如果这条边外部有三角形且未加入map的话，map中插入这对键值

						//LT.insert(PLT(Poly[i], Triangle(p.Num, Poly[i].p2, Poly[i].p1)));
						LT[Edge(Poly[i].p2, Poly[i].p1)] = Triangle(p.Num, Poly[i].p1, Poly[i].p2);

					}

					continue;
				}
				else {												//p位于圆内

					Triangle tt(tri.a, tri.b, tri.c);

					//需要删除这个邻接三角形，并在map中去除这对键值(注意键值对的删除，这个三角形的每条边的两对键值都需要删，因为边去掉后，两边的邻接三角形都应该删除
					vector<Triangle>::iterator iter = find(_triangles.begin(), _triangles.end(), tt);
					_triangles.erase(iter);

					LT.erase(Poly[i]);
					LT.erase(Edge(Poly[i].p2, Poly[i].p1));
					int pp1 = Poly[i].p1, pp2 = Poly[i].p2;
					//先插入新边，再删除原边，注意顺序
					if (pp1 == tri.a) {
						Poly.insert(Poly.begin() + i, Edge(tri.a, tri.b));
						Poly.insert(Poly.begin() + i + 1, Edge(tri.b, tri.c));
						for (int j = 0; j < Poly.size(); j++) {
							if (Poly[j].p1 == pp1 && Poly[j].p2 == pp2) {
								Poly.erase(Poly.begin() + j);
							}
						}
					}
					else if (pp1 == tri.b) {
						Poly.insert(Poly.begin() + i, Edge(tri.c, tri.a));
						Poly.insert(Poly.begin() + i + 1, Edge(tri.b, tri.c));
						for (int j = 0; j < Poly.size(); j++) {
							if (Poly[j].p1 == pp1 && Poly[j].p2 == pp2) {
								Poly.erase(Poly.begin() + j);
							}
						}
					}
					else if (pp1 == tri.c) {
						Poly.insert(Poly.begin() + i, Edge(tri.c, tri.a));
						Poly.insert(Poly.begin() + i + 1, Edge(tri.a, tri.b));
						for (int j = 0; j < Poly.size(); j++) {
							if (Poly[j].p1 == pp1 && Poly[j].p2 == pp2) {
								Poly.erase(Poly.begin() + j);
							}
						}
					}
					//使下轮循环求的是插入的新边						
					i--;
				}
			}
			else {					//无邻接三角形，不用管
				continue;
			}
		}

		vector<Edge> Polyf;
		Polyf.push_back(Poly[0]);

		//给Poly排个序，使Poly首尾相接
		for (int i = 1; i < Poly.size(); i++) {
			int t = Polyf[i - 1].p2;
			int j, id1, id2;
			for (j = 0; j < Poly.size(); j++) {
				if (Poly[j].p1 == t)
				{
					t = j;
					break;
				}

			}
			Polyf.push_back(Poly[t]);
		}
		for (int i = 0; i < Polyf.size(); i++) {
			Poly[i].p1 = Polyf[i].p1;
			Poly[i].p2 = Polyf[i].p2;
		}

		cout << "该点的影响域Poly创建完成！" << endl;
		for (int i = 0; i < Poly.size(); i++) {
			cout << "    " << Poly[i].p1 << "     " << Poly[i].p2 << endl;
		}

	}

	//查询ab边右边的点
	int PointId(int a, int b) {	
		int p = -1;
		Triangle tri = LT[Edge(b, a)];
		if (tri.a == a && tri.b == b)
			p = tri.c;
		if (tri.b == a && tri.c == b)
			p = tri.a;
		if (tri.c == a && tri.a == b)
			p = tri.b;
		return p;
	}

	//局部优化，通过交换对角线，是否必须得使用递归？
	void  LocalOptimization(Triangle tri) {
		
		stack<Edge> Poly;
		Poly.push(Edge(tri.c, tri.a));
		Poly.push(Edge(tri.b, tri.c));
		Poly.push(Edge(tri.a, tri.b));
		

		while (!Poly.empty()) {
			Edge e = Poly.top();
			Poly.pop();
			int pr_id = PointId(e.p1, e.p2);
			int pl_id = PointId(e.p2, e.p1);
			if (pr_id != -1 && pl_id != -1) {	//必须这两点同时存在时，才可以交换对角线
				Point* pr =new Point(_points[pr_id]);				//边右边的点
				Point* pl =new Point(_points[pl_id]);				//边左边的点
				if (JudgeCricumcircle_Det3(LT[e], *pr) == true) {	//如果p在外接圆内的话，把外三角的其余两边加入栈中

					//把新的两边push进poly,注意顺序
					int oid = PointId(e.p2, e.p1);
					Poly.push(Edge(oid, e.p2));
					Poly.push(Edge(e.p1, oid));

					//删除原来的两个三角形，删除这两个三角形相关的LT，两个顶点的VT、VP
					Triangle tri1 = LT[e];
					Triangle tri2 = LT[Edge(e.p2, e.p1)];
					vector<Triangle>::iterator iter1 = remove(_triangles.begin(), _triangles.end(), tri1);
					_triangles.pop_back();
					vector<Triangle>::iterator iter2 = remove(_triangles.begin(), _triangles.end(), tri2);
					
					//_triangles.erase(iter2,_triangles.end());
					_triangles.pop_back();

					LT.erase(Edge(tri1.b, tri1.a));
					LT.erase(Edge(tri1.c, tri1.b));
					LT.erase(Edge(tri1.a, tri1.c));
					LT.erase(Edge(tri2.b, tri2.a));
					LT.erase(Edge(tri2.c, tri2.b));
					LT.erase(Edge(tri2.a, tri2.c));

					vector<Triangle>::iterator iter_t1 = remove(_points[tri1.a].VT.begin(), _points[tri1.a].VT.end(), tri1);
					_points[tri1.a].VT.pop_back();
					vector<Triangle>::iterator iter_t2 = remove(_points[tri1.b].VT.begin(), _points[tri1.b].VT.end(), tri1);
					_points[tri1.b].VT.pop_back();
					vector<Triangle>::iterator iter_t3 = remove(_points[tri1.c].VT.begin(), _points[tri1.c].VT.end(), tri1);
					_points[tri1.c].VT.pop_back();

					vector<Triangle>::iterator iter_tt1 = remove(_points[tri2.a].VT.begin(), _points[tri2.a].VT.end(), tri2);
					_points[tri2.a].VT.pop_back();
					vector<Triangle>::iterator iter_tt2 = remove(_points[tri2.b].VT.begin(), _points[tri2.b].VT.end(), tri2);
					_points[tri2.b].VT.pop_back();
					vector<Triangle>::iterator iter_tt3 = remove(_points[tri2.c].VT.begin(), _points[tri2.c].VT.end(), tri2);
					_points[tri2.c].VT.pop_back();

					vector<Point>::iterator iter_p1 = remove(_points[e.p1].VP.begin(), _points[e.p1].VP.end(), _points[e.p2]);
					vector<Point>::iterator iter_p2 = remove(_points[e.p2].VP.begin(), _points[e.p2].VP.end(), _points[e.p1]);
					_points[e.p1].VP.pop_back();
					_points[e.p2].VP.pop_back();

					//建立两个新的三角形，建立五条边的LT，四个顶点的VT、VP		(和BuildNewTri作用不同）
					Triangle* new_tri1 = new Triangle(pl->Num, e.p2, pr->Num);
					Triangle* new_tri2 = new Triangle(pl->Num, pr->Num, e.p1);
					_triangles.push_back(*new_tri1);
					_triangles.push_back(*new_tri2);

					LT[Edge(pl->Num, pr->Num)] = *new_tri1;
					LT[Edge(pr->Num, pl->Num)] = *new_tri2;
					LT[Edge(pr->Num, e.p2)] = *new_tri1;
					LT[Edge(e.p2, pl->Num)] = *new_tri1;
					LT[Edge(pl->Num, e.p1)] = *new_tri2;
					LT[Edge(e.p1, pr->Num)] = *new_tri2;

					_points[pl->Num].VT.push_back(*new_tri1), _points[pl->Num].VT.push_back(*new_tri2);
					_points[pr->Num].VT.push_back(*new_tri1), _points[pr->Num].VT.push_back(*new_tri2);
					_points[e.p1].VT.push_back(*new_tri2);
					_points[e.p2].VT.push_back(*new_tri1);

					_points[pl->Num].VP.push_back(_points[pr->Num]);
					_points[pr->Num].VP.push_back(_points[pl->Num]);

					delete new_tri1;
					delete new_tri2;
					delete pr;
					delete pl;
				}
				
			}
			


		}

		//确认新加入的三角形之间的邻接关系		
		int tlen = _triangles.size() - 3;
		for (int i = 0; i < _triangles.size() - tlen; i++) {
			if (i < _triangles.size() - tlen - 1) {
				Edge e(_triangles[i].a, _triangles[i].b);
				LT.insert(PLT(e, _triangles[i + 1]));
			}
			if (i >= 1) {
				Edge e1(_triangles[i].c, _triangles[i].a);
				LT.insert(PLT(e1, _triangles[i - 1]));
			}
		}
		int len = _triangles.size();
		Edge e(_triangles[len - tlen - 1].a, _triangles[len - tlen - 1].b);
		LT.insert(PLT(e, _triangles[0]));

		Edge e1(_triangles[0].c, _triangles[0].a);
		LT.insert(PLT(e1, _triangles[len - tlen - 1]));
		//更新映射关系
	}
	

	void Triangulation(const vector<Point>& points) {

		Init(points);

		int cnt = 4;

		//_points.resize(10000);
		//_triangles.resize(10000);

		for (auto p : points) {
			
			//把points中的点逐个添加到_points中
			Point pp(p.x, p.y, p.z, cnt);
			int pid = pp.Num;
			_points.push_back(pp);

			//找到距离pi最近的点pk
			int kid = FindNearestPoint(pid);

			//后将当前点插入跳表，先在除pid外的点形成的跳表中找到距离pid最近的点
			int k = (p.y - min_y) / ht_w;
			ht[k].add(p.x, cnt);
			cnt++;
			this->Print();


			//在与pk相连的所有点中，找到（pk,pi) · (pk,pq)点积最大的q点
			int qid = FindOtherPoint(pid, kid);

			//通过查询pk和pq所在三角形的交集，确定pi所在的三角形
			Triangle tri = PointInTri(pid, kid, qid);

			//建立三个新的三角形，更新映射关系
			BuildNewTri(pid, tri);

			//局部优化，使三角形满足空圆准则
			//LocalOptimization(tri);

			cout << _points.size() << "  " << _triangles.size() << endl;
			cout << endl;

		}
	}

	void Print() {

		cout << "ht_m : " << ht_m << endl;
		cout << "ht_w : " << ht_w << endl;
		for (int k = 0; k <= ht_m; k++) {
			ht[k]._Print();
			cout << endl << endl;
		}

	}

	~Delaunay_skiplist() {
		cout << " Delaunay_skiplist destroyed ! " << endl;
	}

};