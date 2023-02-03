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

	map<int,Skiplist> ht;				//��������
	double ht_w;						//�������
	int ht_m;							//��������±꣬�±�Ϊ0~ht_m
	
	double min_y, min_x, max_x, max_y;	//�㼯�еı߽�����	
	const int len = 600;

	vector<Triangle> _triangles;		//��������Ƭ����
	set<Edge> _edges;					//�߼���
	vector<Point> _points;

	//ӳ���ϵ����
	map<Edge, Triangle> LT;				//��->����ߵ�������


public:

	double DistTwoPoints(int qid, int pid) {
		if (pid == -1 || qid == -1)
			return 1e11;
		return sqrt(pow(_points[pid].x - _points[qid].x, 2) + pow(_points[pid].y - _points[qid].y, 2));
	}

	double JudgePointLine_Cross(Point p1, Point p2, Point p) {
		double temp = (p2.x - p1.x) * (p.y - p1.y) - (p.x - p1.x) * (p2.y - p1.y);		//ʹ������������б���

		return temp;
	}

	//��������ʽ
	double Det3(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3) {
		return a1 * (b2 * c3 - c2 * b3) + a2 * (b3 * c1 - b1 * c3) + a3 * (b1 * c2 - c1 * b2);
	}

	//P�Ƿ������������Բ�ڵ��б�׼��,true��ʾԲ��
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

		const Point circumCenter(rx, ry, 0);		//Բ��
		const double circumRadius = _points[Tri.a].Distance2(circumCenter);		//a�㵽Բ�ĵľ���
		const double dist = P.Distance2(circumCenter);	//p�㵽Բ�ĵľ���

		//cout << "dist:  " << sqrt(dist) << "   R:  " << sqrt(circumRadius) << "   Բ�ģ� " << circumCenter.x << "  " << circumCenter.y << "  Point: " << P.Num << "  Tri: " << Tri.a << " " << Tri.b << " " << Tri.c << endl;

		return dist < circumRadius;

	}


	bool PointInTri_Cross(Point p, Triangle tri) {
		Point p1 = _points[tri.a];
		Point p2 = _points[tri.b];
		Point p3 = _points[tri.c];

		int flag = 0;

		if (JudgePointLine_Cross(p1, p2, p) <= 0) {	//λ�������ߵ��Ҳ�
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

	//һ����ht��_points�в����ʼ���ε��ĸ����㣬����ӳ���ϵ��ʼ��
	void Init(const vector<Point>& points) {

		 max_x = 0, max_y = 0;
		 min_y = 1e10, min_x = 1e10;

		//ȷ����ֵ
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
		double Threshold = 5.5 * sqrt(n);	//����ʽ�㷨ȷ������ֵ
		ht_m = n / Threshold;				//��������
		ht_w = (max_y-min_y) / ht_m;		//�������

		//�����ĸ��¶��㼰���������Σ���push��ȥ
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
		//�ĸ�����push��ht��
		for (auto p : _points){
			int k = (p.y - min_y) / ht_w;
			cout << "k: " << k << endl;
			ht[k].add(p.x, p.Num);
		//	this->Print();
		}
		

		//����LTӳ���ϵ��������init_tri																		
		LT.insert(PLT(Edge(0, 2), tri1));
		LT.insert(PLT(Edge(2, 0), tri2));

		//����STӳ���ϵ
		_points[0].VT.push_back(tri1), _points[0].VT.push_back(tri2);
		_points[1].VT.push_back(tri1);
		_points[2].VT.push_back(tri1), _points[2].VT.push_back(tri2);
		_points[3].VT.push_back(tri2);

		//����SPӳ���ϵ
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

	//�ҵ��ڱ�����ht[k]�о���p����ĵ㣬�����ظõ������
	int NearestPointStrip(int pid, int k) {	
		int cnt = _points[pid].Num;
		//int k = (_points[pid].y - min_y) / ht_w;
		ht[k].add(_points[pid].x, cnt);
		double dmin = 1e9;
		int id_pre = -1, id_cur = -1;
		SkiplistNode* tmp = ht[k].search(_points[pid].x);	//�ҵ�pid���ڵ�λ��
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

	
	//ȫ��Χ��Ѱ�Ҿ���pid����ĵ�,����pid��
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
		vector<Triangle> result;			//������㹲ͬ���ڵ�����������
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

	
	
	//����p�������ζ��㣬���������µ������Σ��Ƿ������������ eg�����ڱ���
	void BuildNewTri(int pid, Triangle tri) {

		//����LT�����������µ�������
		Triangle tri1(_points[pid].Num, tri.a, tri.b);
		Triangle tri2(_points[pid].Num, tri.b, tri.c);
		Triangle tri3(_points[pid].Num, tri.c, tri.a);
		_triangles.push_back(tri1);
		_triangles.push_back(tri2);
		_triangles.push_back(tri3);

		//ɾ���������β�ɾ����֮��ص�LTӳ��
		vector<Triangle>::iterator iter = find(_triangles.begin(), _triangles.end(), tri);
		_triangles.erase(iter);
		LT.erase(Edge(tri.b, tri.a));
		LT.erase(Edge(tri.c, tri.b));
		LT.erase(Edge(tri.a, tri.c));

		//����LTӳ���ϵ,��Ȧ���ԣ��������		
		LT[Edge(tri.b, tri.a)] = tri1;
		LT[Edge(tri.c, tri.b)] = tri2;
		LT[Edge(tri.a, tri.c)] = tri3;
		LT[Edge(pid, tri.a)] = tri3, LT[Edge(tri.a, pid)] = tri1;
		LT[Edge(tri.b, pid)] = tri2, LT[Edge(pid, tri.b)] = tri1;
		LT[Edge(tri.c, pid)] = tri3, LT[Edge(pid, tri.c)] = tri2;

		

		//ɾ��ST���ɵĶ���-������ӳ���ϵ
		vector<Triangle>::iterator iter1 = find(_points[tri.a].VT.begin(), _points[tri.a].VT.end(), tri);
		_points[tri.a].VT.erase(iter1);
		vector<Triangle>::iterator iter2 = find(_points[tri.b].VT.begin(), _points[tri.b].VT.end(), tri);
		_points[tri.b].VT.erase(iter2);
		vector<Triangle>::iterator iter3 = find(_points[tri.c].VT.begin(), _points[tri.c].VT.end(), tri);
		_points[tri.c].VT.erase(iter3);

		
		//����ST������-������ӳ���ϵ
		_points[pid].VT.push_back(tri1);				//��pΪ����
		_points[pid].VT.push_back(tri2);
		_points[pid].VT.push_back(tri3);
		_points[tri.a].VT.push_back(tri1);				//��a��b��cΪ����
		_points[tri.a].VT.push_back(tri3);
		_points[tri.b].VT.push_back(tri1);
		_points[tri.b].VT.push_back(tri2);
		_points[tri.c].VT.push_back(tri2);
		_points[tri.c].VT.push_back(tri3);

		

		//����SP��������ӹ�ϵ
		_points[pid].VP.push_back(_points[tri.a]);		//��pΪ����
		_points[pid].VP.push_back(_points[tri.b]);
		_points[pid].VP.push_back(_points[tri.c]);
		_points[tri.a].VP.push_back(_points[pid]);		//��a��b��cΪ����
		_points[tri.b].VP.push_back(_points[pid]);
		_points[tri.c].VP.push_back(_points[pid]);

		
	}

	void CreatePoly(Point p, vector<Edge>& Poly) {

		for (int i = 0; i < Poly.size(); i++) {
			if (LT.find(Poly[i]) != LT.end()) {			//������������ڽ������Σ��ſ����Ƿ��������Բ��

				Triangle tri = LT.at(Poly[i]);		//tri��ʾ�������ⲿ���ڽ�������

				/***************************************************************************************�������޸Ŀ����Բ����ж�***************/
				if (JudgeCricumcircle_Det3(tri, p) == false) {	//pλ��Բ��,��Ҫ����ָ���������ڵ��ǶԼ�ֵ

					if (LT.find(Poly[i]) != LT.end()) {		//����������ⲿ����������δ����map�Ļ���map�в�����Լ�ֵ

						//LT.insert(PLT(Poly[i], Triangle(p.Num, Poly[i].p2, Poly[i].p1)));
						LT[Edge(Poly[i].p2, Poly[i].p1)] = Triangle(p.Num, Poly[i].p1, Poly[i].p2);

					}

					continue;
				}
				else {												//pλ��Բ��

					Triangle tt(tri.a, tri.b, tri.c);

					//��Ҫɾ������ڽ������Σ�����map��ȥ����Լ�ֵ(ע���ֵ�Ե�ɾ������������ε�ÿ���ߵ����Լ�ֵ����Ҫɾ����Ϊ��ȥ�������ߵ��ڽ������ζ�Ӧ��ɾ��
					vector<Triangle>::iterator iter = find(_triangles.begin(), _triangles.end(), tt);
					_triangles.erase(iter);

					LT.erase(Poly[i]);
					LT.erase(Edge(Poly[i].p2, Poly[i].p1));
					int pp1 = Poly[i].p1, pp2 = Poly[i].p2;
					//�Ȳ����±ߣ���ɾ��ԭ�ߣ�ע��˳��
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
					//ʹ����ѭ������ǲ�����±�						
					i--;
				}
			}
			else {					//���ڽ������Σ����ù�
				continue;
			}
		}

		vector<Edge> Polyf;
		Polyf.push_back(Poly[0]);

		//��Poly�Ÿ���ʹPoly��β���
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

		cout << "�õ��Ӱ����Poly������ɣ�" << endl;
		for (int i = 0; i < Poly.size(); i++) {
			cout << "    " << Poly[i].p1 << "     " << Poly[i].p2 << endl;
		}

	}

	//��ѯab���ұߵĵ�
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

	//�ֲ��Ż���ͨ�������Խ��ߣ��Ƿ�����ʹ�õݹ飿
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
			if (pr_id != -1 && pl_id != -1) {	//����������ͬʱ����ʱ���ſ��Խ����Խ���
				Point* pr =new Point(_points[pr_id]);				//���ұߵĵ�
				Point* pl =new Point(_points[pl_id]);				//����ߵĵ�
				if (JudgeCricumcircle_Det3(LT[e], *pr) == true) {	//���p�����Բ�ڵĻ����������ǵ��������߼���ջ��

					//���µ�����push��poly,ע��˳��
					int oid = PointId(e.p2, e.p1);
					Poly.push(Edge(oid, e.p2));
					Poly.push(Edge(e.p1, oid));

					//ɾ��ԭ�������������Σ�ɾ����������������ص�LT�����������VT��VP
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

					//���������µ������Σ����������ߵ�LT���ĸ������VT��VP		(��BuildNewTri���ò�ͬ��
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

		//ȷ���¼����������֮����ڽӹ�ϵ		
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
		//����ӳ���ϵ
	}
	

	void Triangulation(const vector<Point>& points) {

		Init(points);

		int cnt = 4;

		//_points.resize(10000);
		//_triangles.resize(10000);

		for (auto p : points) {
			
			//��points�еĵ������ӵ�_points��
			Point pp(p.x, p.y, p.z, cnt);
			int pid = pp.Num;
			_points.push_back(pp);

			//�ҵ�����pi����ĵ�pk
			int kid = FindNearestPoint(pid);

			//�󽫵�ǰ������������ڳ�pid��ĵ��γɵ��������ҵ�����pid����ĵ�
			int k = (p.y - min_y) / ht_w;
			ht[k].add(p.x, cnt);
			cnt++;
			this->Print();


			//����pk���������е��У��ҵ���pk,pi) �� (pk,pq)�������q��
			int qid = FindOtherPoint(pid, kid);

			//ͨ����ѯpk��pq���������εĽ�����ȷ��pi���ڵ�������
			Triangle tri = PointInTri(pid, kid, qid);

			//���������µ������Σ�����ӳ���ϵ
			BuildNewTri(pid, tri);

			//�ֲ��Ż���ʹ�����������Բ׼��
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