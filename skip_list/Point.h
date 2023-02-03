#pragma once

#ifndef POINT
#define POINT

#include"Point.h"
#include"Edge.h"
#include"Triangle.h"
#include"skiplist.h"

#include <math.h>
#include<set>
#include<vector>


class Point
{

public:
	int Num;						//�������
	double x;
	double y;
	double z;
	vector<Triangle>VT;		//����->�Ըõ�������������μ���	
	vector<Point>VP;				//��->��õ������ĵ�

	double dist;					//�����ĵ�ľ���

	Point() = default;				//ʹ��Ĭ�ϵĹ��캯��
	Point(const Point& point) = default;		//ʹ��Ĭ�ϵĿ������캯��

	Point(const double point_x, const double point_y, const double point_z) //һ�㹹�캯��
		:x(point_x), y(point_y), z(point_z)
	{
	}

	Point(const double point_x, const double point_y, const double point_z,int num) //һ�㹹�캯��
		:x(point_x), y(point_y), z(point_z),Num(num)
	{
	}


	double Distance2(const Point& point) const {//
		//����֮������ƽ��
		const double dx = this->x - point.x;
		const double dy = this->y - point.y;
		//const double dz = this->z - point.z;
		//return dx * dx + dy * dy + dz * dz;
		return dx * dx + dy * dy;
	}

	//double Distance(const Point& point)const {

		//return sqrt(Distance2(point));
	//}

	double Norm2() const {
		//�㵽ԭ��ľ���ƽ��
		return x * x + y * y + z * z;
	}



	Point& operator= (const Point&) = default;	//ʹ��Ĭ�ϵĿ������캯��
	Point& operator= (Point&&) = default;
	bool operator == (const Point& point) const {
		//����==���ţ������ж��������Ƿ���ͬ
		return (this->x == point.x) && (this->y == point.y) && (this->z == point.z);
	}
	bool operator > (const Point& point) const {
		
		return this->dist >= point.dist;
	}
	bool operator < (const Point& point) const {

		return this->dist < point.dist;
	}

	~Point()
	{
		//cout << this->Num << " destruct !" << endl;
	}
};



#endif // !POINT


