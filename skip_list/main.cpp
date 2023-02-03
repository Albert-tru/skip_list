#include <malloc.h>
#include <vld.h>

#include"Skiplist.h"
#include"Edge.h"
#include"Triangle.h"
#include"Point.h"
#include"Delaunay_skiplist.h"

#include <fstream>
#include <sstream>
#include<algorithm>


#define M_TRIM_THRESHOLD    -1
#define M_TOP_PAD           -2
#define M_MMAP_THRESHOLD    -3
#define M_MMAP_MAX          -4
#define M_CHECK_ACTION      -5
#define M_PERTURB           -6
#define M_ARENA_TEST        -7
#define M_ARENA_MAX         -8

#ifndef M_MXFAST
# define M_MXFAST  1    /* maximum request size for "fastbins" */
#endif
#ifndef M_NLBLKS
# define M_NLBLKS  2    /* UNUSED in this malloc */
#endif
#ifndef M_GRAIN
# define M_GRAIN   3    /* UNUSED in this malloc */
#endif
#ifndef M_KEEP
# define M_KEEP    4    /* UNUSED in this malloc */
#endif



double arr1[4001][3];
vector<Point>points;


void ReadFile() {
	ifstream infile;
	infile.open("200组随机数.txt");
	if (!infile) cout << "error" << endl;
	double t1, t2, t3, t4, t5, t6, t7;
	int i = 0;
	int j = 0;
	while (infile >> t1 >> t2) {
		cout << t1 << "         " << t2 << endl;
		arr1[i][0] = t1;
		arr1[i][1] = t2;
		arr1[i][2] = 0;
		j++;
		//cout << arr1[i][0] << " " << arr1[i][1] << " " << arr1[i][2] << " " << j << endl;
		i++;
	}
	int N = 4000;

	double centx = 50, centy = 50;
	for (int i = 0; i < j; i++)
	{
		Point tp; tp.Num = i;
		//float x = (rand() % (999 + 1) / (float)(999 + 1)) * 100;	//产生（1,100）之间的随机数，保留三位小数
		//float y = (rand() % (999 + 1) / (float)(999 + 1)) * 100;	
		//arr1[i][0] = x, arr1[i][1] = y;
		//cout << arr1[i][0] << " " << arr1[i][1] << " " << arr1[i][2] << " " << endl;
		tp.x = arr1[i][0], tp.y = arr1[i][1], tp.z = 0;
		double dis = sqrt(pow(tp.x - centx, 2) + pow(tp.y - centy, 2));
		//cout << "iiiii  " << i << endl;
		tp.dist = dis;
		points.insert(points.end(), tp);
		//ps1->InsertNextPoint(arr1[i][0], arr1[i][1], arr1[i][2]);
	}

}

void RandNumber() {
	//范围： 
	int Num = 100;
	for (int i = 0; i < Num; ++i)
	{
		//temp.push_back(i + 1);
	}

	//random_shuffle(temp.begin(), temp.end());
}

void TestAdd() {
	

	//sp.search(0);

}

int main() {

	// 禁止malloc调用mmap分配内存
	//mallopt(M_TRIM_THRESHOLD, 0); // 禁止内存缩进，sbrk申请的内存释放后不会归还给操作系统


	Skiplist sp;

	//sp.add(12.1,3);

	//ReadFile();


	Delaunay_skiplist dsl;

	//dsl.DataPreparation(points);
	//dsl.Init();
	dsl.Triangulation(points);

	//dsl.Print();

	//dsl.FindNearestPoint(0);

	for (auto t : points) {
		sp.erase(t.x, t.Num);
	}
	//delete dsl;
	
	return 0;
}
