#pragma once
#include<iostream>
#include<vector>
#include<ctime>
#include <random>

#include"Edge.h"
#include"Triangle.h"
#include"Point.h"


struct SkiplistNode {
    double val;
    int id;
    std::vector<SkiplistNode*>_nextV;           //�ýڵ������ָ��
    std::vector<SkiplistNode*>_frontV;          //�ýڵ������ǰ���ڵ�
    SkiplistNode(double _val, int leve, int _id)
        :val(_val), _nextV(leve, nullptr), _frontV(leve, nullptr), id(_id)
    {}
};

class Skiplist {
    typedef SkiplistNode Node;
private:
    Node* _head;//ͷ�ڵ�
    size_t _maxleve = 32;//ÿ���ڵ����Ĳ���
    double _p = 0.25;//��ÿ���ڵ��һ��ĸ���Ϊp

    int _randomLeve() {
        //�����ڵ���������
        int leve = 1;
        while (rand() < (RAND_MAX * _p) && leve < _maxleve) {
            //rand()<(RAND_MAX*_p�ĸ���Ϊ0.25
            leve += 1;
        }
        return leve;
    }

    //����ĳ���ڵ������ǰһ��ָ��
    std::vector<Node*>_findPrevNode(double num) {
        //����Ҫ���Ҳ����λ�ã�ǰһ���ڵ�ָ��
        Node* cur = _head;
        int curleve = cur->_nextV.size() - 1;
        std::vector<Node*>prev(curleve + 1, _head);//����ǰһ���ڵ��ָ��
        while (curleve >= 0) {//�����û���ҵ��ڵ�����һ��ʱ��Ҫ����ѭ��
            //�ڵ������ƶ�ʱ������ǰһ���ڵ�ָ������
            if (cur->_nextV[curleve] && cur->_nextV[curleve]->val < num) {
                //���������ߣ�������һ���ڵ�
                //�����������һ���ڵ�Ϊ��ʱ��Ҫ������ҲҪ�����ƶ�
                cur = cur->_nextV[curleve];
            }
            else if (cur->_nextV[curleve] == nullptr || cur->_nextV[curleve]->val >= num) {
                //����prev���飬����������
                prev[curleve] = cur;
                curleve -= 1;
            }
        }
        return prev;
    }

public:
    Skiplist() {
        //ͷ�ڵ�ֵΪ-1,��ʼΪ��һ��
        _head = new Node(-1, 1, -1);                                             //ͷ���idΪ-1��
        srand(time(0));//���������
    }

    Node* search(double target) {
        Node* cur = _head;
        int curleve = cur->_nextV.size() - 1;
        double dmin = 1e9;
        Node* tmp = nullptr;
        while (curleve >= 0) {//�����û���ҵ��ڵ�����һ��ʱ��Ҫ����ѭ��
            if (cur->_nextV[curleve] && cur->_nextV[curleve]->val < target) {
                //���������ߣ�������һ���ڵ�
                //�����������һ���ڵ�Ϊ��ʱ��Ҫ������ҲҪ�����ƶ�
                cur = cur->_nextV[curleve];
            }
            else if (cur->_nextV[curleve] == nullptr || cur->_nextV[curleve]->val > target) {
                //����������
                curleve -= 1;
            }
            else {
                //�ҵ�������ڵ�
                return cur->_nextV[curleve];
            }
        }
        return nullptr;
    }

    int searchmin(double target) {
        Node* cur = _head;
        int curleve = cur->_nextV.size() - 1;
        double dmin = 1e9;
        Node* tmp = nullptr;
        while (curleve >= 0) {//�����û���ҵ��ڵ�����һ��ʱ��Ҫ����ѭ��
            if (cur->_nextV[curleve] && cur->_nextV[curleve]->val < target) {
                //���������ߣ�������һ���ڵ�
                //�����������һ���ڵ�Ϊ��ʱ��Ҫ������ҲҪ�����ƶ�
                cur = cur->_nextV[curleve];
                if (dmin > abs(target - cur->val)) {
                    tmp = cur;
                    dmin = abs(target - cur->val);
                }
                    
            }
            else if (cur->_nextV[curleve] == nullptr || cur->_nextV[curleve]->val > target) {
                //����������
                curleve -= 1;
            }
            else {
                //�ҵ�������ڵ�
                if(dmin > abs(target - cur->_nextV[curleve] -> val))
                    return cur->_nextV[curleve]->id;
                return tmp->id;
            }
        }
        int tid = -1;
        if (tmp) {
            tid = tmp->val;
        }
        if (tmp && tmp->_nextV[curleve+1] && abs(tmp->val - target) > abs(tmp->_nextV[curleve+1]->val - target))
            tid = tmp->_nextV[curleve]->id;
        return tid;
    }

    bool searchp(double target) {
        Node* cur = _head;
        int curleve = cur->_nextV.size() - 1;
        while (curleve >= 0) {//�����û���ҵ��ڵ�����һ��ʱ��Ҫ����ѭ��
            if (cur->_nextV[curleve] && cur->_nextV[curleve]->val < target) {
                //���������ߣ�������һ���ڵ�
                //�����������һ���ڵ�Ϊ��ʱ��Ҫ������ҲҪ�����ƶ�
                cur = cur->_nextV[curleve];
            }
            else if (cur->_nextV[curleve] == nullptr || cur->_nextV[curleve]->val > target) {
                //����������
                curleve -= 1;
            }
            else {
                //�ҵ�������ڵ�
                return true;
            }
        }
        return false;
    }

    void add(double num,int _id) {
        //�������ʱ��������
        std::vector<Node*>prev = _findPrevNode(num);
        int newLeve = _randomLeve();
        Node* newNode = new Node(num, newLeve, _id);
        //���newLeve����ͷ�ڵ������������ѡ������head����
        if (newLeve > _head->_nextV.size()) {
            _head->_nextV.resize(newLeve);//���������ڵ㵼��ͷ�ڵ��������
            _head->_frontV.resize(newLeve);
            prev.resize(newLeve, _head);//����Ĳ���ָ��ͷ�ڵ㡣
        }

        //����ǰ��ڵ�
        for (size_t i = 0; i < newLeve; i++) {
            if (prev[i]->_nextV[i] == nullptr ) {
                /*  ��˫������ֻ��һ��ͷ���ʱ */
                newNode->_nextV[i] = prev[i]->_nextV[i];		//  ����ָ������
                prev[i]->_nextV[i] = newNode;
                newNode->_frontV[i] = prev[i];			//  ���������front ָ��ͷ���
            }
            else {
                prev[i]->_nextV[i]->_frontV[i] = newNode;
                newNode->_frontV[i] = prev[i];
                newNode->_nextV[i] = prev[i]->_nextV[i];
                prev[i]->_nextV[i] = newNode;
            }
        }
    }

    //���Դ�ӡÿһ������
    void _Print() {
        int leve = _head->_nextV.size();
        for (int i = leve - 1; i >= 0; i--) {
            Node* cur = _head;
            Node* p = _head;
            //��ӡÿһ�������
            while (cur != nullptr) {
                std::cout << cur->val << "->";
                cur = cur->_nextV[i];
            }
            std::cout << "nullptr" << std::endl<<std::endl;

            /*
            while ((cur->_frontV[i]) != nullptr) {
                std::cout << cur->val << "->";
                cur = cur->_frontV[i];
            }
            std::cout << "nullptr" << std::endl << std::endl;
            */
        }

    }

    bool erase(double num, int _id) {
        //ɾ���ڵ㣬�ҵ����ڵ��ǰһ��ָ�룬�ڰ�Ҫɾ���ڵ��������һ��ָ����������
        std::vector<Node*>prev = _findPrevNode(num);
        if (prev[0]->_nextV[0] == nullptr || prev[0]->_nextV[0]->val != num) {
            //�������²㶼�Ҳ�������ڵ㣬˵���Ҳ����ڵ�
            return false;
        }
        //�ҵ��˶�Ӧ�ڵ�
        Node* del = prev[0]->_nextV[0];
        for (size_t i = 0; i < del->_nextV.size(); i++) {
            //����del�ڵ��ÿһ��ǰ��ָ��
            prev[i]->_nextV[i] = del->_nextV[i];
        }
        delete del;

        //���ɾ�����Լ�������ĸ߶ȣ��������Ч��
        int _headLeve = _head->_nextV.size() - 1;
        while (_headLeve >= 0) {
            if (_head->_nextV[_headLeve] == nullptr) {
                _headLeve -= 1;
            }
            else {
                break;
            }
        }
        _head->_nextV.resize(_headLeve + 1);
        return true;
    }
};

/**
 * Your Skiplist object will be instantiated and called as such:
 * Skiplist* obj = new Skiplist();
 * bool param_1 = obj->search(target);
 * obj->add(num);
 * bool param_3 = obj->erase(num);
 */
