#include <iostream>
#include "heap.h"


int main() {
	vector<edge> p;
	p.push_back(make_pair(0, 123123));
	p.push_back(make_pair(3, 12));
	p.push_back(make_pair(2, 3));
	p.push_back(make_pair(5, 437));
	p.push_back(make_pair(4, 23));
	p.push_back(make_pair(1, 0));
	Heap h(p);
	cout << h.find_min().first << " " << h.find_min().second<< endl;
	h.delete_min();
	h.decrease_key(0, 123122);
	cout << h.contain(1) << endl;
	cout << h.find_min().first << " " << h.find_min().second<< endl;
	h.delete_min();
	cout << h.contain(0) << endl;
	cout << h.find_min().first << " " << h.find_min().second<< endl;
	h.delete_min();
	cout << h.contain(2) << endl;
	cout << h.find_min().first << " " << h.find_min().second<< endl;
	h.delete_min();
	cout << h.contain(3) << endl;
	cout << h.find_min().first << " " << h.find_min().second<< endl;
	h.delete_min();
	cout << h.contain(4) << endl;
	cout << h.find_min().first << " " << h.find_min().second<< endl;
	h.delete_min();
	return 0;
}
