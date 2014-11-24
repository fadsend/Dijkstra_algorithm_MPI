#include "heap.h"
#include <iostream>

int Heap::left_(int i) {
	return 2 * i;
}

int Heap::right_(int i) {
	return 2 * i + 1;
}

int Heap::parent_(int i) {
	return i / 2;
}

void Heap::build_heap_() {
	for (int i = elements_.size() / 2; i >= 0; --i) {
		min_heapify_(i);
	}
}

void Heap::min_heapify_(int i) {
	int l = left_(i);
	int r = right_(i);
	int smallest;
	if (l < elements_.size() && elements_[l].second < elements_[i].second) {
		smallest = l;
	} else {
		smallest = i;
	}
	if (r < elements_.size() && elements_[r].second < elements_[smallest].second) {
		smallest = r;
	}
	if (smallest != i) {
		swap(elements_[smallest], elements_[i]);
		swap(map_[elements_[smallest].first], map_[elements_[i].first]);
		min_heapify_(smallest);
	}
}

Heap::Heap(vector<edge> elements) : elements_(elements), map_(elements.size()) {
	for (size_t i = 0; i < map_.size(); ++i) {
		map_[elements_[i].first] = i;
	}
	build_heap_();
}

Heap::Heap(edge* elements, int size) : map_(size), elements_(elements, elements + size) {
	for (size_t i = 0; i < size; ++i) {
		map_[elements_[i].first] = i;
	}
	build_heap_();
}
void Heap::delete_min() {
//	cout << "min: " << elements_[0].first << " " <<elements_[0].second << endl; 
	map_[elements_[0].first] = -1;
	map_[elements_[elements_.size() - 1].first] = 0;
//	cout << "deleting : " << elements_[0].first << ", " << elements_[0].second << endl;
	swap(elements_[0], elements_[elements_.size() - 1]);
	//swap(map_[elements_[0].first], map_[elements_[elements_.size() - 1].first]);
	elements_.pop_back();
	min_heapify_(0);
//	cout << "start: ";
//	for(auto i : elements_) 
//		cout << "( " << i.first << ", " << i.second << ")";
//	cout << endl;
}

void Heap::decrease_key(int i, int value) {
	int location = map_[i];
	elements_[location].second -= value;
	if (elements_[location].second < 0) {
//cout << i << " v = " <<  value <<  "el = " << elements_[location].first << "old = " << elements_[location].second << endl;
//		cout << "loc = " << location << endl;
//		cout << elements_[location].first << " " << elements_[location].second << endl;
		throw "Negative edge weight";
	}
	while (location > 0 and elements_[location].second < elements_[parent_(location)].second) {
		swap(elements_[location], elements_[parent_(location)]);
		swap(map_[elements_[location].first], map_[elements_[parent_(location)].first]);
		location = parent_(location);
	}
//	cout << "loc = " << location << endl;
}

void Heap::remove(int i) {
	int location = map_[i];
	map_[elements_[location].first] = -1;
	map_[elements_[elements_.size() - 1].first] = location;
	swap(elements_[location], elements_[elements_.size() - 1]);
	elements_.pop_back();
	min_heapify_(0);
}

int Heap::get_value(int i) {
	return elements_[map_[i]].second;
}

int Heap::size() {
	return elements_.size();
}



vector<edge> Heap::get_elements() {
	return elements_;
}

int Heap::get_top() {
	return elements_[0].second;
}
