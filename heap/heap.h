#include <utility>
#include <vector>

using namespace std;

typedef pair<int, int> edge;

class Heap { 
private:
	const int INF = 10000;
	vector<edge> elements_;
	vector<int> map_;
	inline int left_(int i);
	inline int right_(int i);
	inline int parent_(int i);
	void min_heapify_(int i);
	void build_heap_();
public:
	Heap(edge* elements, int size);
	Heap(vector<edge> elemetns);
	void decrease_key(int i, int value);
	edge find_min();
	void delete_min();
	void remove(int i);
	int get_value(int i);
	int size();
	bool contain(int i);
	vector<edge> get_elements();
	int get_top();
};


//Индекс и позицию хранить в отдельном массиве
