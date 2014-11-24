#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <utility>
#include <random>
#include <limits>
#include <fstream>
#include <vector>
#include <queue>
#include <set>
#include "heap/heap.h"
#include "mpi.h"

using namespace std;

static const int MAX_WEIGHT = 100;
static const double EDGES_PROBABILITY = 1;
static const int INF = MAX_WEIGHT * 100;
static const int START = 0;
static const int END = 1;

int proc_rank, size;

template<typename T>
void p(T a) { cout << a << endl; }

void write_graph_to_file(const int* adjacency_matrix, int num_of_vertices) {
	/*
		This function write the structure of the graph to the file 
		in dot language.
	*/
	ofstream out_file;
	out_file.open("out.dot");
	out_file << "graph my_graph {\n";
	bool* writed = new bool[num_of_vertices];
	for (int i = 0; i < num_of_vertices; ++i)
		writed[i] = false;
	for (int i = 0; i < num_of_vertices; ++i) {
		for (int j = i; j < num_of_vertices; ++j) {
			if (adjacency_matrix[i * num_of_vertices + j] != INF && adjacency_matrix[i * num_of_vertices + j] != 0) {
				out_file << i << " -- " << j  << "[label=" << adjacency_matrix[i * num_of_vertices + j] << "];\n";
				writed[i] = writed[j] = true;
			}
		}
	}
	for (int i = 0; i < num_of_vertices; ++i) {
		if (!writed[i]) {
			out_file << i << " ;\n";
		}
	}
	out_file << "}\n";
	out_file.close();
	delete writed;
}

int* generate_graph(int num_of_vertices) {
	/*
		This function creates a random graph which represents 
		by an adjacency matrix.
	*/
	
	//Allocating memory for the matrix
	int* adjacency_matrix = new int[num_of_vertices * num_of_vertices]; 

	//Returns true with specified probability 
	auto has_edge = [] () {
		return (rand() % (int) (1 / EDGES_PROBABILITY)) == 0;
	};

	//Returns random weigth
	auto get_weight = []() {
		return (rand() % MAX_WEIGHT) + 1; 
	};
	
	srand(time(0));
	
	//Fills the matrix with random values
	for (int i = 0; i < num_of_vertices; ++i) {
		for (int j = i; j < num_of_vertices; ++j) {
			if (i == j) {
				adjacency_matrix[i * num_of_vertices + j] = 0;
				continue;
			} 
			if (has_edge()) {
				adjacency_matrix[i * num_of_vertices + j] = adjacency_matrix[j * num_of_vertices + i] = get_weight();
			} else {
				adjacency_matrix[i * num_of_vertices + j] = adjacency_matrix[j * num_of_vertices + i] = INF;
			}
		}
	}

	//Writing structure of the graph to a file
//	write_graph_to_file(adjacency_matrix, num_of_vertices);
	
	//Test print
	/*cout << "Adjacency matrix for the graph:" << endl;
	for (int i = 0; i < num_of_vertices; ++i) {
		for (int j = 0; j < num_of_vertices; ++j) {
			if (adjacency_matrix[i * num_of_vertices + j] == INF) {
				printf("  inf ");
			} else {
				printf("%5d ", adjacency_matrix[i * num_of_vertices + j]);
			}
		}
		cout << endl;
	}*/

	return adjacency_matrix;
}

void serial_dijkstras_algorithm(int* adjacency_matrix, int num_of_vertices, int source, bool del = true) {
	//	This function does dijkstra's algorithm for a random graph.
	
	double time1 = MPI_Wtime();
	//Initialize destination and visited arrays
	int* destination = new int[num_of_vertices];
	
	//Initialize arrays
	for (int i = 0; i < num_of_vertices; ++i) {
		destination[i] = INF;
	}
	destination[source] = 0;

	//Forming a priority queue
	vector<edge> vec(num_of_vertices);
	for (int i = 0; i < num_of_vertices; ++i) {
		vec[i] = make_pair(i, destination[i]);
	}
	Heap priority_queue(vec);
	vec.clear();
	
	//Start time
	
	int iteration_count = 0;
	
	while(iteration_count != num_of_vertices) {
	
		//Find not selected vertex with minimum distance
		auto min_vertex = priority_queue.find_min();
		priority_queue.delete_min();	
		
		int current = min_vertex.first;
		int current_value = min_vertex.second;

		//Fix distances from adjacent to the current vertices
		for (int j = 0; j < num_of_vertices; ++j) {
			if (adjacency_matrix[current * num_of_vertices + j] != INF) {
				if (priority_queue.get_value(j) > current_value + adjacency_matrix[current * num_of_vertices + j]) {
					priority_queue.decrease_key(j, priority_queue.get_value(j) - (current_value + adjacency_matrix[current * num_of_vertices + j])); 	
					destination[j] = current_value + adjacency_matrix[current * num_of_vertices + j];
				}
			}
		}
		
		iteration_count += 1;
	}

	//End time
	double time2 = MPI_Wtime();
	
	//Printing results
	cout << "Serial\n"; //from " << source << ": " <<  endl;
	for (int i = 0; i < num_of_vertices; ++i) {
		if (destination[i] == INF) {
			cout << "inf ";
		} else {
			cout << destination[i] << " ";
		}
	}
	cout << endl;
	cout << "Time: " << time2 - time1 << endl;;
	
	//Free allocated memory
	delete[] destination;
	if (del)
		delete[] adjacency_matrix;
}

int* split_vertices(int num_of_vertices) {
	int* result = new int[size];
	bool is_divide = (num_of_vertices % size) == 0;
	int step = (is_divide) ? num_of_vertices / size : num_of_vertices / size + 1;
	int index = 0;
	for (int i = 0; i < num_of_vertices; i += step) {
		if (!is_divide && ((num_of_vertices - i) == (step - 1) * (size - index))) {
			step -= 1;
			is_divide = true;
		}
		result[index] = i;
		index += 1;
	}
	return result;
}

void parallel_dijkstras_algorithm(int num_of_vertices, int source) {
	int* adjacency_matrix = nullptr;//d
	int* destination = nullptr;//d
	int* splits = nullptr;//d
	int* lengths = nullptr;//d
	int* displacements = nullptr;//d
	int length = 0;
	int start = 0, end = 0;
	int* matrix_disp = nullptr;//d
	int* current_vertex = new int[num_of_vertices];//d
	int* matrix_lengths = nullptr;//d
	double time1, time2;
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int* id = new int[size + 1];//d
	if (proc_rank == 0) {	
		adjacency_matrix = generate_graph(num_of_vertices);
		serial_dijkstras_algorithm(adjacency_matrix, num_of_vertices, source, false);
		cout.flush();
		splits = split_vertices(num_of_vertices);
		lengths = new int[size];
		matrix_lengths = new int[size];
		displacements = new int[size];
		matrix_disp = new int[size];
		destination = new int[num_of_vertices];
		for (int i = 0; i < num_of_vertices; ++i) {
			destination[i] = INF;
		}
		destination[source] = 0;
		for (int i = 0; i < size - 1; ++i) {
			lengths[i] = splits[i + 1] - splits[i];
			displacements[i] = splits[i];
			id[i] = splits[i];
		}
		id[size - 1] = splits[size - 1];
		id[size] = num_of_vertices;
		displacements[size - 1] = splits[size - 1];
		lengths[size - 1] = num_of_vertices - splits[size - 1];
		for (int i = 0; i < size; ++i) {
			matrix_disp[i] = displacements[i] * num_of_vertices;
			matrix_lengths[i] = lengths[i] * num_of_vertices;
		}
	}
	MPI_Bcast(id, size + 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(splits, 1, MPI_INT, &start, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(lengths, 1, MPI_INT, &length, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	int* local_dest = new int[length];
	int* partition = new int[length * num_of_vertices];

	MPI_Scatterv(adjacency_matrix, matrix_lengths, matrix_disp, MPI_INT, partition, length * num_of_vertices, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(destination, lengths, splits, MPI_INT, local_dest, length, MPI_INT, 0, MPI_COMM_WORLD);
	if (proc_rank == 0) {
		delete[] matrix_disp;
		delete[] matrix_lengths;
		delete[] adjacency_matrix;
		for (int i = 0; i < size - 1; ++i) {
			splits[i] = splits[i + 1];
		}
		splits[size - 1] = num_of_vertices;
	}
	MPI_Scatter(splits, 1, MPI_INT, &end, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	//Construct a priority_queue
	vector<edge> vec;
	for (int i = 0; i < end - start; ++i) {
		vec.push_back(make_pair(i, local_dest[i]));
	}
	Heap priority_queue(vec);
	vec.clear();

	if (proc_rank == 0) {
		delete[] splits;
	}
	int iteration_count = 0;
		time1 = MPI_Wtime();
	while(iteration_count != num_of_vertices) {
		
		//Find not selected vertex with minimum distance
		pair<int, int> min_pair;
		int local_min[2];
		int global_min[2];
		if (priority_queue.size() != 0) {
			min_pair  = priority_queue.find_min();
			local_min[0] = min_pair.second; 
			local_min[1] = min_pair.first + start;
		} else {
			local_min[0] = INF;
			local_min[1] = -1;
		}

		//Reduce all local local_min to global one
		MPI_Allreduce(local_min, global_min, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
		
		//Find a process what has the next vertex
		int curr_proc;
		for (int s = 0; s < size + 1; ++s) {
			if (global_min[1] >= id[s] && global_min[1] < id[s + 1]) {
				curr_proc = s;
				break;
			}
		}
		
		//Copy the current raw from the process
		if (proc_rank == curr_proc) {
			priority_queue.delete_min();
			for (int i = 0; i < num_of_vertices; ++i) {
				current_vertex[i] = partition[(global_min[1] - start)* num_of_vertices + i];
			}
		}

		//Broadcast the current vertex raw
		MPI_Bcast(current_vertex, num_of_vertices, MPI_INT, curr_proc, MPI_COMM_WORLD);
		
		//Relax all adjacent to the current vertices
		for (int j = 0; j < end - start; ++j) {
			if (current_vertex[j + start] != INF) {
				if (priority_queue.get_value(j) > global_min[0] + current_vertex[ j + start]) {
					int tmp = local_dest[j];
					local_dest[j] = global_min[0] + current_vertex[ j + start];
					if (priority_queue.contain(j)) {
						priority_queue.decrease_key(j, tmp - (global_min[0] + current_vertex[j + start])); 	
					} 
				} 
			}
		}
				
		iteration_count += 1;
	}
	if (proc_rank == 0) {
		time2 = MPI_Wtime();
	}

	//Gather all local destination to the final destination array
	MPI_Gatherv(local_dest, length, MPI_INT, destination, lengths, displacements, MPI_INT, 0, MPI_COMM_WORLD);
	
	//Printing results
	if (proc_rank == 0) {
		//time2 = MPI_Wtime();
		cout << "Parallel:\n";// from " << source << ": " <<  endl;
		for (int i = 0; i < num_of_vertices; ++i) {
			if (destination[i] == INF) {
				cout << "inf ";
			} else {
				cout << destination[i] << " ";
			}
		}
		cout << endl;
		cout << "Time: " << time2 - time1 << endl;;
	}

	//Free allocated memory
	delete[] id;
	delete[] local_dest;
	if (proc_rank == 0) {
		delete[] destination;
		delete[] displacements;
		delete[] lengths;
		//delete[] adjacency_matrix;
	}
	delete[] current_vertex;
	delete[] partition;
	//delete[] selected;
}

int main(int argc, char* argv[])
{
	if (argc < 3) {
		cout << "Error! Too few arguments\n";
		return 1;
	}
	int num_of_vertices = atoi(argv[1]);
	int source = atoi(argv[2]);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (size == 1) {
		int* adjacency_matrix = generate_graph(num_of_vertices);
		serial_dijkstras_algorithm(adjacency_matrix, num_of_vertices, source);
	} else {
		parallel_dijkstras_algorithm(num_of_vertices, source);
	}
	MPI_Finalize();
	return 0;
}
