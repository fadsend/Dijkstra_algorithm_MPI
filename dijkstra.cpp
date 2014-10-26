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
	write_graph_to_file(adjacency_matrix, num_of_vertices);
	
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
	/*
		This function does dijkstra's algorithm for a random graph.
	*/
	
	//Generate adjacency matrix for graph
	
	//Initialize destination and visited arrays
	int* destination = new int[num_of_vertices];
	
	//Initialize arrays
	for (int i = 0; i < num_of_vertices; ++i) {
		destination[i] = INF;
	}
	destination[source] = 0;

	//Provides relaxation
	auto relax = [destination, adjacency_matrix, num_of_vertices] (auto i, auto j) {
		if (destination[j] > destination[i] + adjacency_matrix[i * num_of_vertices + j]) {
			destination[j] = destination[i] + adjacency_matrix[i * num_of_vertices + j];
		}
	};
	
	//Init and fill selected array
	bool* selected = new bool[num_of_vertices];
	for (int i = 0; i < num_of_vertices; ++i) {
		selected[i] = false;
	}

	//Search for the least not selected vertex
	auto all_selected = [selected, num_of_vertices] () {
		for (int i = 0; i < num_of_vertices; ++i) {
			if (!selected[i]) {
				return false;
			}
		}
		return true;
	};
	
	auto find_min = [selected, num_of_vertices, destination, source] () {
		int min_dest = INF;
		int idx = 0;
		for (int i = 0; i < num_of_vertices; ++i) {
			if (destination[i] < min_dest && !selected[i]) {
				min_dest = destination[i];
				idx = i;
			}
		}
		return idx;
	};

	double time1 = MPI_Wtime();
	int iteration_count = 0;

	//Repeat until all vertices will be selected
	while (!all_selected()) {
		
		//Find not selected vertex with minimum distance
		int current = find_min();
					
		//Relax all adjacent to the current vertices
		for (int j = 0; j < num_of_vertices; ++j) {
			if (adjacency_matrix[current * num_of_vertices + j] != INF)
				relax(current, j);
		}
		
		//Label current vertex as selected
		selected[current] = true;
				
		iteration_count += 1;

		if (iteration_count > num_of_vertices)
			break;
	}
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
	delete[] selected;
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
	//Test print
	/*for(int i = 0; i< size; ++i) {
		cout << result[i] << " ";
	}
	cout << endl;*/
	
	return result;
}

void parallel_dijkstras_algorithm(int num_of_vertices, int source) {
	int* adjacency_matrix = nullptr;
	int* destination = nullptr;
	bool* selected = new bool[num_of_vertices];
	int* splits = nullptr;
	int* local_dest = nullptr;
	int* lengths = nullptr;
	int* displacements = nullptr;
	int length = 0;
	int start = 0, end = 0;
	int* partition = nullptr;
	int* matrix_disp = nullptr;
	int* current_vertex = new int[num_of_vertices];
	int* matrix_lengths = nullptr;
	for (int i = 0; i < num_of_vertices; ++i) {
		selected[i] = false;
	}
	double time1, time2;
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (proc_rank == 0) {	
		adjacency_matrix = generate_graph(num_of_vertices);
		serial_dijkstras_algorithm(adjacency_matrix, num_of_vertices, source, false);
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
		}
		displacements[size - 1] = splits[size - 1];
		lengths[size - 1] = num_of_vertices - splits[size - 1];
		for (int i = 0; i < size; ++i) {
			matrix_disp[i] = displacements[i] * num_of_vertices;
			matrix_lengths[i] = lengths[i] * num_of_vertices;
		}
	} else {
	//	adjacency_matrix = new int[num_of_vertices * num_of_vertices];
	}
	//MPI_Bcast(adjacency_matrix, num_of_vertices * num_of_vertices, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(splits, 1, MPI_INT, &start, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(lengths, 1, MPI_INT, &length, 1, MPI_INT, 0, MPI_COMM_WORLD);
	local_dest = new int[length];
	partition = new int[length * num_of_vertices];

	MPI_Scatterv(adjacency_matrix, matrix_lengths, matrix_disp, MPI_INT, partition, length * num_of_vertices, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(destination, lengths, splits, MPI_INT, local_dest, length, MPI_INT, 0, MPI_COMM_WORLD);
	if (proc_rank == 0) {
		//delete[] adjacency_matrix;
		for (int i = 0; i < size - 1; ++i) {
			splits[i] = splits[i + 1];
		}
		splits[size - 1] = num_of_vertices;
	}
	MPI_Scatter(splits, 1, MPI_INT, &end, 1, MPI_INT, 0, MPI_COMM_WORLD);
	auto relax = [local_dest, current_vertex, num_of_vertices, start] (auto min_value, auto i , auto j) {
		if (local_dest[j] > min_value + current_vertex[ j + start]) {
			local_dest[j] = min_value + current_vertex[j + start];
		}
	};
	auto all_selected = [selected, num_of_vertices] () {
		for (int i = 0; i < num_of_vertices; ++i) {
			if (!selected[i]) {
				return false;
			}
		}
		return true;
	};
	
	auto find_min = [selected, num_of_vertices, destination, source, start, end, local_dest] () {
		int min_dest = INF;
		int idx = 0;
		for (int i = 0; i < end - start; ++i) {
			if (local_dest[i] < min_dest && !selected[i + start]) {
				min_dest = local_dest[i];
				idx = i;
			}
		}
		return make_pair(min_dest, idx);
	};

	int iteration_count = 0;

	time1 = MPI_Wtime();
	//Repeat until all vertices will be selected
	
	while (!all_selected()) {
		
		//Find not selected vertex with minimum distance
		auto min_pair  = find_min();
		int mins[2]; mins[0] = min_pair.first; mins[1] = min_pair.second + start;
		int all_mins[2];
		MPI_Allreduce(mins, all_mins, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
		/*int curr_proc = -1;
		if (all_mins[1] >= start && all_mins[1] < end) {
			curr_proc = proc_rank;
			//cout << start << endl;
			//int tmp_vertex = new int[num_of_vertices];
			for (int i = 0; i < num_of_vertices; ++i) {
				current_vertex[i] = partition[(all_mins[1] - start)* num_of_vertices + i];
			}
			if (proc_rank != 0) {
				if (curr_proc != -1)
				MPI_Send(current_vertex, num_of_vertices, MPI_INT, 0, 0, MPI_COMM_WORLD);
			}
		//		cout << "( " << partition[(all_mins[1] - start)* num_of_vertices + i] << ", " << adjacency_matrix[all_mins[1] * num_of_vertices + i] << ")" << endl;
		} 
		if (proc_rank == 0) {
			MPI_Status status;
			if (curr_proc != -1)
			MPI_Recv(current_vertex, num_of_vertices, MPI_INT, curr_proc, 0, MPI_COMM_WORLD, &status);
		//	curr_proc = 0;*/
		if (proc_rank == 0) {
			for(int i = 0; i < num_of_vertices; ++i) {
				current_vertex[i] = adjacency_matrix[all_mins[1] * num_of_vertices + i];
			}
		}
		//cout << curr_proc << endl;
		MPI_Bcast(current_vertex, num_of_vertices, MPI_INT, 0, MPI_COMM_WORLD);
		//Relax all adjacent to the current vertices
		for (int j = 0; j < end - start; ++j) {
			if (current_vertex[j + start] != INF)
				relax(all_mins[0], all_mins[1], j);
		}
		//Label current vertex as selected
		selected[all_mins[1]] = true;
				
		iteration_count += 1;

		if (iteration_count > num_of_vertices)
			break;
	}

	MPI_Gatherv(local_dest, length, MPI_INT, destination, lengths, displacements, MPI_INT, 0, MPI_COMM_WORLD);
	//Printing results
	if (proc_rank == 0) {
		time2 = MPI_Wtime();
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
	delete[] destination;
	if (proc_rank != 0) {
		//delete[] adjacency_matrix;
		delete[] splits;
	}
	delete[] selected;
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
