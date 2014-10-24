#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <utility>
#include <random>
#include <limits>
#include <fstream>
#include "mpi.h"

using namespace std;

static const int MAX_WEIGHT = 100;
static const double EDGES_PROBABILITY = 0.5;
static const int INF = MAX_WEIGHT * 100;

int proc_rank, size;

void write_graph_to_file(const int* adjacency_matrix, size_t num_of_vertices) {
	/*
		This function write the structure of the graph to the file 
		in dot language.
	*/
	ofstream out_file;
	out_file.open("out.dot");
	out_file << "graph my_graph {\n";
	bool* writed = new bool[num_of_vertices];
	for (size_t i = 0; i < num_of_vertices; ++i)
		writed[i] = false;
	for (size_t i = 0; i < num_of_vertices; ++i) {
		for (size_t j = i; j < num_of_vertices; ++j) {
			if (adjacency_matrix[i * num_of_vertices + j] != INF && adjacency_matrix[i * num_of_vertices + j] != 0) {
				out_file << i << " -- " << j  << "[label=" << adjacency_matrix[i * num_of_vertices + j] << "];\n";
				writed[i] = writed[j] = true;
			}
		}
	}
	for (size_t i = 0; i < num_of_vertices; ++i) {
		if (!writed[i]) {
			out_file << i << " ;\n";
		}
	}
	out_file << "}\n";
	out_file.close();
	delete writed;
}

int* generate_graph(size_t num_of_vertices) {
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
	for (size_t i = 0; i < num_of_vertices; ++i) {
		for (size_t j = i; j < num_of_vertices; ++j) {
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
	cout << "Adjacency matrix for the graph:" << endl;
	for (size_t i = 0; i < num_of_vertices; ++i) {
		for (size_t j = 0; j < num_of_vertices; ++j) {
			if (adjacency_matrix[i * num_of_vertices + j] == INF) {
				printf("  inf ");
			} else {
				printf("%5d ", adjacency_matrix[i * num_of_vertices + j]);
			}
		}
		cout << endl;
	}

	return adjacency_matrix;
}

void serial_dijkstras_algorithm(size_t num_of_vertices, size_t source) {
	/*
		This function does dijkstra's algorithm for a random graph.
	*/
	
	//Generate adjacency matrix for graph
	const int* adjacency_matrix = generate_graph(num_of_vertices);
	
	//Initialize destination and visited arrays
	int* destination = new int[num_of_vertices];
	
	//Initialize arrays
	for (size_t i = 0; i < num_of_vertices; ++i) {
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
	for (size_t i = 0; i < num_of_vertices; ++i) {
		selected[i] = false;
	}

	//Search for the least not selected vertex
	auto all_selected = [selected, num_of_vertices] () {
		for (size_t i = 0; i < num_of_vertices; ++i) {
			if (!selected[i]) {
				return false;
			}
		}
		return true;
	};
	
	auto find_min = [selected, num_of_vertices, destination, source] () {
		int min_dest = destination[0];
		size_t idx = source;
		for (size_t i = 0; i < num_of_vertices; ++i) {
			if (destination[i] < min_dest && !selected[i]) {
				min_dest = destination[i];
				idx = i;
			}
		}
		return idx;
	};

	size_t iteration_count = 0;

	//Repeat until all vertices will be selected
	while (!all_selected()) {
		
		//Find not selected vertex with minimum distance
		size_t current = find_min();
					
		//Relax all adjacent to the current vertices
		for (size_t j = 0; j < num_of_vertices; ++j) {
			if (adjacency_matrix[current * num_of_vertices + j] != INF)
				relax(current, j);
		}

		//Label current vertex as selected
		selected[current] = true;
				
		iteration_count += 1;

		if (iteration_count > num_of_vertices)
			break;
	}

	//Printing results
	cout << "Destinations from " << source << ": " <<  endl;
	for (size_t i = 0; i < num_of_vertices; ++i) {
		if (destination[i] == INF) {
			cout << "inf ";
		} else {
			cout << destination[i] << " ";
		}
	}
	cout << endl;

	//Free allocated memory
	delete[] destination;
	delete[] adjacency_matrix;
	delete[] selected;
}

void parallel_dijkstras_algorithm() {
	/*
		This function does parallel dijkstra's algorithm for a random graph.
	*/
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
	if (proc_rank == 0) {
		//size_t num_of_vertices, source;
		cout << "Please choose the number of vertices and source vertex:\n";
		//cin >> num_of_vertices >> source;
		//int* adjacency_matrix = generate_graph(num_of_vertices);
	}
}
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (size == 1) {
		size_t num_of_vertices, source;
		cout << "Please choose the number of vertices and source vertex:\n";
		cin >> num_of_vertices >> source;
		double time1 = MPI_Wtime();
		serial_dijkstras_algorithm(num_of_vertices, source);
		double time2 = MPI_Wtime();
		cout << "Time: " << time2 - time1 << endl;;
	} else {
		parallel_dijkstras_algorithm();
	}
	MPI_Finalize();
	return 0;
}
