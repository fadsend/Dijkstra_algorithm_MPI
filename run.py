#!/usr/bin/python3

import sys
import os

def run(num_of_vertices, source, num_of_proc):
	output = os.popen("mpirun -n " + str(num_of_proc) + " ./dijkstra " + str(num_of_vertices) + " " + str(source)).read()
	out_list = output.split("\n")
	parallel_out, serial_out = "", ""
	parallel_time, serial_time = "", ""
	for idx in range(len(out_list)):
		if out_list[idx].find("Parallel") == 0:
			parallel_out = out_list[idx + 1]
			parallel_time = out_list[idx + 2];
		if out_list[idx].find("Serial") == 0:
			serial_out = out_list[idx + 1]
			serial_time = out_list[idx + 2]
	parallel_out = parallel_out.split(" ")
	serial_out = serial_out.split(" ")
	is_error = False
	if parallel_out == [] or serial_out == []:
		print("Error!")
		is_error = True
	for i in range(min(len(parallel_out), len(serial_out))):
		if parallel_out != serial_out:
			print("Error!")
			is_error = True
			break
	if (not is_error):
		print("OK")
		print("Parallel time : " + parallel_time[6:])
		print("Serial time : " + serial_time[6:])
	else:
		print("Parallel : " + str(parallel_out))
		print("Serial   : " + str(serial_out))

if __name__ == "__main__":
	if (len(sys.argv) < 4):
		print("Format: ./run.py <num_of_vertices> <source> <num_of_proc>")
	else:
		run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
