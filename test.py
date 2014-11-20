#!/usr/bin/python3

import random
import sys
import os

def run(num_of_vertices, source, num_of_proc, out = False):
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
	if out:
		return output
	else:
		return None

def main():
	if (len(sys.argv) < 2):
		print("Run with ./test.py <number of runs> [-o]")
		return
	if (len(sys.argv) == 3 and sys.argv[2] == "-o"):
		out = True
	else:
		out = False
	arg = int(sys.argv[1])
	for i in range(arg):
		print("-------------------------------------------------------")
		print("Run #" + str(i) + ":")
		num_of_vertices = random.randrange(10000) + 10
		source = random.randrange(num_of_vertices)
		num_of_proc = min(random.randrange(10) + 5, num_of_vertices) 
		print("Number of vertices: " + str(num_of_vertices))
		print("Source: " + str(source))
		print("Number of processses: " + str(num_of_proc))
		print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print("Result:")
		output = run(num_of_vertices, source, num_of_proc, out)
		print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		if (output):
			print("********************************************************")
			print("Program output:")
			print(output)
			print("********************************************************")

if __name__ == "__main__":
	main()
