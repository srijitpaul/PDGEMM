from scipy import stats
import numpy as np
import pandas as pd
import pdb
from matplotlib import pyplot as plt



def gen_performance_plot(data_path, plot_path):

	""" Generate performance plots for a given performance input file (e.g. ../data/performance_seq_s2000_m5.tsv)

	Args:
		data_path = "../data/performance_seq_s2000_m5.tsv"
		plot_path = "../data/performance_vs_matsize.png"
	Returns:
		...
	Usage:
		gen_performance_plot("../data/performance_seq_s2000_m5.tsv", "../data/performance_vs_matsize.png")

	"""
	# 
	# Read in data
	#
	size_to_sample_to_times = {}

	count = 0

	fi = open(data_path,"r")

	headers = fi.readline().split("\t") #	Sample No.	Initialization_Time	Matrix_Fill_Time	GEMM_Time

	while True:
		count += 1 
		line = fi.readline() # e.g. 50	1	1.6047e-05s	3.2989e-05s	0.000558953s

		if line=="":
			break

		line_split = line.split("\t")

		if not 'matrix_size' in locals():
			matrix_size_start = int(line_split[0])

		#pdb.set_trace()

		matrix_size 		= int(line_split[0])
		sample_no 			= int(line_split[1])
		t_initialization 	= float(line_split[2].replace("s",""))
		t_matrixfill 		= float(line_split[3].replace("s",""))
		t_gemm 				= float(line_split[4].rstrip().replace("s",""))

		if matrix_size==(matrix_size_start*2):
			matrix_size_step = matrix_size - matrix_size_start

		try:
			size_to_sample_to_times[matrix_size][sample_no] = {
					"Initialization_Time": t_initialization,
					"Matrix_Fill_Time":    t_matrixfill,
					"GEMM_Time":           t_gemm
				}
		except KeyError:
			size_to_sample_to_times[matrix_size] = {
				sample_no:{
					"Initialization_Time": t_initialization,
					"Matrix_Fill_Time":    t_matrixfill,
					"GEMM_Time":           t_gemm
				}
			}
	fi.close()


	matrix_size_end = matrix_size
	sample_size 	= sample_no

	#
	# Statistics per matrix size
	# 
	size_to_t_initialization 	= np.zeros((sample_no,(matrix_size_end/matrix_size_step)))
	size_to_t_matrixfill 		= np.zeros((sample_no,(matrix_size_end/matrix_size_step)))
	size_to_t_gemm 				= np.zeros((sample_no,(matrix_size_end/matrix_size_step)))

	for matrix_size in size_to_sample_to_times.keys():
		for sample in size_to_sample_to_times[matrix_size].keys():
			size_to_t_initialization[sample-1, (matrix_size/matrix_size_step)-1] 	= size_to_sample_to_times[matrix_size][sample]['Matrix_Fill_Time']
			size_to_t_matrixfill[sample-1, (matrix_size/matrix_size_step)-1] 		= size_to_sample_to_times[matrix_size][sample]['Initialization_Time']
			size_to_t_gemm[sample-1, (matrix_size/matrix_size_step)-1] 				= size_to_sample_to_times[matrix_size][sample]['GEMM_Time']

	# save to file
	#np.savetxt("../data/tmp/scrap.tsv", size_to_t_gemm, delimiter="\t")

	t_gemm_means = []
	t_initialization_means = []
	t_matrixfill_means = []

	for i in range(0,len(size_to_t_gemm[1,:])):

		t_gemm_means.append(np.mean(size_to_t_gemm[:,i]))
		t_initialization_means.append(np.mean(size_to_t_initialization[:,i]))
		t_matrixfill_means.append(np.mean(size_to_t_matrixfill[:,i]))

	matrix_sizes = range(matrix_size_start,matrix_size_end+matrix_size_step,matrix_size_step)

	t_initialization_errors =  stats.sem(size_to_t_initialization,axis=0)
	t_matrixfill_errors 	=  stats.sem(size_to_t_matrixfill,axis=0)
	t_gemm_errors 			= stats.sem(size_to_t_gemm,axis=0)

	#
	# Plotting
	#
	# pd.read_csv(fpath,sep="\t") // DataFrame.from_csv

	f, (ax1, ax2, ax3) = plt.subplots(3, sharex=False, sharey=False)

	plt.subplots_adjust(hspace=.4)

	ax1.set_title("Initialization", loc="left")
	ax1.set_ylabel('Time (s)')
	
	ax1.errorbar(matrix_sizes, t_initialization_means, yerr=t_initialization_errors)

	ax2.set_title("Matrix Filling", loc="left")
	ax2.set_ylabel('Time (s)')
	
	ax2.errorbar(matrix_sizes, t_matrixfill_means, yerr=t_matrixfill_errors)

	ax3.set_title("DGEMM", loc="left")
	ax3.set_ylabel('Time (s)')
	ax3.set_xlabel('Matrix Size', fontsize=12)

	
	ax3.errorbar(matrix_sizes, t_gemm_means, yerr=t_gemm_errors)
	plt.savefig(plot_path)

	plt.close()

def gen_performance_stats(data_path):

	""" Generate performance plots for a given performance input file (e.g. ../data/performance_seq_s2000_m5.tsv)

	Args:
		data_path = "../data/performance_seq_s2000_m5.tsv"
	Returns:
		...
	Usage:
		matrix_sizes, t_initialization_means, t_matrixfill_means, t_gemm_means, t_initialization_errors, t_matrixfill_errors, t_gemm_errors = gen_performance_stats("../data/performance_seq_s2000_m5.tsv")

	"""
	# 
	# Read in data
	#
	size_to_sample_to_times = {}

	count = 0

	fi = open(data_path,"r")

	headers = fi.readline().split("\t") #	Sample No.	Initialization_Time	Matrix_Fill_Time	GEMM_Time

	while True:
		count += 1 
		line = fi.readline() # e.g. 50	1	1.6047e-05s	3.2989e-05s	0.000558953s

		if line=="":
			break

		line_split = line.split("\t")

		if not 'matrix_size' in locals():
			matrix_size_start = int(line_split[0])

		#pdb.set_trace()

		matrix_size 		= int(line_split[0])
		sample_no 			= int(line_split[1])
		t_initialization 	= float(line_split[2].replace("s",""))
		t_matrixfill 		= float(line_split[3].replace("s",""))
		t_gemm 				= float(line_split[4].rstrip().replace("s",""))

		if matrix_size==(matrix_size_start*2):
			matrix_size_step = matrix_size - matrix_size_start

		try:
			size_to_sample_to_times[matrix_size][sample_no] = {
					"Initialization_Time": t_initialization,
					"Matrix_Fill_Time":    t_matrixfill,
					"GEMM_Time":           t_gemm
				}
		except KeyError:
			size_to_sample_to_times[matrix_size] = {
				sample_no:{
					"Initialization_Time": t_initialization,
					"Matrix_Fill_Time":    t_matrixfill,
					"GEMM_Time":           t_gemm
				}
			}
	fi.close()


	matrix_size_end = matrix_size
	sample_size 	= sample_no

	#
	# Statistics per matrix size
	# 
	size_to_t_initialization 	= np.zeros((sample_no,(matrix_size_end/matrix_size_step)))
	size_to_t_matrixfill 		= np.zeros((sample_no,(matrix_size_end/matrix_size_step)))
	size_to_t_gemm 				= np.zeros((sample_no,(matrix_size_end/matrix_size_step)))

	for matrix_size in size_to_sample_to_times.keys():
		for sample in size_to_sample_to_times[matrix_size].keys():
			size_to_t_initialization[sample-1, (matrix_size/matrix_size_step)-1] 	= size_to_sample_to_times[matrix_size][sample]['Matrix_Fill_Time']
			size_to_t_matrixfill[sample-1, (matrix_size/matrix_size_step)-1] 		= size_to_sample_to_times[matrix_size][sample]['Initialization_Time']
			size_to_t_gemm[sample-1, (matrix_size/matrix_size_step)-1] 				= size_to_sample_to_times[matrix_size][sample]['GEMM_Time']

	# save to file
	#np.savetxt("../data/tmp/scrap.tsv", size_to_t_gemm, delimiter="\t")

	t_gemm_means 			= []
	t_initialization_means 	= []
	t_matrixfill_means 		= []

	for i in range(0,len(size_to_t_gemm[1,:])):

		t_gemm_means.append(np.mean(size_to_t_gemm[:,i]))
		t_initialization_means.append(np.mean(size_to_t_initialization[:,i]))
		t_matrixfill_means.append(np.mean(size_to_t_matrixfill[:,i]))

	matrix_sizes = range(matrix_size_start,matrix_size_end+matrix_size_step,matrix_size_step)

	t_initialization_errors =  stats.sem(size_to_t_initialization,axis=0)
	t_matrixfill_errors 	=  stats.sem(size_to_t_matrixfill,axis=0)
	t_gemm_errors 			=  stats.sem(size_to_t_gemm,axis=0)

	return matrix_sizes, t_initialization_means, t_matrixfill_means, t_gemm_means, t_initialization_errors, t_matrixfill_errors, t_gemm_errors

def plot_all(performance_versions, performance_data, plot_path):

	""" Plots all DGEMM timing segments: initialization, matrix filling, matrix multiplication, allowing comparison between different versions of the program

	Args:
		performance_versions = [performance_v1_name, performance_v2_name]
			performance_v1_name = "Sequential"
			performance_v2_name = "Parallel OpenMP"

		performance_data = [performance_v1, performance_v2, ... ]
			performance_v1 = [seq_mat_sizes, seq_t_init_means, seq_t_mfill_means, seq_t_gemm_means, seq_t_init_errors, seq_t_mfill_errors, seq_t_gemm_errors]
			performance_v2 = [omp_mat_sizes, omp_t_init_means, omp_t_mfill_means, omp_t_gemm_means, omp_t_init_errors, omp_t_mfill_errors, omp_t_gemm_errors]

		plot_path = "../data/performance.png"
	Returns:
		...
	Usage:
		plot_all([[seq_mat_sizes, seq_t_init_means, seq_t_mfill_means, seq_t_gemm_means, seq_t_init_errors, seq_t_mfill_errors, seq_t_gemm_errors],[omp_mat_sizes, omp_t_init_means, omp_t_mfill_means, omp_t_gemm_means, omp_t_init_errors, omp_t_mfill_errors, omp_t_gemm_errors]],"../data/performance.png")

	"""

	#
	# Plotting
	#
	# pd.read_csv(fpath,sep="\t") // DataFrame.from_csv

	#plot_path = "../data/performance.png"

	colors = ["b","g","r","c","m","y","k","w"]

	f, (ax1, ax2, ax3) = plt.subplots(3, sharex=False, sharey=False)

	plt.subplots_adjust(hspace=.4)

	for i,version in enumerate(performance_data):

		seq_mat_sizes 		= version[0]

		seq_t_init_means	= version[1]
		seq_t_mfill_means 	= version[2]
		seq_t_gemm_means	= version[3]

		seq_t_init_errors	= version[4]
		seq_t_mfill_errors	= version[5]
		seq_t_gemm_errors	= version[6]

		ax1.set_title("Initialization", loc="left")
		ax1.set_ylabel('Time (s)')
		ax1.errorbar(seq_mat_sizes, seq_t_init_means, yerr=seq_t_init_errors)
		#ax1.errorbar(omp_mat_sizes, omp_t_init_means, yerr=omp_t_init_errors)

		ax2.set_title("Matrix Filling", loc="left")
		ax2.set_ylabel('Time (s)')
		ax2.errorbar(seq_mat_sizes, seq_t_mfill_means, yerr=seq_t_mfill_errors)
		#ax2.errorbar(omp_mat_sizes, omp_t_mfill_means, yerr=omp_t_mfill_errors)

		ax3.set_title("DGEMM", loc="left")
		ax3.set_ylabel('Time (s)')
		ax3.set_xlabel('Matrix Size', fontsize=12)
		seq = ax3.errorbar(seq_mat_sizes, seq_t_gemm_means, yerr=seq_t_gemm_errors, label="_nolegend_", color=colors[i])
		#omp = ax3.errorbar(omp_mat_sizes, omp_t_gemm_means, yerr=omp_t_gemm_errors, label="_nolegend_", color="blue")

		seq2 = ax3.plot(seq_mat_sizes, seq_t_gemm_means, label=performance_versions[i], color=colors[i])
		#omp2 = ax3.plot(omp_mat_sizes, omp_t_gemm_means, label="Parallel OpenMP", color="blue")

		#ax1.legend(handles=[seq, omp], loc="upper left")
		ax3.legend(loc="upper left", fontsize=10)

		#plt.legend(handler_map={line1: HandlerLine2D(numpoints=4)})

	plt.savefig(plot_path)

	plt.close()

def plot_gemm(performance_versions, performance_data, plot_path):

	""" Plots all DGEMM timing segments: initialization, matrix filling, matrix multiplication, allowing comparison between different versions of the program

	Args:
		performance_versions = [performance_v1_name, performance_v2_name]
			performance_v1_name = "Sequential"
			performance_v2_name = "Parallel OpenMP"

		performance_data = [performance_v1, performance_v2, ... ]
			performance_v1 = [seq_mat_sizes, seq_t_gemm_means, seq_t_gemm_errors]
			performance_v2 = [omp_mat_sizes, omp_t_gemm_means, omp_t_gemm_errors]

		plot_path = "../data/performance.png"
	Returns:
		...
	Usage:
		plot_all([[seq_mat_sizes, seq_t_gemm_means, seq_t_gemm_errors],[omp_mat_sizes, omp_t_gemm_means, omp_t_gemm_errors]],"../data/performance.png")

	"""

	#
	# Plotting
	#
	# pd.read_csv(fpath,sep="\t") // DataFrame.from_csv

	#plot_path = "../data/performance.png"

	colors = ["b","g","r","c","m","y","k","w"]

	f, (ax1) = plt.subplots(1, sharex=False, sharey=False)

	plt.subplots_adjust(hspace=.4)

	for i,version in enumerate(performance_data):

		seq_mat_sizes 		= version[0]
		seq_t_gemm_means	= version[1]
		seq_t_gemm_errors	= version[2]

		ax1.set_title("DGEMM", loc="left")
		ax1.set_ylabel('Time (s)')
		ax1.set_xlabel('Matrix Size', fontsize=12)
		#seq = ax1.errorbar(seq_mat_sizes, seq_t_gemm_means, yerr=seq_t_gemm_errors, label="_nolegend_", color=colors[i], fmt='o', ms=1)
		seq = ax1.errorbar(seq_mat_sizes, seq_t_gemm_means, yerr=seq_t_gemm_errors, label=performance_versions[i], color=colors[i], fmt='o', ms=1)
		#seq2 = ax1.(seq_mat_sizes, seq_t_gemm_means, label=performance_versions[i], color=colors[i], lw=1, ms=1)

		ax1.legend(loc="upper left", fontsize=10)

		#plt.legend(handler_map={line1: HandlerLine2D(numpoints=4)})

	plt.savefig(plot_path)

	plt.close()



#
# 1) Sequenial code
#
seq_mat_sizes, seq_t_init_means, seq_t_mfill_means, seq_t_gemm_means, seq_t_init_errors, seq_t_mfill_errors, seq_t_gemm_errors = gen_performance_stats("../data/performance_seq_s2000_m5_v1.tsv")

#
# 2) Parallel Code: OpenMP 
#
omp_mat_sizes, omp_t_init_means, omp_t_mfill_means, omp_t_gemm_means, omp_t_init_errors, omp_t_mfill_errors, omp_t_gemm_errors = gen_performance_stats("../data/performance_omp_s2000_m5.tsv")

#
# Plot to compare performance: all time segments
#
plot_all(["Sequential", "Parallel OpenMP"], [[seq_mat_sizes, seq_t_init_means, seq_t_mfill_means, seq_t_gemm_means, seq_t_init_errors, seq_t_mfill_errors, seq_t_gemm_errors],[omp_mat_sizes, omp_t_init_means, omp_t_mfill_means, omp_t_gemm_means, omp_t_init_errors, omp_t_mfill_errors, omp_t_gemm_errors]], "../data/performance.png")

#
# Plot to compare performance: gemm only
#
plot_gemm(["Sequential", "Parallel OpenMP"], [[seq_mat_sizes, seq_t_gemm_means, seq_t_gemm_errors],[omp_mat_sizes, omp_t_gemm_means, omp_t_gemm_errors]],"../data/performance_gemm.png")



