#!/bin/bash
# Setup

exe="mat_mpi_prag1_simd.out"
prefix1="mpi"  # e.g. "mpi" in "performance_<<mpi>>${tasks}_omp${threads}_s${maxsize}_m${samplesize}.tsv"
prefix2="simd" # e.g. "omp" in "performance_mpi${tasks}_<<omp>>${threads}_s${maxsize}_m${samplesize}.tsv"

samplesize=3
#maxsize=5040
minsize=500

#mintasks=16  
#maxtasks=4
cpu_x_tasks=96  # i.e. --cpus-per-task x --ntasks
#maxcpus=12
threads=1

for square in 16 25 36 49 64 81
do
	mintasks=${square} # $tasks must be a square number, so each iteration of $tasks is: $tasks, $tasks^2, $tasks^2^2, $tasks^2^2^2, etc. 

	#
	# Sanity checks
	#
	sqrt_mintasks=`echo "sqrt (${mintasks})" | bc`
	
	if [ "$(($sqrt_mintasks * $sqrt_mintasks))" -ne "$mintasks" ]; then
		echo "ERROR: the mintasks variable must be a square number!"
		exit	
	fi
	
	#
	# Determine sequence of matrix size parameterizations
	#
	
	echo "Calculating valid matrix size parameterizations..."
	while [ "$((${minsize} % ${sqrt_mintasks} ))" -ne "0" ]; do
	        minsize=$(($minsize + 1))
	done
	maxsize=$(($minsize*10))
	
	echo "Loading modules..."
	module purge; module load GCC; module load iimpi 
	
	#
	# Main Loop for running performance benchmarks
	#
	
	echo "Running performance analysis..."
	
	for ((tasks=${mintasks}; $((${threads}*${tasks})) <= ${cpu_x_tasks}; tasks=$((${tasks}*${tasks})) )); do
	
		echo -e "\tNo. TASKS used: $tasks"
	
		for ((threads=${threads}; $((${threads}*${tasks})) <= ${cpu_x_tasks}; threads+=1)); do 
		
			export OMP_NUM_THREADS=${threads}
		
			echo -e "\t\tNo. THREADS used: $OMP_NUM_THREADS"
		
			script="performance_$prefix1${tasks}_$prefix2${threads}_s${maxsize}_m${samplesize}.tsv" 
			#^ e.g. performance_mpi16_omp24_s5120_m3.tsv
			# i.e. 16 mpi processes, 24 OpenMP threads, maximum matrix size of 5120 and repitiions of 3 per matrix size
		
			wd=${PWD##*/} 
		
			# Jump to correct directory
			if [ $wd != "src" ]; then
				echo "\t\t\tLooking for ParallelDGEMM..."
				target=`find ~ -type d -name "ParallelDGEMM"`
				cd $target/src
			elif [ $wd == "ParallelDGEMM" ]; then
				cd $target/src
			fi
		
			# Add a data directory if absent
			if [ -d "../data/" ]; then
				if [ -f "../data/$script" ]; then	
					echo -e "\t\t\tPerformance data found in ../data/$script, removing..."
					rm ../data/$script
				else
					echo -e "\t\t\tPerformance data NOT found, creating ../data/$script"
				fi
			else
				echo -e "\t\t\tData directory NOT found, creating ../data/$script"
				mkdir ../data/
				mkdir ../data/tmp/
			fi
		
			# If no executable is available, try to compile
			if [ -f ./$exe ]; then
				echo -e "\t\t\tExecutable found... ./$exe"
			else
				echo -e "\t\t\tExecutables NOT found, compiling... "
				mpicxx -std=c++11 -o ./$exe mpi-matmul-timed.cpp gemm.cpp mpi-helpers.cpp matrix_t.cpp -fopenmp -O3
			fi
		
			# performance
			echo -e "Matrix_Size\tSample No.\tInitialization_Time\tMatrix_Fill_Time\tGEMM_Time" >> ../data/$script
		
			for ((size=${minsize}; size <= ${maxsize}; size+=${minsize})); do
				echo -e "\t\t\t matrix size = $size"
				for ((nsample=1; nsample <=${samplesize}; nsample++ )); do
					echo -e "\t\t\t\t sample = $nsample"
					printf '%s\t%s\t' "$size" "$nsample" >> ../data/$script
					#./mat_omp.out $size >> ../data/$script   
					OMP_NUM_THREADS=$threads srun -n ${tasks} --cpus-per-task=${threads} ./$exe $size >> ../data/$script
				done
			done
		done
	
	done
	
	echo "Complete!"
	echo -e "\t Find tab-delimited performance data in: ../data/$script"
done

