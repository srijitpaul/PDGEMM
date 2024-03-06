#!/bin/bash
# Setup

maxsize=5000
samplesize=3

module purge; module load intel;  

echo "Running performance analysis..."
t=48

for ((threads=6; threads<=6; threads+=6)); do 

	#threads=4
	export OMP_NUM_THREADS=${threads}


	echo -e "\tNo. threads used: $OMP_NUM_THREADS"

	script="performance_blas${threads}_s${maxsize}_m${samplesize}.tsv"

	wd=${PWD##*/} 

	# Jump to correct directory
	if [ $wd != "src" ]; then
		echo "\t\tLooking for ParallelDGEMM..."
		target=`find ~ -type d -name "ParallelDGEMM"`
		cd $target/src
	elif [ $wd == "ParallelDGEMM" ]; then
		cd $target/src
	fi

	# Add a data directory if absent
	if [ -d "../data/" ]; then
		if [ -f "../data/$script" ]; then	
			echo -e "\t\tPerformance data found in ../data/$script, removing..."
			rm ../data/$script
		fi
		# if ![ -d "../data/tmp" ]; then
		# 	mkdir ../data/tmp/	
		# fi
	else
		echo -e "\t\tData directory NOT found, creating ../data/"
		mkdir ../data/
		mkdir ../data/tmp/
	fi

	# If no executable is available, try to compile
	if [ -f "./mat_blas.out" ]; then
		echo -e "\t\tExecutable found... ./mat_blas.out"
	else
		echo -e "\t\tExecutables NOT found, compiling..."
		mpicxx -std=c++11 -o mat_blas.out blas-matmul.cpp cblas-dgemm.cpp helpers.cpp matrix_t.cpp -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
	fi

	# performance
	echo -e "Matrix_Size\tSample No.\tInitialization_Time\tMatrix_Fill_Time\tGEMM_Time" >> ../data/$script

	for ((size=500; size <= ${maxsize}; size+=500)); do
		echo -e "\t\t matrix size = $size"
		for ((nsample=1; nsample <=${samplesize}; nsample++ )); do
			echo -e "\t\t\t sample = $nsample"
			printf '%s\t%s\t' "$size" "$nsample" >> ../data/$script
			#./mat_omp.out $size >> ../data/$script   
			OMP_NUM_THREADS=$threads srun -n 1 ./mat_blas.out $size >> ../data/$script
		done
	done
done

echo "Complete!"
echo -e "\t Find tab-delimited performance data in: ../data/$script"
