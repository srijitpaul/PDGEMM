#!/bin/bash
# Setup

maxsize=5000
samplesize=3
t=48
exe="./mat_prag1_simd.out"

module purge; module load GCC 

echo "Running performance analysis..."

for ((threads=1; threads<=$t; threads+=1)); do 

	export OMP_NUM_THREADS=${threads}

	echo -e "\tNo. threads used: $OMP_NUM_THREADS"

	script="performance_prag1_simd${threads}_s${maxsize}_m${samplesize}.tsv"

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
	if [ -f $exe ]; then
		echo -e "\t\tExecutable found... ./$exe"
	else
		echo -e "\t\tExecutables NOT found, compiling..."
		g++ -std=c++11 -o $exe matmul.cpp gemm.cpp helpers.cpp matrix_t.cpp -fopenmp -O3
	fi

	# performance
	echo -e "Matrix_Size\tSample No.\tInitialization_Time\tMatrix_Fill_Time\tGEMM_Time" >> ../data/$script

	for ((size=5000; size <= ${maxsize}; size+=500)); do
		echo -e "\t\t matrix size = $size"
		for ((nsample=1; nsample <=${samplesize}; nsample++ )); do
			echo -e "\t\t\t sample = $nsample"
			printf '%s\t%s\t' "$size" "$nsample" >> ../data/$script
			#./mat_omp.out $size >> ../data/$script   
			OMP_NUM_THREADS=$threads srun -n 1 ./$exe $size >> ../data/$script
		done
	done
done

echo "Complete!"
echo -e "\t Find tab-delimited performance data in: ../data/$script"
