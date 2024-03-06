
# salloc --reservation=hpc-leap-5 --time=60 --nodes=1

# Set PWD
wd=${PWD##*/} 

if [ $wd != "src" ]; then
	echo "\t\tLooking for ParallelDGEMM..."
	target=`find ~ -type d -name "ParallelDGEMM"`
	cd $target/src
elif [ $wd == "ParallelDGEMM" ]; then
	cd $target/src
fi

# Compile if executable missing...
if [ -d ./mat_omp.out ]; then 
	echo "executable found: ./mat_omp.out: no need to compile..."
else
	echo "executable NOT found, compiling: ./mat_omp.out"
	g++ -std=c++11 -o mat_omp.out matmul.cpp gemm.cpp helpers.cpp matrix_t.cpp -fopenmp
fi

# Load

echo "Loading modules..."

module purge;module load GCC/5.1.0

echo "Now allocate resources: "
echo -e "\tsalloc --reservation=hpc-leap-5 --time=600 --nodes=1"

