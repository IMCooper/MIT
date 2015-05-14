# Deal generated:
# cube:
for p in 0 1 2
do
./polynomial_benchmark -p ${p} -i ../input_files/dealcube.prm -o test
done

# distorted cube:
for p in 0 1 2
do
./polynomial_benchmark -p ${p} -i ../input_files/dealcube_distorted.prm -o test
done

# cylinder:
for p in 0 1 2
do
./polynomial_benchmark -p ${p} -i ../input_files/dealcyl.prm -o test
done

# sphere
for p in 0 1 2
do
./polynomial_benchmark -p ${p} -i ../input_files/dealsphere.prm -o test
done

# CUBIT generated:
# cube:
for p in 0 1 2
do
./polynomial_benchmark -p ${p} -i ../input_files/cubitcube.prm -o test
done

# cylinder
for p in 0 1 2
do
./polynomial_benchmark -p ${p} -i ../input_files/cubitcyl.prm -o test
done
