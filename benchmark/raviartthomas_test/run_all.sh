make clean
make
for p in 0 1 2 3
do
./rt_test -p ${p}
done
