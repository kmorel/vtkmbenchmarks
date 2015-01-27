#bin/bash

./BenchmarkCuda --file=/temp/test.nhdr --pipeline=1 > cuda_threshold.csv
./BenchmarkCuda --file=/temp/test.nhdr --pipeline=2 > cuda_mc.csv

./BenchmarkSerial --file=/temp/test.nhdr --pipeline=1 > serial_threshold.csv
./BenchmarkSerial --file=/temp/test.nhdr --pipeline=2 > serial_mc.csv
