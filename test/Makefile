testStatic: testStatic.cu
	nvcc -m64 -arch=sm_20 -o testStatic testStatic.cu

test: testStatic
	./testStatic
