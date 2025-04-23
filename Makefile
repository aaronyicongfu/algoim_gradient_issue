ALGOIM_DIR = /Users/fyc/git/algoim/algoim
LAPACK_DIR = /opt/homebrew/Cellar/lapack/3.12.0
default:
	c++ -std=c++17 -I${ALGOIM_DIR} -I${LAPACK_DIR}/include -L${LAPACK_DIR}/lib -llapacke test_algoim_gradient_2d.cpp -o test && ./test

clean:
	rm -f test
