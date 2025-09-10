#include <stdio.h>
#include <algorithm>
#include <math.h>
#include "CMU418intrin.h"
#include "logger.h"
using namespace std;


void absSerial(float* values, float* output, int N) {
    for (int i=0; i<N; i++) {
	float x = values[i];
	if (x < 0) {
	    output[i] = -x;
	} else {
	    output[i] = x;
	}
    }
}

// implementation of absolute value using 15418 instrinsics
void absVector(float* values, float* output, int N) {
    __cmu418_vec_float x;
    __cmu418_vec_float result;
    __cmu418_vec_float zero = _cmu418_vset_float(0.f);
    __cmu418_mask maskAll, maskIsNegative, maskIsNotNegative;

    //  Note: Take a careful look at this loop indexing.  This example
    //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
    //  Why is that the case?
    for (int i=0; i<N; i+=VECTOR_WIDTH) {

	// All ones
	maskAll = _cmu418_init_ones();

	// All zeros
	maskIsNegative = _cmu418_init_ones(0);

	// Load vector of values from contiguous memory addresses
	_cmu418_vload_float(x, values+i, maskAll);               // x = values[i];

	// Set mask according to predicate
	_cmu418_vlt_float(maskIsNegative, x, zero, maskAll);     // if (x < 0) {

	// Execute instruction using mask ("if" clause)
	_cmu418_vsub_float(result, zero, x, maskIsNegative);      //   output[i] = -x;

	// Inverse maskIsNegative to generate "else" mask
	maskIsNotNegative = _cmu418_mask_not(maskIsNegative);     // } else {

	// Execute instruction ("else" clause)
	_cmu418_vload_float(result, values+i, maskIsNotNegative); //   output[i] = x; }

	// Write results back to memory
	_cmu418_vstore_float(output+i, result, maskAll);
    }
}

// Accepts an array of values and an array of exponents
// For each element, compute values[i]^exponents[i] and clamp value to
// 4.18.  Store result in outputs.
// Uses iterative squaring, so that total iterations is proportional
// to the log_2 of the exponent
void clampedExpSerial(float* values, int* exponents, float* output, int N) {
    for (int i=0; i<N; i++) {
	float x = values[i];
	float result = 1.f;
	int y = exponents[i];
	float xpower = x;
	while (y > 0) {
	    if (y & 0x1) {
			result *= xpower;
		}
	    xpower = xpower * xpower;
	    y >>= 1;
	}
	if (result > 4.18f) {
	    result = 4.18f;
	}
	output[i] = result;
    }
}

void clampedExpVector(float* values, int* exponents, float* output, int N) {
    // Implement your vectorized version of clampedExpSerial here
    //  ...

    __cmu418_vec_float x;
    __cmu418_vec_int y;
    __cmu418_vec_float result;
    __cmu418_vec_float xpower;
    __cmu418_vec_float ones = _cmu418_vset_float(1.f);
    __cmu418_vec_float clamp_val = _cmu418_vset_float(4.18f);
    __cmu418_vec_int int_ones = _cmu418_vset_int(1);
    __cmu418_vec_int zeros = _cmu418_vset_int(0);
    __cmu418_mask maskAll, maskActive, maskOddBit, maskGreaterThanClamp;

    for (int i = 0; i < N; i += VECTOR_WIDTH) {
        int remaining = N - i;
        if (remaining >= VECTOR_WIDTH) {
            maskAll = _cmu418_init_ones();
        } else {
            maskAll = _cmu418_init_ones(remaining);
        }
        
        _cmu418_vload_float(x, values + i, maskAll);
        _cmu418_vload_int(y, exponents + i, maskAll);
        
        result = _cmu418_vset_float(1.f);
        xpower = x;
        
        _cmu418_vgt_int(maskActive, y, zeros, maskAll);
        
        while (_cmu418_cntbits(maskActive) > 0) {
            maskOddBit = _cmu418_init_ones(0);
            
            __cmu418_vec_int temp_and;
            _cmu418_vbitand_int(temp_and, y, int_ones, maskAll);
            _cmu418_vgt_int(maskOddBit, temp_and, zeros, maskAll);
            
            maskOddBit = _cmu418_mask_and(maskOddBit, maskActive);
            
            _cmu418_vmult_float(result, result, xpower, maskOddBit);
            
            _cmu418_vmult_float(xpower, xpower, xpower, maskActive);
            
            _cmu418_vshiftright_int(y, y, int_ones, maskActive);
            
            _cmu418_vgt_int(maskActive, y, zeros, maskAll);
        }
        
        _cmu418_vgt_float(maskGreaterThanClamp, result, clamp_val, maskAll);
        _cmu418_vset_float(result, 4.18f, maskGreaterThanClamp);
        
        _cmu418_vstore_float(output + i, result, maskAll);
    }

}


float arraySumSerial(float* values, int N) {
    float sum = 0;
    for (int i=0; i<N; i++) {
	sum += values[i];
    }

    return sum;
}

// Assume N % VECTOR_WIDTH == 0
// Assume VECTOR_WIDTH is a power of 2
float arraySumVector(float* values, int N) {
	// Implement your vectorized version here
    //  ...
    __cmu418_vec_float vec_sum = _cmu418_vset_float(0.f);
    __cmu418_vec_float vec_values;
    __cmu418_mask maskAll = _cmu418_init_ones();
    
    for (int i = 0; i < N; i += VECTOR_WIDTH) {
        _cmu418_vload_float(vec_values, values + i, maskAll);
        _cmu418_vadd_float(vec_sum, vec_sum, vec_values, maskAll);
    }
    
    float result = 0.f;
    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result += vec_sum.value[i];
    }
    
    return result;
}
