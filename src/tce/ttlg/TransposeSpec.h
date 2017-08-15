/*
 * TransposeSpec.h
 *
 *  Created on: Jan 30, 2017
 *      Author: arjun
 */

#ifndef TRANSPOSESPEC_H_
#define TRANSPOSESPEC_H_

#include "ourproj.h"


class TransposeSpec {
	int ndim;
	 int* sizes;
	 int* permutation;
	TensorType dataType;
	 double beta = 0.0;
	 double alpha = 1.0;

public:
	int* getSizes() {
		return sizes;
	}
	void setSizes(int* sizes) {
		this->sizes = sizes;
	}
	unsigned long int getVolume()
	{
		unsigned long vol = 1;
		for(int i = 0; i < ndim; i++)
		{
			vol *= sizes[i];
		}
		return vol;
	}

	int* getPermutation() {
		return permutation;
	}
	void setPermutation(int* permutation) {
		this->permutation = permutation;
	}
	TensorType getDataType() {
		return dataType;
	}
	void setDataType(TensorType dataType) {
		this->dataType = dataType;
	}
	double getBeta() {
		return beta;
	}
	void setBeta(double beta) {
		this->beta = beta;
	}
	TransposeSpec(int ndim, int *sizes, int *perm, TensorType dataType)
	{
		this->ndim = ndim;
		this->sizes = new int[ndim];
		this->permutation = new int[ndim];
		for(int i = 0; i < ndim; i++)
		{
			this->sizes[i] = sizes[i];
			this->permutation[i] = perm[i];
		}
		this->dataType = dataType;
	}
TransposeSpec() {
		// TODO Auto-generated constructor stub

	}


	//virtual
~TransposeSpec(){};

	int getNdim() const {
		return ndim;
	}

	void setNdim(int ndim) {
		this->ndim = ndim;
	}
double getAlpha() const
{
    return alpha;
}

void setAlpha(double alpha = 1.0)
{
    this->alpha = alpha;
}
}
;

#endif /* TRANSPOSESPEC_H_ */
