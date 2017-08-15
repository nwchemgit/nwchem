/*
 * Parameters.h
 *
 *  Created on: Jan 30, 2017
 *      Author: arjun
 */
#include <string>
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

class Parameters {

	 float warpEfficiency;
	 float occupancy;
	 int sharedMemSize1;
	 int sharedMemSize2;
	 int tbSize;
	 int tileSize;
	 int tileSize1;
	 int numBlocksPerSM;
	 int paddingSize;
		int blockAIndex;
		int blockBIndex;
	 int numElementsProcessedPerBlock;
	 int numElementsProcessedPerBlock1;
	 bool isParetoOptimal;
	 int remElements;
double bandwidth;
double etime;
	public:
	int getRemElements() {
		return remElements;
	}

	void setRemElements(int remElements) {
		this->remElements = remElements;
	}

	bool getIsParetoOptimal() {
		return isParetoOptimal;
	}

	void setIsParetoOptimal(bool isParetoOptimal) {
		this->isParetoOptimal = isParetoOptimal;
	}

	int getSharedMemSize2() {
		return sharedMemSize2;
	}

	void setSharedMemSize2(int sharedMemSize2) {
		this->sharedMemSize2 = sharedMemSize2;
	}

	int getNumElementsProcessedPerBlock() {
		return numElementsProcessedPerBlock;
	}

	void setNumElementsProcessedPerBlock(int numElementsProcessedPerBlock) {
		this->numElementsProcessedPerBlock = numElementsProcessedPerBlock;
	}
	int getNumElementsProcessedPerBlock1() {
		return numElementsProcessedPerBlock1;
	}

	void setNumElementsProcessedPerBlock1(int numElementsProcessedPerBlock1) {
		this->numElementsProcessedPerBlock1 = numElementsProcessedPerBlock1;
	}

	int getPaddingSize() {
		return paddingSize;
	}

	void setPaddingSize(int paddingSize) {
		this->paddingSize = paddingSize;
	}
	int getBlockAIndex() {
		return blockAIndex;
	}

	void setBlockAIndex(int blockAIndex) {
		this->blockAIndex = blockAIndex;
	}
	int getBlockBIndex() {
		return blockBIndex;
	}

	void setBlockBIndex(int blockBIndex) {
		this->blockBIndex = blockBIndex;
	}

	 float getWarpEfficiency() {
		return warpEfficiency;
	}

	void setWarpEfficiency(float warpEfficiency) {
		this->warpEfficiency = warpEfficiency;
	}

	 float getOccupancy() {
		return occupancy;
	}

	void setOccupancy(float occupancy) {
		this->occupancy = occupancy;
	}

	 /*string getSharedMemSize() {
		if(sharedMemSize2 != 0)
			return sharedMemSize1+"]["+sharedMemSize2;
		else
			return sharedMemSize1+"";
	}
*/
	void setSharedMemSize1(int sharedMemSize) {
		this->sharedMemSize1 = sharedMemSize;
	}


	int getTbSize() {
		return tbSize;
	}

	void setTbSize(int tbSize) {
		this->tbSize = tbSize;
	}

	int getTileSize() {
		return tileSize;
	}
	int getTileSize1() {
		return tileSize1;
	}

	int getSharedMemSize1() {
		return sharedMemSize1;
	}

	void setTileSize(int tileSize) {
		this->tileSize = tileSize;
	}
	void setTileSize1(int tileSize) {
		this->tileSize1 = tileSize;
	}

	int getNumBlocksPerSM() {
		return numBlocksPerSM;
	}
	double getBW() {
		return bandwidth;
	}
	void setBW(double bw) {
		bandwidth = bw;
	}
	double getTime() {
		return etime;
	}
	void setTime(double esttime) {
		etime = esttime;
	}

	void setNumBlocksPerSM(int numBlocksPerSM) {
		this->numBlocksPerSM = numBlocksPerSM;
	}


	int hashCode() {
		const int prime = 31;
		int result = 1;
		return result;
	/*	result = prime * result + numBlocksPerSM;
		result = prime * result + Float.floatToIntBits(occupancy);
		result = prime * result + paddingSize;
		result = prime * result + sharedMemSize1;
		result = prime * result + tbSize;
		result = prime * result + tileSize;
		result = prime * result + Float.floatToIntBits(warpEfficiency);
		return result;*/
	}

/*
	 bool equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Parameters other = (Parameters) obj;
		if (numBlocksPerSM != other.numBlocksPerSM)
			return false;
		if (Float.floatToIntBits(occupancy) != Float.floatToIntBits(other.occupancy))
			return false;
		if (paddingSize != other.paddingSize)
			return false;
		if (sharedMemSize1 != other.sharedMemSize1)
			return false;
		if (tbSize != other.tbSize)
			return false;
		if (tileSize != other.tileSize)
			return false;
		if (Float.floatToIntBits(warpEfficiency) != Float.floatToIntBits(other.warpEfficiency))
			return false;
		return true;
	}
*/
/*	string print() {
		int shm1 =1, shm2 =1;
		if(sharedMemSize1 != 0)
			shm1=sharedMemSize1;
		else if(sharedMemSize2 != 0 )
			shm2 =sharedMemSize2;
		return "Parameters [warpEfficiency=" + warpEfficiency + ", occupancy=" + occupancy + ", sharedMemSize="
		+ shm1*shm2 + ", tbSize=" + tbSize + ", tileSize=" + tileSize + ", numBlocksPerSM=" + numBlocksPerSM
		+ ", paddingSize=" + paddingSize +", RatioOfWork-CompleteToRemainder=" + numElementsProcessedPerBlock+":"+remElements + ", IsParetoOptimal=" + isParetoOptimal+"]";

	}


	string toString() {

		int shm1 =1, shm2 =1;
		if(sharedMemSize1 != 0)
			shm1=sharedMemSize1;
		else if(sharedMemSize2 != 0 )
			shm2 =sharedMemSize2;

		return warpEfficiency +" "+ occupancy + " "+ (shm1*shm2 )+" "+tbSize + " " + tileSize + " " + numBlocksPerSM
				+ " " + paddingSize + " " + numElementsProcessedPerBlock+":"+remElements + " " + isParetoOptimal;

	}*/


	Parameters();
	virtual ~Parameters();
};

#endif /* PARAMETERS_H_ */
