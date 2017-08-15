/*
 * ParameterTuner.cpp
 *
 *  Created on: Jan 30, 2017
 */

#include "ParameterTuner.h"
#include "ourproj.h"
#include <math.h>
#include "Parameters.h"
#include "BlockingCase.h"
#include "TransposeSpec.h"
#include <iostream>

#include "ourinclude.h"
using namespace std;

class ParameterTuner {
	int sharedMemLimitPerSM;// = 6144; // (48 * 1024)/8 words per SM
	int numThreadsLimitPerSM;// =  2048; //per SM
	int numThreadBlocksLimitPerSM;// =  16; //per SM
	int threadBlocKSizeLimit;// = 1024;
	int numSMs;// = 15;
	int blockFactor;// = 4
	BlockingCase caseId;

	public:
	unsigned int getTBSize(unsigned int shm)
	{
		return (unsigned) ceil((double)(numThreadsLimitPerSM/floor(32.0*(unsigned)(sharedMemLimitPerSM/ shm)))) * 32;
	}
	BlockingCase getCaseId()
	{
		return caseId;
	}
	Parameters& tune(TransposeSpec &spec) {
		int *sizes = spec.getSizes();
		if((spec.getPermutation()[0] != 0) || (sizes[0] * sizes[1] < 32) || (sizes[spec.getPermutation()[0]] * sizes[spec.getPermutation()[1]] < 32) )
		{
			return tuneFastestVaryingNotMatchingCase(spec);
		}
		else if ( spec.getSizes()[0] >= 32)
		{	
			caseId = BlockingCase::FVI_MATCH_AND_GREATERT32;
			return tuneFastestVaryingMatchingCaseWithoutBlocking(spec);
		}
		else
		{
			caseId = BlockingCase::FVI_MATCH_AND_LESST32;
			return tuneFastestVaryingMatchingCaseWithBlocking(spec);
		}
	}

	Parameters& tuneFastestVaryingNotMatchG32Case(TransposeSpec &spec){
		Parameters *parameters = new Parameters();
		int *sizes = spec.getSizes();
		int* permutation = spec.getPermutation();
		unsigned numElements = 32*32;
		unsigned sharedMemSize = (32 *33);
		unsigned paddingSize = 1;
		unsigned tbSize = 256;
		unsigned numBlocksPerSM = sharedMemLimitPerSM/sharedMemSize;
		parameters->setNumElementsProcessedPerBlock(numElements);
		parameters->setPaddingSize(paddingSize);
		parameters->setOccupancy(100.0f);
		parameters->setTbSize(tbSize);
		parameters->setNumBlocksPerSM(numBlocksPerSM);
		parameters->setTileSize(32);
		parameters->setSharedMemSize1(32);
		parameters->setSharedMemSize2(33);
		double eff = getEfficiency_nomatchg32(sizes[0], sizes[permutation[0]]);
		double bandwidth = getBW_nomatchg32(eff);
		parameters -> setBW(bandwidth);
		unsigned long vol = spec.getVolume();

		//cout <<"cc0= "<< spec.getVolume() << " eff = "<<eff<<" bw = "<<bandwidth<< " time = "<<getTime(bandwidth, spec.getVolume())<<"\t";
		//cout <<"cc1= "<< vol << " eff = "<<eff<<" bw = "<<bandwidth<< " time = "<<getTime(bandwidth, vol)<<"\t";
		parameters -> setTime(getTime(bandwidth, spec.getVolume()));
		parameters->setWarpEfficiency(eff);
		return *parameters;
		//parametersList.add(parameters);
	}
	Parameters& tuneConflictCase(TransposeSpec &spec){
		Parameters *parameters = new Parameters();
		int *sizes = spec.getSizes();
		int* permutation = spec.getPermutation();
		int blockA = 1, blockB = 1;
		int sharedMemSize = 1;
		int sharedMemSize1 = 1;
		int sharedMemSize2 = 1;
		int tbSize = 32;
		int numElements = 0;
		int paddingSize= 0;
		int csize, asize, bsize, bonlysize, pad;//sizes[0];
		int repeat = 0, rlimit, alimit, blimit,count = 0;
		unsigned SHMLIMIT = 1056;
		//unsigned SHMLIMIT = 1400;
		const int limit = 32;//starts from 32 and goes 64, 128...
		int limiti, limito, nlimit;
		unsigned long int volume = spec.getVolume(); 
		int minnumblocks = numSMs * sharedMemLimitPerSM/(33*32); 
		nlimit = sqrt(volume/(blockFactor * minnumblocks* 32*32)); 
		if(nlimit == 0) nlimit = 1;
		double besteff = 0;
		for(int limiti = 0; limiti < nlimit; limiti++)
		{

			rlimit = 32 + 32*limiti;
			int i;
			blockA = 1, blockB = 1, csize = 1;
			for(i = 0; i <= spec.getNdim(); i++)
			{
				if(csize == rlimit)
				{
					break;
				}
				if(csize > rlimit)
				{
					csize/=sizes[i-1];
					//blockA = (rlimit+csize-1)/csize;
						blockA = (rlimit)/csize;
					//cyyout << "blockA = "<<blockA;
					if(blockA != 1){
						//csize/=sizes[i-1];
						csize*=blockA;
					}
					else i--;
					if(blockA == sizes[i-1]) { blockA = 1;}
					break;
				}
			if(i < spec.getNdim())	
				csize*= sizes[i];
			}
				if(i == spec.getNdim() + 1)
					i--;
			//if(blockA == 1 && i < spec.getNdim()) i++;
			alimit = i-1;
			asize = csize;
#ifdef printd
			cout <<  "asize == "<<asize<<" ablock = "<<blockA<<" rlimit= "<<rlimit<<"\n";
#endif
			for(int limito = 0; limito < nlimit; limito++)
			{
				i = 0;
				bonlysize = 1;
				int	limit = 32 + limito * 32;//1024/(asize*2);
				//int	limit = 32;//rlimit;//1024/(asize*2);
				csize = 1;//sizes[permutation[i]];

				for(; i <= spec.getNdim(); i++)
				{
					if(csize == limit)
					{
						break;
					}
					if(csize > limit)
					{
						if(i > 0){
							csize/=sizes[permutation[i-1]];
						//	blockB =  (limit+csize-1)/csize;
								blockB =  (limit)/csize;
#ifdef printd
							cout << "\ncsize = "<<csize;
							cout << "\nblockB = "<<blockB;
#endif
							if( permutation[i-1] > alimit) 
							{
								bonlysize /= sizes[permutation[i-1]];
								if(blockB != 1){
									bonlysize*= blockB;

								}
							}
							if(blockB != 1){
								csize*= blockB;
								if(blockB == sizes[permutation[i-1]]) {blockB = 1;}

							}
							else i--;


						}

						break;
					}	
					if(i == spec.getNdim()) break;
					//if(permutation[i] < alimit) continue;
					//if(i > 0 && ((blockA != 1) && (permutation[i] == alimit)))
					if(i < spec.getNdim())	
					csize*= sizes[permutation[i]];
					if(permutation[i] > alimit)
					{
						bonlysize *= sizes[permutation[i]];
					} 
#ifdef printd
					cout<<"\ni = "<<i<<"bsize = "<<csize<<"\n";
#endif
				}
#ifdef printd
				cout <<  "\nbsize == "<<csize<<"\n";
#endif
				if(i == spec.getNdim() + 1)
					i--;
				//if((blockB == 1) && (i < spec.getNdim())) i++;
				blimit = i-1;
				//cout<<" alimitp = "<<alimit<<" blimitp = "<<blimit<<"\n";
				bsize = csize;
				int n = spec.getNdim();
				int rperm[20];
				for(int i = 0; i < n; i++)
				{
					for(int j = 0; j < n; j++)
					{
						if(permutation[i] == j)
						{
							rperm[j] = i;
						}
					}
				}



				if(blockA > 1)//checking for inner dimensions in output which gets blocked in input
				{
					for(int i = 0; i < blimit; i++)
					{
						if(permutation[i] == alimit)
						{
							asize /= blockA;
							asize *= sizes[permutation[i]]; 
							blockA = 1;
						}
					}
				}
				//	cout <<" blockA = "<<blockA <<" blockB = "<<blockB;
				if(blockB > 1)//checking for inner dimensions in input which gets blocked in output
				{

					//	cout <<"A smaller in B\n";
					for(int i = 0; i < alimit; i++)
					{
						if(rperm[i] == blimit)
						{
							bsize /= blockB;
							bsize *= sizes[i]; 
							//	cout <<"\nNew bsize "<<bsize<<"\n";
							blockB = 1;
						}
					}
				}

				//we need to change the blocksize in case alimit and blimit are same dimension
				if((alimit == permutation[blimit]) && (blockA > 1 || blockB > 1))
				{
#ifdef printd
					cout <<"Blocking dimensions same, bloackA = "<<blockA<<" blockB = "<<blockB<<"\n"; 
#endif
					if(((blockA > blockB) && (blockB != 1)) || (blockA == 1))// > bsize)
					{
						if(blockA == 1)
							bsize = (bsize/blockB)*sizes[alimit], blockB = blockA;
						else
							bsize = (bsize/blockB)*blockA, blockB = blockA;

					}
					else 
					{
						if(blockB == 1)
							asize = (asize/blockA) * sizes[permutation[blimit]], blockA = blockB;
						else
							asize = (asize/blockA) * blockB, blockA = blockB;
					}
				}
				pad =  ((asize %2)+1)%2;
				sharedMemSize1 = asize + pad;
				sharedMemSize2 = bonlysize;
				sharedMemSize = sharedMemSize1 * sharedMemSize2;//csize;//blockB * sizes[permutation[0]];
				//tbSize = max(64, (sharedMemSize/64)*32);
				//tbSize = 320;// max(64, (sharedMemSize/64)*32);
				tbSize = getTBSize(sharedMemSize);
				if(sharedMemSize > SHMLIMIT && besteff > 0)
				{
				//	limiti = nlimit;
				       	break;
				}	
				double 	eff = getEfficiency_overlap(asize, bsize, sizes[alimit], sizes[permutation[blimit]], blockA, blockB);
				//cout <<"\t"<<eff<<"\t"<<asize<<"\t"<<bsize<<"\n";
				if(eff > besteff)
				{
					parameters->setNumElementsProcessedPerBlock(asize);
					parameters->setNumElementsProcessedPerBlock1(bsize);
					parameters->setPaddingSize(rlimit);
					parameters->setTbSize(tbSize);
					parameters->setTileSize(blockA);
					parameters->setSharedMemSize1(sharedMemSize1);
					parameters->setSharedMemSize2(sharedMemSize2);
					parameters->setTileSize1(blockB);
					parameters -> setBlockAIndex(alimit);
					parameters -> setBlockBIndex(blimit);
					besteff = eff;

				}

#ifdef printd
				cout<<" alimit = "<<alimit<<" blimit = "<<blimit<<"\n";
				cout<<" asize = "<<asize<<" bsize = "<<bsize<<"\n";
				cout<<" blockA = "<<blockA<<" blockB = "<<blockB<<"\n";
				cout<<" SM1 = "<<sharedMemSize<<"\n";// blockB = "<<blockB<<"\n";
				cout <<" TBSize = "<<tbSize<<"\n";
#endif
				count++;
				repeat++;
			}
		}// while(sharedMemSize <= 512);
		//double eff = getEfficiency_overlap(asize, bsize, sizes[alimit], sizes[permutation[blimit]], blockA, blockB );
		double bandwidth = getBW_overlap(besteff);
		parameters -> setBW(bandwidth);
		parameters -> setTime(getTime(bandwidth, spec.getVolume()));
		parameters->setWarpEfficiency(besteff);
		return *parameters;

	}
	Parameters& tuneFastestVaryingNotMatchingCase(TransposeSpec &spec){
		int *sizes = spec.getSizes();
		int* permutation = spec.getPermutation();
		TensorType mytype = spec.getDataType();
		int blockAIndex = 0;
		int blockBIndex = 0;
		Parameters &parameters1 = tuneFastestVaryingNotMatchG32Case(spec);
		Parameters &parameters2 = tuneNonConflictCase(spec);
		Parameters &parameters3 = tuneConflictCase(spec);
		double b1, b2, b3;
		b1 = getBW_nomatchg32(parameters1 . getWarpEfficiency());
		//	if(b1 < 0.7)
		//		parameters1.setTbSize(256);
		double eff2 = parameters2 . getWarpEfficiency();
		//if(eff2 < 0.7) eff2*=0.9;
		b2 = getBW_nooverlap(eff2);
		b3 = getBW_overlap(parameters3 . getWarpEfficiency());
		//if(parameters1 . getWarpEfficiency() > parameters2 . getWarpEfficiency())
#ifdef printd
		cout<<"\t"<<parameters1 . getWarpEfficiency()<<"\t"<<parameters2 . getWarpEfficiency()<<"\t"<<parameters3 . getWarpEfficiency();
		cout <<"\t"<<b1<<"\t"<<b2<<"\t"<<b3<<"\t";
#endif

		if(b1 >= b2)
		{
			//if(parameters1 . getWarpEfficiency() > parameters3 . getWarpEfficiency())
			if(b1 >= b3)
			{
				caseId = BlockingCase::FVI_NOMATCH_AND_GREATERT32;
				return parameters1;
			}
			else
			{
				caseId = BlockingCase::FVI_NOMATCH_GENERAL_OVERLAP;
				return parameters3;
			}
		}
		else if(b2 >= b3)
		{
			caseId = BlockingCase::FVI_NOMATCH_GENERAL;
			return parameters2;
		}
		else
		{
			caseId = BlockingCase::FVI_NOMATCH_GENERAL_OVERLAP;
			return parameters3;
		}


	}
	Parameters& tuneNonConflictCase(TransposeSpec &spec){
		Parameters *parameters = new Parameters();
		int *sizes = spec.getSizes();
		int* permutation = spec.getPermutation();
		int sharedMemSize = 1;
		int sharedMemSize1 = 1;
		int sharedMemSize2 = 1;
		int numElements = 0;
		int paddingSize= 0;
		int n = spec.getNdim();
		int rperm[20];
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				if(permutation[i] == j)
				{
					rperm[j] = i;
				}
			}
		}
		int csize = 1, asize, bsize;//sizes[0];
		int rlimit = 32;
		//int tbSize = 512;
		int tbSize = 352;
		int alimit = 0, blimit = 0, irlimit = 0;
		int i, blockA = 1, blockB = 1;
		irlimit += 32;
		int limiti, limito, nlimit;
		unsigned long int volume = spec.getVolume(); 
		int minnumblocks = numSMs * sharedMemLimitPerSM/(33*32); 
		nlimit = sqrt(volume/(blockFactor * minnumblocks* 32*32)); 
		//nlimit = 4;
		if(nlimit == 0) nlimit = 1;
#ifdef printd
		cout <<"\nnlimit = "<<nlimit;
#endif
		double besteff = 0;

		for(limiti = 0; limiti <  nlimit; limiti++)
		{
			int limitir;
				limitir	= 32+ limiti*32;	
			for(limito = 0; limito < nlimit; limito++)
			{
				int limitor;
					limitor = 32+ limito*32;	
				bool conflict = false;
				csize = 1, blockA = 1, blockB = 1;
				for(i = 0; i < spec.getNdim(); i++)
				{
					if(csize == limitir)
					{
						break;
					}
					if(csize > limitir)
					{
						csize/=sizes[i-1];
						//blockA = (limitir+csize-1)/csize;
						blockA = (limitir)/csize;
						//cyyout << "blockA = "<<blockA;
						if(blockA != 1){
							//csize/=sizes[i-1];
							csize*=blockA;
						}
						else i--;
						if(blockA == sizes[i-1]) { blockA = 1;}
						break;
					}	
					csize*= sizes[i];
				}
				if(i == spec.getNdim() + 1)
					i--;
				//if(blockA == 1 && i < spec.getNdim()) i++;
				alimit = i-1;
				asize = csize;
				i = 0;
				int conflicti = -1;
				//rlimit = 64;
				csize = 1;//sizes[permutation[i]];
				for(; i <= spec.getNdim(); i++)
				{
					if (i > 0 && permutation[i-1] < alimit)
					{
						limito = nlimit;
						conflict = true;
#ifdef printd
						cout <<"caseid changed "<< i<<"\n";
#endif
						//caseId = BlockingCase::FVI_NOMATCH_GENERAL_OVERLAP;
						//break;

					}
					if(csize == limitor)
					{
						break;
					}
					if(csize> limitor)
					{

						csize/=sizes[permutation[i-1]];
						//blockB =  (limitor+csize-1)/csize;
						blockB =  (limitor)/csize;
#ifdef printd
						cout << "\ncsize = "<<csize;
						cout << "\nblockB = "<<blockB;
#endif
						if(blockB != 1){
							//if (permutation[i-1] < alimit || ((blockA != 1) && (permutation[i] == alimit)))
							if (i > 0 && permutation[i-1] <= alimit )
							{
								limito = nlimit;
								conflict = true;
#ifdef printd
								cout <<"caseid changed "<< i<<"\n";
#endif
								//		caseId = BlockingCase::FVI_NOMATCH_GENERAL_OVERLAP;
								//		break;
							}

							//csize/=sizes[permutation[i-1]];
							csize*= blockB;
							if(blockB == sizes[permutation[i-1]]) { blockB = 1;}

						}
						else{ i--;
							if(conflicti == i)
								conflict = false;
						}
						break;
					}	
					if(i == spec.getNdim()) break;
					//if(permutation[i] < alimit) continue;
					//if(i > 0 && ((blockA != 1) && (permutation[i] == alimit)))
					if((permutation[i] <= alimit))
					{
						limito = nlimit;
#ifdef printd
						cout <<"caseid changed "<< i<<"\n";
#endif
						conflict = true;
						if(conflicti == -1)
							conflicti = i;
						//caseId = BlockingCase::FVI_NOMATCH_GENERAL_OVERLAP;
						//break;
					}
					if(i < spec.getNdim())
					csize*= sizes[permutation[i]];
#ifdef printd
					cout<<"\ni = "<<i<<"bsize "<<csize<<"\n";
#endif
				}
				if(i == spec.getNdim() + 1)
					i--;
				//if((blockB == 1) && (i < spec.getNdim())) i++;
				blimit = i-1;
				bsize = csize;
#ifdef printd
				cout <<"\nAsize = "<<asize<<" Bsize = "<<bsize<<" alimit = "<<alimit<<" blimit = "<<blimit<<" blockA = "<<blockA<<" blockB = "<<blockB;

#endif
				//double eff = getEfficiency_nooverlap(asize, bsize, sizes[alimit], sizes[permutation[blimit]], blockA, blockB );
				double eff = 0;
				if(!conflict)
					eff       = getEfficiency_nooverlap(asize, bsize, sizes[alimit], sizes[permutation[blimit]], blockA, blockB );
#ifdef printd
				cout <<"\nEff = "<<eff<<"\n";
#endif
				if(eff >= besteff)
				{
					parameters->setNumElementsProcessedPerBlock(asize);
					parameters->setNumElementsProcessedPerBlock1(bsize);
					parameters -> setBlockAIndex(alimit);
					parameters -> setBlockBIndex(blimit);
					parameters->setTileSize(blockA);
					parameters->setTileSize1(blockB);
					besteff = eff;
				}	       
			}//while(rlimit <= maxlimit);
		}
		sharedMemSize1 = 33;
		sharedMemSize2 = 32;
		sharedMemSize = sharedMemSize1 * sharedMemSize2;//csize;//blockB * sizes[permutation[0]];
		paddingSize = 32;
		parameters->setPaddingSize(paddingSize);
		parameters->setWarpEfficiency(besteff);
		double bandwidth = getBW_overlap(besteff);
		parameters -> setBW(bandwidth);
		parameters -> setTime(getTime(bandwidth, spec.getVolume()));
		tbSize = getTBSize(sharedMemSize);
		parameters->setTbSize(tbSize);
		parameters->setSharedMemSize1(sharedMemSize1);
		parameters->setSharedMemSize2(sharedMemSize2);
		return *parameters;


	}



	Parameters& tuneFastestVaryingMatchingCaseWithoutBlocking(TransposeSpec &spec){
		Parameters *parameters = new Parameters();
		// warp efficiency
		float maxWarpEfficiency = 0.0f;
		int tbSize = 32;
		int numMoves = spec.getSizes()[0];
		for(int i= 32; i < 2048 && i <= spec.getSizes()[0]; i+=32){
			double totalThreadsinActiveWarps =  ceil((float)numMoves/(float)32) * 32;
			float warpEfficiency = (float) ((float)spec.getSizes()[0]/(float)(totalThreadsinActiveWarps))*100;
			if(warpEfficiency > maxWarpEfficiency){
				maxWarpEfficiency = warpEfficiency;
				tbSize = (int)totalThreadsinActiveWarps;
			}
			//			std::cout<<"warpefficiency : tbsize = "<<warpEfficiency<< " : "<<tbSize;
		}
		tbSize= 128; // fixing TODO optimal value?
		parameters->setWarpEfficiency(maxWarpEfficiency);
		parameters->setTbSize(tbSize);
		if(spec.getSizes()[0] < 1024 && spec.getSizes()[0] > 256){
			parameters->setTileSize(2);
		}else if(spec.getSizes()[0] <= 256){
			parameters->setTileSize(4);
		}
		double bandwidth = getBW_matchg32();
		parameters -> setBW(bandwidth);
		parameters -> setTime(getTime(bandwidth, spec.getVolume()));
		//cout <<"\nMatching >= 32\n";
		//HashSet<Parameters> returnVal = new HashSet<Parameters>();
		//returnVal.add(parameters);
		parameters->setTileSize(1);
		return *parameters;//returnVal;
	}

	Parameters& tuneFastestVaryingMatchingCaseWithBlocking(TransposeSpec &spec){
		//warp efficiency
		// occupancy
		// Indexing overhead
		int *sizes = spec.getSizes();
		int* permutation = spec.getPermutation();
		int blockA;
		blockA  = (32+sizes[0]-1)/sizes[0];
		int sharedMemSize = 1;
		int tbSize = 32;
		Parameters *parameters = new Parameters();
		int planeSize, numElements, paddingSize;
		//for(i = 1; i < sizes.length; i += 2){ //TODO the blocking happens only at 1 and 2 indices (starting from 0)
		double mintime = 0;
		float occupancy, warpEfficiency;
		int maxPossibleBlocksPerSM;
		float best = 0; int bblock = 1;
		double bf;
		/*do {
		  numElements = blockA*blockA*sizes[0];
		  planeSize = blockA * sizes[0];
		  paddingSize = (32 - (planeSize % 32) + sizes[0])%32;
		//cout <<"here "<<paddingSize<<"\n";
		sharedMemSize = (planeSize + paddingSize)* blockA;
		if(sharedMemSize > sharedMemLimitPerSM/6)

		break;

		maxPossibleBlocksPerSM = sharedMemLimitPerSM/sharedMemSize;

		//if(maxPossibleBlocksPerSM > numThreadBlocksLimitPerSM)
		//	continue;

		// occupancy
		// warp efficiency

		tbSize = blockA * min(32, blockA * sizes[0]); //blockA warps
		if(tbSize <= threadBlocKSizeLimit){
		if (numThreadsLimitPerSM/tbSize > maxPossibleBlocksPerSM) //which ever is minimum, use that : numthreadLimit or sharedMemLimit
		occupancy = ((float)(tbSize * maxPossibleBlocksPerSM) / (float) numThreadsLimitPerSM) * 100;
		else
		occupancy = ((float)(tbSize * (numThreadsLimitPerSM/tbSize)) /(float) numThreadsLimitPerSM) * 100;
		double index =  ceil((float)planeSize/(float)32) ;
		double totalThreadsinActiveWarps =  index * 32;
		const int remainder1 = sizes[1] % blockA;
		const int remainder2 = sizes[permutation[1]] % blockA;

		const int ilimit = remainder1 * sizes[0];
		const int olimit = remainder2 * sizes[0];
		const int plain = blockA * sizes[0];
		double f1, f2, f3, f4, f;
		int minlimit = min(ilimit, olimit);
		f1 =  ((plain/32)  + (double)(plain%32) /32)/ (int)((plain+31)/32);
		f2 =  ((ilimit/32)  + (double)(ilimit%32) /32)/ (int)(max(1,(plain+31)/32));
		f3 =  ((olimit/32)  + (double)(olimit%32) /32)/ (int)(max(1,(plain+31)/32));
		f4 =  ((minlimit/32)  + (double)(minlimit%32) /32)/ (int)(max(1,(plain+31)/32));
		//printf("\tf1=%lf\t", f1 =  ((plain/32)  + (double)(plain%32) /32)/ (int)((plain+31)/32));
		//	printf("\tf2=%lf\t", f2 =  ((ilimit/32)  + (double)(ilimit%32) /32)/ (int)(max(1,(plain+31)/32)));
		//	printf("\tf3=%lf\t", f3 =  ((olimit/32)  + (double)(olimit%32) /32)/ (int)(max(1,(plain+31)/32)));
		//	printf("\tf4=%lf\t", f4 =  ((minlimit/32)  + (double)(minlimit%32) /32)/ (int)(max(1,(plain+31)/32)));
		int asize = sizes[1];
		int bsize = sizes[permutation[1]];
		//	printf("\t%d\t%d\t%d\t%d\t", asize/blockA, asize%blockA, bsize/blockA,bsize%blockA );
		//int amax = min(blockA, 32);
		//int bmax = min(blockB, 32);
		int amax = blockA;
		int bmax = blockA;
		//printf("\tf=%lf\t", f = ((asize/amax) * (bsize/bmax) *f1 + (double)(asize/amax) * (bsize%bmax > 0) *f3+ (double)(asize%amax>0) * (bsize/bmax)*f2 + (double)(asize%amax > 0) * (bsize%bmax > 0) *f4 )/ (int)(((asize+amax-1)/amax) * ((bsize+bmax-1)/bmax)));
		f = ((asize/amax) * (bsize/bmax) *f1 + (double)(asize/amax) * (bsize%bmax > 0) *f3+ (double)(asize%amax>0) * (bsize/bmax)*f2 + (double)(asize%amax > 0) * (bsize%bmax > 0) *f4 )/ (int)(((asize+amax-1)/amax) * ((bsize+bmax-1)/bmax));
		//cout <<"f = "<<f<<" blbock = "<<blockA<<" ";

		warpEfficiency = f;//((float)planeSize / (float)totalThreadsinActiveWarps)*100;
		if(warpEfficiency >= best) {best = warpEfficiency; bblock = blockA; bf = f;}
		}
		blockA++;
		}
		while((blockA < sizes[1]) && (blockA < sizes[permutation[1]]));
		*/	bblock = 8;
		int mul;
		if(sizes[0] <=8) mul = 16;
		else if(sizes[0] <= 16) mul = 8;
		else mul = 4;
		bblock = min(min(mul, sizes[1]), sizes[permutation[1]]);	
		tbSize = bblock * min(32, bblock * sizes[0]); //blockA warps
		//cout<<"\t"<<bf << "\t"<<bblock<<"\t"<<tbSize<<"\t";
		numElements = bblock*bblock*sizes[0];
		planeSize = bblock * sizes[0];
		paddingSize = (32 - (planeSize % 32) + sizes[0])%32;
		//tbSize = bblock * 32; //blockA warps
		if (numThreadsLimitPerSM/tbSize > maxPossibleBlocksPerSM) //which ever is minimum, use that : numthreadLimit or sharedMemLimit
			occupancy = ((float)(tbSize * maxPossibleBlocksPerSM) / (float) numThreadsLimitPerSM) * 100;
		else
			occupancy = ((float)(tbSize * (numThreadsLimitPerSM/tbSize)) /(float) numThreadsLimitPerSM) * 100;
		parameters->setNumElementsProcessedPerBlock(numElements);
		int rem = (sizes[permutation[1]] % bblock) * (sizes[1] % bblock) * sizes[0];
		sharedMemSize = (planeSize + paddingSize)* bblock;
		parameters->setRemElements(rem);
		parameters->setPaddingSize(paddingSize);
		parameters->setOccupancy(occupancy);
		parameters->setWarpEfficiency(best);
		parameters->setTbSize(tbSize);
		parameters->setNumBlocksPerSM(maxPossibleBlocksPerSM);
		parameters->setTileSize(bblock);
		parameters->setSharedMemSize1(sharedMemSize);
		double eff = getEfficiency_matchl32(sizes[0], sizes[1],sizes[permutation[1]], bblock);
		//cout <<" Eff = "<<eff <<"\t";
		double bandwidth = getBW_matchl32(eff, bblock);
		parameters -> setBW(bandwidth);
		parameters -> setTime(getTime(bandwidth, spec.getVolume()));


		return *parameters;
	}

	ParameterTuner() {
		// TODO Auto-generated constructor stub
		sharedMemLimitPerSM = 6144; // (48 * 1024)/8 words per
		numThreadsLimitPerSM =  2048; //per SM
		numThreadBlocksLimitPerSM =  16; //per SM
		threadBlocKSizeLimit = 1024;
		numSMs = 15;
		blockFactor = 6;

	}

	~ParameterTuner() {
		// TODO Auto-generated destructor stub
	}
	};

