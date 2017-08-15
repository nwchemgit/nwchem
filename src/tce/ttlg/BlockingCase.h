/*
 * BlockingCase.h
 *
 *  Created on: Jan 30, 2017
 *      Author: arjun
 */

#ifndef BLOCKINGCASE_H_
#define BLOCKINGCASE_H_

class BlockingCase {
public:
	//enum BlockingCase {
static enum blockingcase {
		FVI_MATCH_AND_LESST32 = 0,
		FVI_MATCH_AND_GREATERT32 = 1,
		FVI_NOMATCH_ALL_DIFFERENT = 2,
		FVI_NOMATCH_COMMON_SECOND_INDEX = 3,
		FVI_NOMATCH_BOTH_SECOND_INDICES_REPEATED_BUT_DIFFERENT = 4,
		FVI_NOMATCH_ONE_SECOND_INDEX_REPEATED = 5,
		FVI_NOMATCH_AND_GREATERT32 = 6,
		FVI_NOMATCH_GENERAL = 7,
		FVI_NOMATCH_GENERAL_OVERLAP = 8
}myblockingcase;

		 int mode;
		 BlockingCase(int mode){
			this->mode = mode;
		}

		 int getMode() {
			return mode;
		}
		 bool operator==(blockingcase A){
			 return mode == A;
		 }

		/* static BlockingCase valueOf(int type) {
			for (BlockingCase vMode : values()) {
				if (vMode.mode == type)
					return vMode;
			}

		}*/

		/* static BlockingCase OfName(String name) {
			for (BlockingCase vMode : BlockingCase.class.getEnumConstants()) {
				if (vMode.name().equals(name))
					return vMode;
			}*/



	BlockingCase();
	virtual ~BlockingCase();
};

#endif /* BLOCKINGCASE_H_ */
