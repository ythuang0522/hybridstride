//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioCorrectionProcess - Self-correction using FM-index walk for PacBio reads
//

#ifndef SeedFeature_H
#define SeedFeature_H

// #include "HashMap.h"
// #include "Util.h"
// #include "SequenceProcessFramework.h"
// #include "SequenceWorkItem.h"
// #include "Metrics.h"
#include "BWTIndexSet.h"
#include "SampledSuffixArray.h"
#include "BWTAlgorithms.h"
#include "KmerDistribution.h"

struct SeedFeature
{
	public:
		SeedFeature(size_t startPos, std::string str, bool repeat, size_t kmerSize, size_t repeatCutoff);

		SeedFeature(){};
		
		// append current seed string with extendedStr
		void append(std::string extendedStr);
		void setBestKmerSize(size_t kmerSize);		
		void estimateBestKmerSize(const BWT* pBWT);
		bool inline isSmall(){return seedLength<=17?true:false ;}
		
		size_t seedStartPos;
		size_t seedEndPos;
		size_t seedLength;
		std::string seedStr;
		bool isRepeat;
		bool isPBSeed;
		bool isNextRepeat;
		
		// estimated by calling estimateBestKmerSize
		size_t startBestKmerSize;
		size_t endBestKmerSize;
		size_t startKmerFreq;
		size_t endKmerFreq;
		
	private:
		size_t freqUpperBound;
		size_t freqLowerBound;
		size_t minKmerSize;
		size_t stepSize;
		//estimate kmer size
		void increaseStartKmerSize(const BWT* pBWT);
		void decreaseStartKmerSize(const BWT* pBWT);

		//estimate kmer size
		void increaseEndKmerSize(const BWT* pBWT);
		void decreaseEndKmerSize(const BWT* pBWT);
};

#endif
