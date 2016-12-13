///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioCorrectionProcess.cpp - Self-correction using FM-index walk for PacBio reads
//
#include "PacBioCorrectionProcess.h"
#include "CorrectionThresholds.h"
#include "HashMap.h"
#include <iomanip>
#include "SAIPBSelfCTree.h"
#include "SAIPBHybridCTree.h"
#include "RollingPBSelfCTree.h"

#include "Timer.h"
using namespace std;


PacBioCorrectionProcess::PacBioCorrectionProcess(const PacBioCorrectionParameters params) : m_params(params)
{
}

PacBioCorrectionProcess::~PacBioCorrectionProcess()
{

}


// PacBio Self Correction by Ya and YTH, v20151202.
// 1. Identify highly-accurate seeds within PacBio reads
// 2. For each pair of seeds, perform kmer extension using local kmer frequency collected by FM-index extension
PacBioCorrectionResult PacBioCorrectionProcess::PBSelfCorrection(const SequenceWorkItem& workItem)
{	
	// std::cout << workItem.read.id << "\n";

	PacBioCorrectionResult result;
	
	std::vector<SeedFeature> seedVec, pacbioCorrectedStrs;
	std::string readSeq = workItem.read.seq.toString();
		
	// find seeds using fixed or dynamic kmers depending on 1st round or not
	seedVec = seedingByPacBio_v2(readSeq);
		
	result.totalSeedNum = seedVec.size();

	// push the first seed into pacbioCorrectedStrs, which will be popped later as source seed
	if(seedVec.size() >= 2)
	{
		result.correctedLen += seedVec.at(0).seedStr.length();
		// if(m_params.isSplit)
			pacbioCorrectedStrs.push_back(seedVec.at(0));
		// else
			// pacbioCorrectedStrs.push_back(readSeq.substr(0, seedVec.at(0).seedEndPos+1));
		
		// if(!m_params.isSplit)
			// pacbioCorrectedStrs.back().seedStr.reserve(readSeq.length()*1.5);
	}
	else
	{
		// give up reads with less than 2 seeds
		result.merge = false;
		return result;
	}
	
	realCorrect(readSeq, seedVec, pacbioCorrectedStrs, result);
		
	result.merge = true;
	result.totalReadsLen = readSeq.length();
	for(size_t result_count = 0 ; result_count < pacbioCorrectedStrs.size() ; result_count++)
		result.correctedPacbioStrs.push_back(pacbioCorrectedStrs[result_count].seedStr);
	
	return result;
}


void PacBioHybridCorrectionProcess::seedingByPacBio_v2(const string& readSeq, std::vector<SeedFeature>& seedVec)
{	
	// computing kmer threshold of various kmer size
	// kmer size	kmer threshold*(PB reads coverage/60)+5
	// 17	8*(PB reads coverage/60)+5
	// 27	7*(PB reads coverage/60)+5
	// 37	6*(PB reads coverage/60)+5
	// 47	5*(PB reads coverage/60)+5
	// 57	4*(PB reads coverage/60)+5
	// 67	3*(PB reads coverage/60)+5
	// 77	2*(PB reads coverage/60)+5
	// 87	1*(PB reads coverage/60)+5
	// 97	0*(PB reads coverage/60)+5
	// regression: kmer threshold=(-0.1*(kmer size)+9.7)*(PB reads coverage/60)+5
	std::vector<float> kmerThreshold;
	kmerThreshold.resize(97+1,5);
	for(size_t kmerSize=0 ; kmerSize<=97 ; kmerSize++)
	{
		float kmerThresholdValue=(-0.1*kmerSize+9.7)*((float)m_params.PBcoverage/60);
		kmerThreshold.at(kmerSize)+=kmerThresholdValue;
	}
	
	// search for solid kmers as seeds
	for(size_t pos=0; pos+m_params.PBKmerLength<readSeq.length(); pos++)
	{
		size_t dynamicKmerSize = m_params.PBKmerLength;

		string kmer=readSeq.substr(pos,dynamicKmerSize);
		BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_params.PBindices.pRBWT, reverse(kmer));
		BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_params.PBindices.pBWT, reverseComplement(kmer));
		size_t kmerFreqs= (fwdInterval.isValid()?fwdInterval.size():0) + (rvcInterval.isValid()?rvcInterval.size():0);
		size_t dynamicKmerThreshold=kmerThreshold.at(dynamicKmerSize);
		
		if(kmerFreqs<dynamicKmerThreshold)
			continue;

		size_t seedStartPos=pos;
		size_t maxKmerFreq=kmerFreqs;
		
		// search for longest solid kmer as one seed if possible
		for(pos=pos+dynamicKmerSize ; pos+dynamicKmerSize<readSeq.length() ; pos++)
		{
			char b=readSeq.at(pos);
			char rcb;
			switch(b)
			{
				case 'A': rcb='T'; break;
				case 'T': rcb='A'; break;
				case 'C': rcb='G'; break;
				case 'G': rcb='C'; break;
			}
			BWTAlgorithms::updateInterval(fwdInterval,b,m_params.PBindices.pRBWT);
			BWTAlgorithms::updateInterval(rvcInterval,rcb,m_params.PBindices.pBWT);
			kmerFreqs = (fwdInterval.isValid()?fwdInterval.size():0) + (rvcInterval.isValid()?rvcInterval.size():0);
			
			dynamicKmerSize++;
			// assert(maxKmerFreq > kmerFreqs);
			if(dynamicKmerSize>=kmerThreshold.size()) break;
			
			dynamicKmerThreshold=kmerThreshold.at(dynamicKmerSize);
			
			if(kmerFreqs>=dynamicKmerThreshold)
				maxKmerFreq=kmerFreqs;
			else
			{
				// break and reset dynamicKmerSize and maxKmerFreq
				dynamicKmerSize--;
				dynamicKmerThreshold=kmerThreshold.at(dynamicKmerSize);
				break;
			}
		}

		// //skip repeat seeds
		// if(maxKmerFreq >= m_params.PBcoverage*2)
			// continue;

		// require sufficient seed length in repeats
		if(maxKmerFreq >= m_params.PBcoverage && dynamicKmerSize-m_params.PBKmerLength <= 4)
			// std::cout << readSeq.substr(seedStartPos, i+dynamicKmerSize-2-seedStartPos+1) << "\n";
			continue;

		size_t seedEndPos = pos-1;

		// super repeat seeds with frequency > 2000 are troublesome, often lead to -3 but no good solution so far, mark first
		bool isSuperRepeat = maxKmerFreq >= m_params.PBcoverage?true:false;
		SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), isSuperRepeat, dynamicKmerSize, m_params.PBcoverage/2);
		newSeed.estimateBestKmerSize(m_params.PBindices.pBWT);

		// if(maxKmerFreq > m_params.PBcoverage/4 /*|| HQkmerFreq > m_params.coverage*2*/)
		// {
			// // super repeat is found in PB index, activate MSA instead
			// // if(maxKmerFreq > m_params.PBcoverage*2 && !seedVec.empty())
			// // {
				// // seedVec.back().isNextRepeat = true;
				// // pos = seedEndPos + pow(HQkmerFreq, 0.3)+10;
			// // }
			
			// continue;
		// }

		// skip low-complexity sequencing errors of PacBio
		// bool isShortAndHighFreq = i-seedStartPos <= 2 && maxKmerFreq > 80;
		if(!isLowComplexity(newSeed.seedStr,0.8))
		{
			// debug
			// std::cout << ">>" << seedStartPos << "." << newSeed.seedLength << ":" << maxKmerFreq << "\t" << dynamicKmerSize << "\t" << HQkmerFreq <<  "\n" << newSeed.seedStr << "\n";
			newSeed.isPBSeed = true;
			seedVec.push_back(newSeed);
		}
	}	
}


// Perform FMindex extension between source and target seeds
// Return FMWalkReturnType
int PacBioCorrectionProcess::extendBetweenSeeds(SeedFeature& source, SeedFeature& target, std::string& rawSeq, std::string& mergedseq, 
												size_t extendKmerSize, size_t dis_between_src_target)
{
	const double maxRatio = 1.1;
	const double minRatio = 0.9;
	const int minOffSet = 30;	//PB159615_16774.fa contains large indels > 30bp

	size_t srcKmerSize = std::max(source.endBestKmerSize, extendKmerSize);

	std::string rawSubseq = rawSeq.substr(target.seedStartPos+1-dis_between_src_target-srcKmerSize, target.seedEndPos-source.seedStartPos+1);	
	SAIPBSelfCorrectTree SAITree(m_params.indices.pBWT, m_params.indices.pRBWT, rawSubseq, m_params.FMWKmerThreshold);
		
	// this occurs when one end is repeat requiring larger extendKmerSize while the other requires small extendKmerSize
	// the source str should be the longest one
	// const int srcMaxLength = maxRatio*(dis_between_src_target+minOffSet) + source.seedLength + extendKmerSize;
	// size_t sourceFreq = SAITree.addHashBySingleSeed(source.seedStr, source.endBestKmerSize, extendKmerSize, srcMaxLength, m_params.isFirst);
	
	//size_t srcKmerSize = std::max(source.endBestKmerSize, extendKmerSize);
	std::string srcStr = source.seedStr.substr(source.seedStr.length()-srcKmerSize);
	const int srcMaxLength = maxRatio*(dis_between_src_target+minOffSet) + srcStr.length() + extendKmerSize;
	size_t sourceFreq = SAITree.addHashBySingleSeed(srcStr, source.endBestKmerSize, extendKmerSize, srcMaxLength, m_params.isFirst);
	
	// Collect local kmer frequency from target upto targetMaxLength
	std::string rvcTargetStr = reverseComplement(target.seedStr);
	const int targetMaxLength = maxRatio*(dis_between_src_target+minOffSet) + rvcTargetStr.length() + extendKmerSize;
	size_t expectedLength = dis_between_src_target + rvcTargetStr.length();
	assert(rvcTargetStr.length()>=extendKmerSize);
	size_t targetFreq = SAITree.addHashBySingleSeed(rvcTargetStr, target.startBestKmerSize, extendKmerSize, targetMaxLength, m_params.isFirst, expectedLength);

	// Estimate upper/lower/expected bounds of search depth
	// int srcMinLength = minRatio*(dis_between_src_target-minOffSet) + source.seedLength + extendKmerSize;
	// if(srcMinLength < 0) srcMinLength = 0;
	// expectedLength = source.seedLength + dis_between_src_target + target.seedLength;	
	// int FMWalkReturnType = SAITree.mergeTwoSeedsUsingHash(source.seedStr, target.seedStr, mergedseq, extendKmerSize, m_params.maxLeaves,
													  // srcMinLength, srcMaxLength, expectedLength);
	
	int srcMinLength = minRatio*(dis_between_src_target-minOffSet) + srcStr.length() + extendKmerSize;
	if(srcMinLength < 0) srcMinLength = 0;
	expectedLength = srcStr.length() + dis_between_src_target + target.seedLength;
	int FMWalkReturnType = SAITree.mergeTwoSeedsUsingHash(srcStr, target.seedStr, mergedseq, extendKmerSize, m_params.maxLeaves,
													  srcMinLength, srcMaxLength, expectedLength);
	
	// retain only extended portion
	if(!mergedseq.empty())
		mergedseq = mergedseq.substr(srcStr.length());
	
	// std::cout << source.seedStartPos << "-" << source.seedStartPos+source.seedLength-1 <<  ":" << source.seedLength << ", " 
	// <<	target.seedStartPos << "-" << target.seedStartPos+target.seedLength-1 <<  ":" << target.seedLength 
	// << ", dis: " << dis_between_src_target << ", " << expectedLength << ", " << srcMaxLength << ", "
	// << FMWalkReturnType << ".\n";
	
	// repeat seeds also lead to -1 or -4, don't reduce kmer and give up this seed
	// PB80779_18480.fa 1,596,115-1,598,135 lead to wrong extension due to over-reduction of repeat source seed
	if( !m_params.isFirst && ((FMWalkReturnType==-1 && sourceFreq > (size_t) m_params.seedKmerThreshold*3) || 
		(FMWalkReturnType==-4 && targetFreq > (size_t) m_params.seedKmerThreshold*3)) )
			FMWalkReturnType = -2;

		
	return FMWalkReturnType;
}
// refine seed interval using larger kmer
std::pair<size_t, size_t> PacBioCorrectionProcess::refineRepeatSeed(const std::string readSeq, size_t& seedStartPos, size_t& seedEndPos)
{
	// initially set to max unisnged int value
	size_t newSeedStartPos = (size_t)-1;
	size_t newSeedEndPos = (size_t)-1;
	size_t startKmerFreq=0, endKmerFreq=0;
	
	const int minRepeatFreq = 40, minFreqDiff = 30;
	
	size_t kmerSize = m_params.kmerLength;
	
	std::string kmer = readSeq.substr(seedStartPos, kmerSize);
	int initKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
	int prevKmerFreq = initKmerFreq;

	// first kmer is a repeat
	if(initKmerFreq > minRepeatFreq)
	{
		newSeedStartPos = seedStartPos;
		startKmerFreq = initKmerFreq;
	}
	
	
	// identify breakpoints of large freq difference between two kmers	
	for(size_t i=seedStartPos+1 ; i+kmerSize-1 <= seedEndPos; i++)
	{
		kmer = readSeq.substr(i, kmerSize);
		int currKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);

		// std::cout << i << ": " << kmer << "\t" << currKmerFreq << "\n";

		// error kmers within repeats often lead to large freq diff
		bool isLargeFreqDiff = currKmerFreq - prevKmerFreq > minFreqDiff;
		
		// PB36993_4517.fa, TTATGTAAGGAGTATTTGAT
		// some error kmers are with moderate frequency and no freq diff can be observed
		// pick up the first repeat kmer as starting point
		bool isRepeatKmer = (newSeedStartPos == (size_t)-1) && (currKmerFreq >= (int)minRepeatFreq);
		if(isLargeFreqDiff || isRepeatKmer)
		{
			// capture the highest repeat start pos
			bool isBetterRepeatKmer = (startKmerFreq!=0 && currKmerFreq > (int)startKmerFreq);
			if(newSeedStartPos == (size_t)-1 || isBetterRepeatKmer)
			{
				newSeedStartPos = i;
				startKmerFreq = currKmerFreq;
			}
		}
			
		// repeat end is reached
		// PB36993_4517.fa, AGGCTTGTCTGTAATCGGG
		if(prevKmerFreq - currKmerFreq > minFreqDiff /*|| currKmerFreq < minFreqDiff*/)
		{
			// do not enter unless start pos was found
			// if(newSeedStartPos != (size_t)-1)
			// {
				newSeedEndPos = i + kmerSize -2;
				endKmerFreq = prevKmerFreq;
				break;
			// }
		}
			
		prevKmerFreq = currKmerFreq;
	}
	
	if(newSeedStartPos == (size_t)-1)
	{
		newSeedStartPos = seedStartPos;
		startKmerFreq = initKmerFreq;
	}
		
	if(newSeedEndPos == (size_t)-1)
	{
		newSeedEndPos = seedEndPos;
		endKmerFreq = prevKmerFreq;
	}
	
	// std::cout << newSeedStartPos << "\t" << newSeedEndPos << "\n";

	seedStartPos = newSeedStartPos;
	seedEndPos = newSeedEndPos;
	return std::make_pair(startKmerFreq, endKmerFreq);
}

// return <0: give up and break
// return 0: retry the same target
// return >0: continue to next target
int PacBioCorrectionProcess::FMWalkFailedActions(size_t& extendKmerSize, size_t& numOfTrials, 
								SeedFeature& source, SeedFeature& target, int FMWalkReturnType, int prevFMWalkReturnType)
{
	numOfTrials++;
	// extension failed due to insufficient kmers, reduce large and small kmer sizes
	if(FMWalkReturnType==-1 || FMWalkReturnType==-4)
	{
		// kmers have been enlarged due to repeats, shrink will lead to infinite loop
		if(prevFMWalkReturnType==-3 )
			return -1;
		
		// PB36993_4517.fa, AGGCTTGTCTGTAATCGGG
		if(m_params.isFirst /*&& (source.isRepeat || target.isRepeat)*/)
			return -1;		
		
		source.endBestKmerSize -=2;
		
		target.startBestKmerSize -=2;
		
		extendKmerSize -= 2;

		// std::cout << source.endBestKmerSize << "\t" << target.startBestKmerSize << "\n";

		
		// don't aggressively reduce kmer in the 1st found where most kmers are errors
		if(m_params.isFirst && (source.endBestKmerSize < 15 || target.startBestKmerSize < 15))
			return -1;
			
		if(source.endBestKmerSize < 11 || target.startBestKmerSize < 11 || extendKmerSize < 9)
			return -1;
			
		return 0;
	}
	
	// increase extendKmerSize for reducing repeats
	else if(FMWalkReturnType==-3)
	{
		if(prevFMWalkReturnType==-4 || prevFMWalkReturnType==-1)
			return -1;

		// exponential growth is required in super large repeats. Otherwise speed is too slow
		source.endBestKmerSize += pow(2, numOfTrials+1);
		target.startBestKmerSize += pow(2, numOfTrials+1);
		extendKmerSize += pow(2, numOfTrials+1);
				
		// bug: PB7017_14581_0_14647.fa
		// extendKmerSize is less than seedLength , dunno why
		if(source.seedLength < source.endBestKmerSize || target.seedLength < target.startBestKmerSize ||
			source.seedLength < extendKmerSize || target.seedLength < extendKmerSize )
			return -1;
		
		return 0;
	}
	else if(FMWalkReturnType==-2)
	{
		// probably chimera, need more observations
		// largeKmerSize = m_params.kmerLength;
		// extendKmerSize = largeKmerSize - 2;
		return 1;
	}
	
	return 1;
}

bool PacBioCorrectionProcess::isLowComplexity (std::string seq, float & GCratio, float threshold)
{
	size_t seqLen = seq.length();
	size_t countG =0 ;
	size_t countC =0 ;
	size_t countT =0 ;
	size_t countA =0 ;

	for (size_t i=0; i<seqLen; i++)
	{
		switch(seq[i]){
			case 'A': countA ++ ;break;
			case 'T': countT ++ ;break;
			case 'C': countC ++ ;break;
			case 'G': countG ++ ;break;
			default:  assert(false);
		}
	}

	GCratio = (float)(countG+countC)/seqLen ;

	if (  ((float) countA/seqLen >= threshold ) || ((float) countT/seqLen >= threshold)
			|| ((float) countC/seqLen >= threshold ) || ((float) countG/seqLen >= threshold) )
		return true;

	return false;

}

//
//
//
PacBioCorrectionPostProcess::PacBioCorrectionPostProcess(std::ostream* pCorrectedWriter,
std::ostream* pDiscardWriter,
const PacBioCorrectionParameters params) :
m_pCorrectedWriter(pCorrectedWriter),
m_pDiscardWriter(pDiscardWriter),
m_params(params),
m_totalReadsLen(0),
m_correctedLen(0),
m_totalSeedNum(0),
m_totalWalkNum(0),
m_correctedNum(0),
m_highErrorNum(0),
m_exceedDepthNum(0),
m_exceedLeaveNum(0),
m_seedDis(0)
{
}

//
PacBioCorrectionPostProcess::~PacBioCorrectionPostProcess()
{
	if(m_totalWalkNum>0 && m_totalReadsLen>0)
	{
		std::cout << std::endl;
		std::cout << "totalReadsLen: " << m_totalReadsLen << ", ";
		std::cout << "correctedLen: " << m_correctedLen << ", ratio: " << (float)(m_correctedLen)/m_totalReadsLen << "." << std::endl;
		std::cout << "totalSeedNum: " << m_totalSeedNum << "." << std::endl;
		std::cout << "totalWalkNum: " << m_totalWalkNum << ", ";
		std::cout << "correctedNum: " << m_correctedNum << ", ratio: " << (float)(m_correctedNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "highErrorNum: " << m_highErrorNum << ", ratio: " << (float)(m_highErrorNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "exceedDepthNum: " << m_exceedDepthNum << ", ratio: " << (float)(m_exceedDepthNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "exceedLeaveNum: " << m_exceedLeaveNum << ", ratio: " << (float)(m_exceedLeaveNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "disBetweenSeeds: " << m_seedDis/m_totalWalkNum << std::endl << std::endl;
	}
}


// Writting results for kmerize and validate
void PacBioCorrectionPostProcess::process(const SequenceWorkItem& item, const PacBioCorrectionResult& result)
{	
	if (result.merge)
	{
		m_totalReadsLen += result.totalReadsLen;
		m_correctedLen += result.correctedLen;
		m_totalSeedNum += result.totalSeedNum;
		m_totalWalkNum += result.totalWalkNum;
		m_correctedNum += result.correctedNum;
		m_highErrorNum += result.highErrorNum;
		m_exceedDepthNum += result.exceedDepthNum;
		m_exceedLeaveNum += result.exceedLeaveNum;
		m_seedDis += result.seedDis;

		//cout << result.correctSequence.toString();
		/*SeqItem mergeRecord;
		stringstream ss;
		ss << item.read.id << "_before_len:" << result.correctSequence.toString().length();
		mergeRecord.id = ss.str();
		mergeRecord.seq = result.correctSequence;
		mergeRecord.write(*m_pCorrectedWriter);*/
		
		for(size_t i = 0 ; i < result.correctedPacbioStrs.size() ; i++)
		{
			SeqItem mergeRecord2;
			std::stringstream ss2;
			ss2 << item.read.id << "_" << i << "_" << result.correctedPacbioStrs[i].toString().length();
			mergeRecord2.id = ss2.str();
			mergeRecord2.seq = result.correctedPacbioStrs[i];
			mergeRecord2.write(*m_pCorrectedWriter);
		}
	}
	else
	{
		// write into discard.fa
		SeqItem mergeRecord2;
		mergeRecord2.id = item.read.id;
		mergeRecord2.seq = item.read.seq;
		mergeRecord2.write(*m_pDiscardWriter);
	}
}
