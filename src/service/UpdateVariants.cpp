//
// Created by song on 8/22/18.
//

#include "UpdateVariants.h"

bool overlapIgnoreStrand(Transcript * t1, Transcript * t2){
    if( t1->getPStart() < t2->getPStart() && t2->getPStart()< t1->getPEnd() ){
        return true;
    } else if ( t1->getPStart() < t2->getPEnd() && t2->getPEnd()< t1->getPEnd() ){
        return true;
    } else if ( t2->getPStart() < t1->getPStart() && t1->getPStart()< t2->getPEnd() ){
        return true;
    } else if ( t2->getPStart() < t1->getPEnd() && t1->getPEnd()< t2->getPEnd() ){
        return true;
    }
    return false;
}

void updateVariants( const std::string & referenceGffFilePath, const std::string & referenceGenomeFilePath,
                     const std::string & sdiFile, const std::string & vcfFix, const size_t & minIntron,
                     std::map<std::string, std::string>& parameters, const std::string & outPutFilePath  ){

    //size_t maxLengthForStructureAlignment = 10000;

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    std::map<std::string, Fasta> referenceGenome;
    readFastaFile(referenceGenomeFilePath, referenceGenome);

    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    std::string regex = get_parameters("cdsParentRegex", parameters);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);
    CheckAndUpdateTranscriptsEnds( referenceTranscriptHashSet, referenceGenome, nucleotideCodeSubstitutionMatrix, minIntron); //the referenceTranscriptHashSet has been sorted

    std::set<std::string> toRemoveChromosomes;
// clean data begin
    for( std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
         it!=referenceTranscriptHashSet.end(); ++it){
        if( referenceGenome.find(it->first)==referenceGenome.end() ){
            toRemoveChromosomes.insert(it->first);
        }
    }
    for( std::map<std::string, Fasta >::iterator it=referenceGenome.begin();
         it!=referenceGenome.end(); ++it){
        if( referenceTranscriptHashSet.find(it->first)==referenceTranscriptHashSet.end() ){
            toRemoveChromosomes.insert(it->first);
        }
    }
    for( std::set<std::string>::iterator it=toRemoveChromosomes.begin();
         it!=toRemoveChromosomes.end(); ++it){
        if( referenceGenome.find(*it) != referenceGenome.end() ){
            referenceGenome.erase(*it);
        }
        if( referenceTranscriptHashSet.find(*it) != referenceTranscriptHashSet.end() ){
            referenceTranscriptHashSet.erase(*it);
        }
    } // clean data end

    // remove ORF shift transcripts from reference accession
    std::cout << "remove ORF shift allele begin" << std::endl;
    // remove ORF shift allele from reference dataset begin
    for ( std::map<std::string, std::vector<Transcript> >::iterator it1=referenceTranscriptHashSet.begin();
          it1!=referenceTranscriptHashSet.end(); ++it1){
        for (int index=0; index < (it1->second.size()); ++index){
            TranscriptUpdateCdsInformation((it1->second)[index], referenceGenome);
            checkOrfState( (it1->second)[index], referenceGenome, nucleotideCodeSubstitutionMatrix, minIntron);
        }
        std::vector<int> transcriptToRemove;
        for (int index=0; index < (it1->second.size()); ++index){
            if ((it1->second)[index].getIfOrfShift()) {
                transcriptToRemove.push_back(index);
            }
        }// loop from end to begin
        for(int index=transcriptToRemove.size()-1; index>=0; --index){
            referenceTranscriptHashSet[it1->first].erase(referenceTranscriptHashSet[it1->first].begin()+transcriptToRemove[index]);
        }
    } // remove ORF shift allele from reference data set end
    std::cout << "remove ORF shift allele done" << std::endl;

    std::map<std::string, std::vector<Variant> > variantsMaps;
    readSdiFile (sdiFile, variantsMaps, vcfFix, referenceGenome);
    std::map<std::string, Fasta> targetGenome;
    getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome, parameters);

    std::map<std::string, Transcript > targetTranscriptsHashMap;
    annotationLiftOver(referenceTranscriptHashSet, targetTranscriptsHashMap, variantsMaps, targetGenome, referenceGenome, parameters, minIntron);
    std::cout << "annotationLiftOver done" << std::endl;
    std::map<std::string, std::vector<Variant> > newVariantsMaps;
    for(std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
        it!=referenceTranscriptHashSet.end(); ++it) {
        if( variantsMaps.find(it->first) != variantsMaps.end() ){
//            std::cout << "line 86" << std::endl;
            FirstLastList sdiRecords;
            for ( int it2 =0; it2<variantsMaps[it->first].size(); ++it2 ){
                Data *data = new Data(variantsMaps[it->first][it2]);
                sdiRecords.insertLast(data);
            }
  //          std::cout << "line 92" << std::endl;

            int lastRefEnd = 0; // since the variants have been update, so the getBasementFromChanged etc. functions do not work well anymore
            int lastIt2=0;
            for (int it2 =0; it2<it->second.size(); ++it2) {
                Transcript *t2 = &targetTranscriptsHashMap[it->second[it2].getName()];
//                std::cout << "line 95" << std::endl;
                if ( (*t2).getIfOrfShift() && ( it2==0|| ( !overlapIgnoreStrand(&targetTranscriptsHashMap[it->second[it2-1].getName()], t2)) && (!overlapIgnoreStrand(&it->second[it2], &it->second[it2-1])) ) ) {
//                    std::cout << "line 99" << std::endl;
                    std::string refGenomeSequence = it->second[it2].getGeneomeSequence();
                    int startTarget = t2->getPStart() - refGenomeSequence.length();
                    if (startTarget < 1) {
                        startTarget = 1;
                    }
//                    std::cout << "line 105" << std::endl;
                    if (it2>0 && startTarget <= targetTranscriptsHashMap[it->second[it2-1].getName()].getPEnd()) { //so this does not affect the structure of last transcript
                        startTarget = targetTranscriptsHashMap[it->second[it2-1].getName()].getPEnd() + 1;
                    }
//                    std::cout << "line 109 " << startTarget << std::endl;
                    int startReference = getBasementFromChanged( it->first, startTarget, variantsMaps );
//                    std::cout << "line 111" << std::endl;
                    if( it2>0 && startReference <= it->second[it2-1].getPEnd() ){
                        startReference = it->second[it2-1].getPEnd()+1;
                    }
                    if( it2>0 && startReference <= lastRefEnd ){
                        startReference = it->second[lastIt2].getPEnd()+1;
                        startTarget = targetTranscriptsHashMap[it->second[lastIt2].getName()].getPEnd() + 1;
                    }else{
                    //std::vector<Variant> overlappedVariants;
//                        std::cout << "line 117" << std::endl;
                        Data *thisData = sdiRecords.getFirst();
                        while ( thisData!= NULL) {
                            while (   (thisData->getMapSingleRecord().getPosition() - 2) <= startReference &&
                                   ( startReference<= (thisData->getMapSingleRecord().getPosition() + thisData->getMapSingleRecord().getReference().length() +
                                     1) )   ) {
                                ++startReference;
                            }
                            thisData=thisData->getNext();
                        }
                        startTarget = getChangedFromBasement( it->first, startReference, variantsMaps );
                    }
//                    std::cout << "line 122" << std::endl;
                    if (startTarget < t2->getPStart() && startReference < it->second[it2].getPStart()) { // else this transcript could not be realigned
                        int endTarget = t2->getPEnd() + refGenomeSequence.length();
                        if (endTarget > targetGenome[it->first].getSequence().length()) {
                            endTarget = targetGenome[it->first].getSequence().length();
                        }
                        int endReference = getBasementFromChanged(it->first, endTarget, variantsMaps);
                        if (endReference > referenceGenome[it->first].getSequence().length()) {
                            endReference = referenceGenome[it->first].getSequence().length();
                        }

                        Data * thisData = sdiRecords.getLast();
//                        std::cout << "line 130" << std::endl;
                        while (thisData != NULL) {
                            while (  (thisData->getMapSingleRecord().getPosition() - 2) <= endReference &&
                                   (endReference <= (thisData->getMapSingleRecord().getPosition() +
                                                    thisData->getMapSingleRecord().getReference().length() +
                                                    1)) ) {
                                --endReference;
                            }
                            thisData = thisData->getPrev();
                        }
                        endTarget = getChangedFromBasement(it->first, endReference, variantsMaps);
//                        std::cout << "line 141" << std::endl;
                        if(endReference > it->second[it2].getPEnd() && endTarget > t2->getPEnd()){
                            std::string dna_b = getSubsequence(targetGenome, it->first, startTarget, endTarget,
                                                               it->second[it2].getStrand());
                            int startCodonPosition = 1;
                            int stopCodonPosition = refGenomeSequence.length() - 2;

                            std::map<int, int> importantPositions;
                            for (std::vector<GenomeBasicFeature>::iterator it4 = it->second[it2].getCdsVector().begin();
                                 it4 != it->second[it2].getCdsVector().end(); ++it4) {
                                importantPositions[(*it4).getStart()] = 0;
                                importantPositions[(*it4).getEnd()] = 0;
                            }
                            std::string databaseAlignmentResult = "";
                            std::string queryAlignmentResult = "";
                            std::vector<SpliceSitePosition> spliceSitePositions;
                            int targetPosition = 0;
                            int referencePosition = 0;
//                            std::cout << "line 159" << std::endl;
                            if (it->second[it2].getStrand() == POSITIVE) {
                                if (it->second[it2].getCdsVector().size() > 1) {
                                    for (size_t i = 1; i < it->second[it2].getCdsVector().size(); ++i) {
                                        SpliceSitePosition spliceSitePosition(
                                                it->second[it2].getCdsVector()[i - 1].getEnd() - it->second[it2].getPStart() + 2,
                                                it->second[it2].getCdsVector()[i].getStart() - it->second[it2].getPStart());
                                        spliceSitePositions.push_back(spliceSitePosition);
                                    }
                                }
                                AlignTranscript nw(refGenomeSequence, dna_b, startCodonPosition,
                                                   stopCodonPosition, spliceSitePositions, parameters,
                                                   nucleotideCodeSubstitutionMatrix);
                                databaseAlignmentResult = nw.getAlignment_d();
                                queryAlignmentResult = nw.getAlignment_q();
//                                                        if (t2->getCdsVector().size() > 1) {
//                                                            std::cout << t2->getName() << " " << "POSITIVE" << std::endl;
//                                                            nw.print_results();
//                                                        }
                                referencePosition = it->second[it2].getPStart() - 1;
                                for (size_t tp = 0; tp < queryAlignmentResult.length(); ++tp) {
                                    if (queryAlignmentResult[tp] != '-') {
                                        ++targetPosition;
                                    }
                                    if (databaseAlignmentResult[tp] != '-') {
                                        ++referencePosition;
                                        if (importantPositions.find(referencePosition) != importantPositions.end()) {
                                            importantPositions[referencePosition] = startTarget + targetPosition - 1;
                                        }
                                    }
                                }
                            } else {
                                if (it->second[it2].getCdsVector().size() > 1) {
                                    for (size_t i = it->second[it2].getCdsVector().size() - 1; i > 0; --i) {
                                        SpliceSitePosition spliceSitePosition(
                                                it->second[it2].getPEnd() - it->second[it2].getCdsVector()[i].getStart() +
                                                2,
                                                it->second[it2].getPEnd() - it->second[it2].getCdsVector()[i - 1].getEnd());
                                        spliceSitePositions.push_back(spliceSitePosition);
                                    }
                                }
                                AlignTranscript nw(refGenomeSequence, dna_b, startCodonPosition,
                                                   stopCodonPosition, spliceSitePositions, parameters,
                                                   nucleotideCodeSubstitutionMatrix);

                                databaseAlignmentResult = nw.getAlignment_d();
                                queryAlignmentResult = nw.getAlignment_q();
//                                                        if (t2->getCdsVector().size() > 1) {
//                                                            std::cout << t2->getName() << " " << "NEGATIVE" << std::endl;
//                                                            nw.print_results();
//                                                        }
                                referencePosition = it->second[it2].getPEnd() + 1;

                                for (size_t tp = 0; tp < nw.getAlignment_d().length(); ++tp) {
                                    if (queryAlignmentResult[tp] != '-') {
                                        ++targetPosition;
                                    }
                                    if (databaseAlignmentResult[tp] != '-') {
                                        --referencePosition;
                                        if (importantPositions.find(referencePosition) != importantPositions.end()) {
                                            importantPositions[referencePosition] = endTarget - targetPosition + 1;
                                        }
                                    }
                                }
                            }
//                            std::cout << "line 234" << std::endl;
                            Transcript targetTranscript(t2->getName(), t2->getChromeSomeName(), it->second[it2].getStrand());
                            for (std::vector<GenomeBasicFeature>::iterator it4 = it->second[it2].getCdsVector().begin();
                                 it4 != it->second[it2].getCdsVector().end(); ++it4) {
                                GenomeBasicFeature cds(importantPositions[it4->getStart()],
                                                       importantPositions[it4->getEnd()]);
                                targetTranscript.addCds(cds);
                            }
                            TranscriptUpdateCdsInformation(targetTranscript, targetGenome);
                            checkOrfState(targetTranscript, targetGenome, nucleotideCodeSubstitutionMatrix, minIntron);
//                            std::cout << "line 244" << std::endl;
                            if ( (!targetTranscript.getIfOrfShift())) {
                                lastRefEnd = endReference;
                                lastIt2 = it2;
                                targetTranscriptsHashMap[targetTranscript.getName()]=targetTranscript;
                                std::string refUpStream = getSubsequence(referenceGenome, it->first, startReference,
                                                                         it->second[it2].getPStart() - 1);
                                std::string targetUpStream = getSubsequence(targetGenome, it->first, startTarget,
                                                                            targetTranscript.getPStart() - 1);
                                std::string _alignment_refUpStream, _alignment_targetUpStream;
                                int langerSeqLength = refUpStream.length();
                                if (langerSeqLength < targetUpStream.length()) {
                                    langerSeqLength = targetUpStream.length();
                                }
//                                std::cout << "line 256" << std::endl;
                                alignSlidingWindow(refUpStream, targetUpStream, _alignment_refUpStream,
                                                   _alignment_targetUpStream,
                                                   langerSeqLength, parameters, nucleotideCodeSubstitutionMatrix);

                                std::string refDownStream = getSubsequence(referenceGenome, it->first,
                                                                           it->second[it2].getPEnd() + 1, endReference);
                                std::string targetDownStream = getSubsequence(targetGenome, it->first,
                                                                              targetTranscript.getPEnd() + 1, endTarget);
                                std::string _alignment_refDownStream, _alignment_targetDownStream;
                                langerSeqLength = refDownStream.length();
                                if (langerSeqLength < targetDownStream.length()) {
                                    langerSeqLength = targetDownStream.length();
                                }
//                                std::cout << "line 270" << std::endl;
                                alignSlidingWindow(refDownStream, targetDownStream, _alignment_refDownStream,
                                                   _alignment_targetDownStream, langerSeqLength, parameters,
                                                   nucleotideCodeSubstitutionMatrix);

                                if (t2->getStrand() == NEGATIVE) {
                                    databaseAlignmentResult = getReverseComplementary(databaseAlignmentResult);
                                    queryAlignmentResult = getReverseComplementary(queryAlignmentResult);
                                }
                                //trim queryAlignmentResult begin
                                int hk = 0;
                                for (int hi = 0; hi < queryAlignmentResult.length(); ++hi) {
                                    if (queryAlignmentResult[hi] != '-') {
                                        ++hk;
                                        if (hk <= targetUpStream.length()) {
                                            queryAlignmentResult[hi] = '-';
                                        } else if (hk > (targetUpStream.length() +
                                                         targetTranscript.getGeneomeSequence().length())) {
                                            queryAlignmentResult[hi] = '-';
                                        }
                                    }
                                }
//                                std::cout << "line 292" << std::endl;
                                //trim queryAlignmentResult end
                                databaseAlignmentResult =
                                        _alignment_refUpStream + databaseAlignmentResult + _alignment_refDownStream;
                                queryAlignmentResult =
                                        _alignment_targetUpStream + queryAlignmentResult + _alignment_targetDownStream;

                                //for debugging begin
                                if( true ){
                                    std::string ad = databaseAlignmentResult;
                                    ad.erase(std::remove(ad.begin(), ad.end(), '-'), ad.end());
                                    std::string aq = queryAlignmentResult;
                                    aq.erase(std::remove(aq.begin(), aq.end(), '-'), aq.end());
                                    std::string d = getSubsequence(referenceGenome, t2->getChromeSomeName(), startReference, endReference);
                                    std::string q = getSubsequence(targetGenome, t2->getChromeSomeName(), startTarget, endTarget);
                                    if( ad.compare(d) != 0 || aq.compare(q)!=0  ){
                                        std::cerr << ad << std::endl << d << std::endl << aq << std::endl << q << std::endl;
                                    }else{
                                        std::cerr << t2->getName() << " " << t2->getStrand() <<  " good " << startReference << " " << endReference << std::endl;
                                    }
                                }
                                //for debugging end

                                //begin re-orginize the chain data structure
                                Data *theOneUpstream = sdiRecords.getFirst(); // todo think more about the extreme case about the variants chain, length, the position of the fisrt one and the last one
                                Data *theOneDownstream = sdiRecords.getLast();
                                Data *theRealLastOne = sdiRecords.getLast();
//                                std::cout << "line 319" << std::endl;
                                while (theOneUpstream->getNext() != NULL &&
                                       theOneUpstream->getMapSingleRecord().getPosition() < startReference) {
                                    theOneUpstream = theOneUpstream->getNext();
                                }
                                //theOneUpstream->getNext() == NULL could not happen here
//                                    if( theOneUpstream->getPrev() == NULL ){
//                                        continue;
//                                    }
                                theOneUpstream = theOneUpstream->getPrev();
                                // if the first is falling into the target region
                                // then theOneUpstream == NULL

                                while (theOneDownstream->getPrev() != NULL &&
                                       theOneDownstream->getMapSingleRecord().getPosition() > endReference) {
                                    theOneDownstream = theOneDownstream->getPrev();
                                }
                                //theOneDownstream->getPrev() == NULL should not happen here
                                if( theOneDownstream->getNext() == NULL ){
                                    continue;
                                }
                                theOneDownstream = theOneDownstream->getNext();
                                // if the last one is falling into the target region
                                // then theOneDownstream == NULL
//                                std::cout << "line 343" << std::endl;
                                if( theOneUpstream!=NULL ){
                                    std::cout << "theOneUpstream " << theOneUpstream->getMapSingleRecord().getPosition() << " " << theOneUpstream->getMapSingleRecord().getChanginglength() << " "
                                          << theOneUpstream->getMapSingleRecord().getReference() << " " << theOneUpstream->getMapSingleRecord().getAlternative() << std::endl;
                                }
                                if( theOneDownstream!=NULL ) {
                                    std::cout << "theOneDownstream "
                                              << theOneDownstream->getMapSingleRecord().getPosition() << " "
                                              << theOneDownstream->getMapSingleRecord().getChanginglength() << " "
                                              << theOneDownstream->getMapSingleRecord().getReference() << " "
                                              << theOneDownstream->getMapSingleRecord().getAlternative() << std::endl;
                                }
                                //clean RAM begin
                                Data *thisone;
                                if( theOneUpstream != NULL){
                                    thisone = theOneUpstream->getNext();
                                }else{
                                    thisone = sdiRecords.getFirst();
                                }
                                while (thisone != theOneDownstream) {
                                    Data *tempData = thisone;
                                    thisone = (thisone->getNext());
                                    std::cout << "delete " << tempData->getMapSingleRecord().getPosition() << " " << tempData->getMapSingleRecord().getChanginglength() << " "
                                              << tempData->getMapSingleRecord().getReference() << " " << tempData->getMapSingleRecord().getAlternative() << std::endl;
                                    delete (tempData);
                                }
                                //clean RAM end
//                                std::cout << "line 365" << std::endl;
                                sdiRecords.setLast(theOneUpstream);
//                                std::cout << "line 367" << std::endl;
                                //end re-orginize the chain data structure
//                                std::cout << "line 368" << std::endl;
                                //then replace all the variants in the vector with all the new variants begin
                                int refLetterNumber = 0;
                                for (int ai = 0; ai < databaseAlignmentResult.length(); ai++) {
                                    if (databaseAlignmentResult[ai] != '-') {
                                        ++refLetterNumber;
                                    }
                                    if (databaseAlignmentResult[ai] != queryAlignmentResult[ai]) {
                                        if (queryAlignmentResult[ai] == '-') {
                                            int position = startReference + refLetterNumber - 1;
                                            std::string ori(1, databaseAlignmentResult[ai]);
                                            std::string result = "-";
                                            Variant mapSingleRecord = Variant(it->first, position, ori, result);
                                            Data *data = new Data(mapSingleRecord);
                                            sdiRecords.insertLast(data);
                                            std::cout << "insert " << data->getMapSingleRecord().getPosition() << " " << data->getMapSingleRecord().getChanginglength() << " "
                                                      << data->getMapSingleRecord().getReference() << " " << data->getMapSingleRecord().getAlternative() << std::endl;
                                        } else if (databaseAlignmentResult[ai] == '-') {
                                            if (sdiRecords.getLast() == NULL) {
                                                int position = startReference + refLetterNumber;
                                                std::string ori = "-";
                                                std::string result(1, queryAlignmentResult[ai]);
                                                Variant mapSingleRecord = Variant(it->first, position, ori, result);
                                                Data *data = new Data(mapSingleRecord);
                                                sdiRecords.insertLast(data);
                                                std::cout << "insert " << data->getMapSingleRecord().getPosition() << " " << data->getMapSingleRecord().getChanginglength() << " "
                                                          << data->getMapSingleRecord().getReference() << " " << data->getMapSingleRecord().getAlternative() << std::endl;
                                            } else {
                                                if (NULL != sdiRecords.getLast() &&
                                                    sdiRecords.getLast()->getMapSingleRecord().getPosition() ==
                                                    (startReference + refLetterNumber)
                                                    && sdiRecords.getLast()->getMapSingleRecord().getChanginglength() > 0 &&
                                                    sdiRecords.getLast()->getMapSingleRecord().getReference().compare(
                                                            "-") == 0) {

                                                    int position = startReference + refLetterNumber;
                                                    std::string ori = "-";
                                                    std::string result =
                                                            sdiRecords.getLast()->getMapSingleRecord().getAlternative()
                                                            + queryAlignmentResult[ai];
                                                    Variant mapSingleRecord = Variant(it->first, position, ori, result);
                                                    Data *data = new Data(mapSingleRecord);
                                                    sdiRecords.deleteLast();
                                                    sdiRecords.insertLast(data);
                                                    std::cout << "insert " << data->getMapSingleRecord().getPosition() << " " << data->getMapSingleRecord().getChanginglength() << " "
                                                              << data->getMapSingleRecord().getReference() << " " << data->getMapSingleRecord().getAlternative() << std::endl;
                                                } else {
                                                    int position = startReference + refLetterNumber;
                                                    std::string ori = "-";
                                                    std::string result(1, queryAlignmentResult[ai]);
                                                    Variant mapSingleRecord = Variant(it->first, position, ori, result);
                                                    Data *data = new Data(mapSingleRecord);
                                                    sdiRecords.insertLast(data);
                                                    std::cout << "insert " << data->getMapSingleRecord().getPosition() << " " << data->getMapSingleRecord().getChanginglength() << " "
                                                              << data->getMapSingleRecord().getReference() << " " << data->getMapSingleRecord().getAlternative() << std::endl;
                                                }
                                            }
                                        } else {
                                            int position = startReference + refLetterNumber - 1;
                                            std::string ori(1, databaseAlignmentResult[ai]);
                                            std::string result(1, queryAlignmentResult[ai]);
                                            Variant mapSingleRecord = Variant(it->first, position, ori, result);
                                            Data *data = new Data(mapSingleRecord);
                                            sdiRecords.insertLast(data);
                                            std::cout << "insert " << data->getMapSingleRecord().getPosition() << " " << data->getMapSingleRecord().getChanginglength() << " "
                                                      << data->getMapSingleRecord().getReference() << " " << data->getMapSingleRecord().getAlternative() << std::endl;
                                        }
                                    }
                                }
//                                std::cout << "line 437" << std::endl;
                                //re-join the chain data structure
                                if (theOneDownstream != NULL) {
                                    theOneDownstream->setPrev(sdiRecords.getLast());
                                    sdiRecords.getLast()->setNext(theOneDownstream);
                                    sdiRecords.setLast(theRealLastOne);
                                }
                            }
                        }
                    }
                }
            }
//            std::cout << "line 449" << std::endl;
            for( int runingCound = 0; runingCound<2; ++runingCound){
//                std::cout << "round " << runingCound << std::endl;
                if( (sdiRecords.getFirst() != NULL)
                    && (sdiRecords.getFirst()->getNext() != NULL ) ){
                    Data* prevOne = sdiRecords.getFirst();
                    Data* currOne = (sdiRecords.getFirst()->getNext());
    //                    std::cout << "2107" << std::endl;
                    while( (currOne != NULL) && (NULL != currOne->getNext())){
    //                        std::cout << "2109" << std::endl;
                        if( sdiRecords.getFirst() == currOne ){
                            prevOne = currOne;
                            currOne = prevOne->getNext();
                        }
                        //std::cout << "793" << std::endl;
                        if (currOne->getMapSingleRecord().getChanginglength()<0 && prevOne->getMapSingleRecord().getChanginglength()<0 &&
                            (prevOne->getMapSingleRecord().getPosition()+abs(prevOne->getMapSingleRecord().getChanginglength()))
                            ==currOne->getMapSingleRecord().getPosition()
                            && prevOne->getMapSingleRecord().getAlternative().compare("-")==0 &&
                            currOne->getMapSingleRecord().getAlternative().compare("-")==0) { // merge to deletions
                            //std::cout << "797 delete prev" << std::endl;
                            int position = prevOne->getMapSingleRecord().getPosition();
                            std::string ori = prevOne->getMapSingleRecord().getReference() +
                                              currOne->getMapSingleRecord().getReference();
                            std::string result = "-";
                            Variant mapSingleRecord2(it->first, position, ori, result);
                            //delete prev one begin
                            if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                                sdiRecords.deleteFirst();
                            } else {
                                prevOne->getPrev()->setNext(currOne);
                                currOne->setPrev(prevOne->getPrev());
                                delete (prevOne);
                            } //delete prev one end
                            currOne->setMapSingleRecord(mapSingleRecord2);
                            prevOne = currOne->getPrev();
                            if( prevOne == NULL ){
                                prevOne = currOne;
                                currOne = prevOne->getNext();
                            }
                            //std::cout << "817 delete prev" << std::endl;
                        } else if (currOne->getMapSingleRecord().getChanginglength()==0 && 0 ==
                                   currOne->getMapSingleRecord().getReference().compare( currOne->getMapSingleRecord().getAlternative())){ // nonsense records
                            //delete current one
                            //std::cout << "820 delete prev" << std::endl;
                            prevOne->setNext( currOne->getNext() );
                            currOne->getNext()->setPrev(prevOne);
                            delete(currOne);
                            currOne = prevOne->getNext();
                            //std::cout << "825 delete prev" << std::endl;
                        } else if ( currOne->getMapSingleRecord().getChanginglength()<0 &&
                                    prevOne->getMapSingleRecord().getChanginglength()>0
                                    && currOne->getMapSingleRecord().getReference().compare(prevOne->getMapSingleRecord().getAlternative())==0 &&
                                    currOne->getMapSingleRecord().getPosition() == prevOne->getMapSingleRecord().getPosition()){ //delete one insertion and next reverse sence deletion
                            //delete current one and prev
                            //std::cout << "830 delete current one and prev" << std::endl;
                            if( ((currOne->getPrev())) == (sdiRecords.getFirst()) ){
                                sdiRecords.deleteFirst();
                                sdiRecords.deleteFirst();
                                prevOne = sdiRecords.getFirst();
                                currOne = (sdiRecords.getFirst()->getNext());
                            }else if ( currOne == sdiRecords.getLast() ){
                                sdiRecords.deleteLast();
                                sdiRecords.deleteLast();
                                currOne = sdiRecords.getLast();
                                prevOne = currOne->getPrev();
                            }else{
                                currOne->getPrev()->getPrev()->setNext(currOne->getNext());
                                currOne->getNext()->setPrev(currOne->getPrev()->getPrev());
                                Data *temp = currOne->getNext();
                                delete(currOne->getPrev());
                                delete(currOne);
                                currOne = temp;
                                prevOne = temp->getPrev();
                            }
                            //std::cout << "844 delete current one and prev" << std::endl;
                        } else if ( currOne->getMapSingleRecord().getChanginglength()>0 && prevOne->getMapSingleRecord().getChanginglength()<0
                                    && currOne->getMapSingleRecord().getAlternative().compare(prevOne->getMapSingleRecord().getReference())==0 &&
                                    (currOne->getMapSingleRecord().getPosition()-1) == prevOne->getMapSingleRecord().getPosition()){
                            //delete current one and prev
                            //std::cout << "850 delete current one and prev" << std::endl;
                            if( ((currOne->getPrev())) == (sdiRecords.getFirst()) ){
                                sdiRecords.deleteFirst();
                                sdiRecords.deleteFirst();
                                prevOne = sdiRecords.getFirst();
                                currOne = (sdiRecords.getFirst()->getNext());
                            }else{
                                currOne->getPrev()->getPrev()->setNext(currOne->getNext());
                                currOne->getNext()->setPrev(currOne->getPrev()->getPrev());
                                Data *temp = currOne->getNext();
                                delete(currOne->getPrev());
                                delete(currOne);
                                currOne = temp;
                                prevOne = temp->getPrev();
                            }
    //                    std::cout << "865 delete current one and prev" << std::endl;
                        } else {
                            prevOne=currOne;
                            currOne=prevOne->getNext();
                        }//std::cout <<  (*itName) << ": link data structure end " << currOne->getMapSingleRecord().getPosition() << std::endl;
                        //std::cout << "2009" << std::endl;
                    }
                }
            }
            //end: merge link data structure
//            std::cout << "line 554" << std::endl;
            std::vector<Variant> sdiRecordsThisOne;
            if( sdiRecords.getFirst()!=NULL ){
                Data* thisone = sdiRecords.getFirst();
                while( thisone!=NULL ){
                    sdiRecordsThisOne.push_back(thisone->getMapSingleRecord());
                    thisone = (thisone->getNext());
                }
            }

            // clear RAM assigned by new Data() begin
            if( sdiRecords.getFirst()!=NULL ){
                Data* thisone = sdiRecords.getFirst();
                while( thisone!=NULL ){
                    Data *tempData = thisone;
                    thisone = (thisone->getNext());
                    delete(tempData);
                }
            }// clear RAM assigned by new Data() end

            std::cout << " begin to sort" << std::endl;
            // transform link to vector and sort and merge nearby records begin
            bool ifChanged = true;
            while(ifChanged){
                std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
                ifChanged = false;
                std::vector<int> sdiRecordsToRomove;
                int oldSize = sdiRecordsThisOne.size();
                for(int j=1; j<oldSize; j++){
                    if(sdiRecordsThisOne[j].getChanginglength()<0 && sdiRecordsThisOne[j-1].getChanginglength()<0 &&
                       (sdiRecordsThisOne[j-1].getPosition()+abs(sdiRecordsThisOne[j-1].getChanginglength()))==sdiRecordsThisOne[j].getPosition()
                       && sdiRecordsThisOne[j-1].getAlternative().compare("-")==0 && sdiRecordsThisOne[j].getAlternative().compare("-")==0 ){
                        int position = sdiRecordsThisOne[j-1].getPosition();
                        std::string ori = sdiRecordsThisOne[j-1].getReference()+sdiRecordsThisOne[j].getReference();
                        std::string result = "-";
                        Variant mapSingleRecord2(it->first, position, ori, result);
                        sdiRecordsToRomove.push_back(j-1);
                        sdiRecordsThisOne[j] = mapSingleRecord2;
                        j++;
                        ifChanged = true;
                    } else if( sdiRecordsThisOne[j].getReference().compare( sdiRecordsThisOne[j].getAlternative())==0){
                        sdiRecordsToRomove.push_back(j); // it does not affect sorting
                    }
                }
                for( int intTpRomoveIndex = sdiRecordsToRomove.size()-1; intTpRomoveIndex>=0 ;--intTpRomoveIndex ){
                    sdiRecordsThisOne.erase(sdiRecordsThisOne.begin()+sdiRecordsToRomove[intTpRomoveIndex]);
                }
            }
            // transform link to vector and sort and merge nearby records end
            newVariantsMaps[it->first]=sdiRecordsThisOne;
        }
    }
    std::ofstream ofile;
    ofile.open(outPutFilePath);
    int itVariantNumebt = 0;
    for( std::map<std::string, std::vector<Variant> >::iterator it=newVariantsMaps.begin();
         it!=newVariantsMaps.end(); ++it ){
        for( std::vector<Variant>::iterator itVariant = it->second.begin();
             itVariant!=it->second.end(); ++itVariant  ){
            if( itVariantNumebt >0 ){
                ofile << std::endl << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" <<
                      itVariant->getChanginglength() << "\t" << itVariant->getReference() << "\t" << itVariant->getAlternative();
            }else{
                ofile << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" <<
                      itVariant->getChanginglength() << "\t" << itVariant->getReference() << "\t" << itVariant->getAlternative();
            }
            ++itVariantNumebt;
        }
    }
    ofile.close();
}
