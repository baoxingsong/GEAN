//
// Created by Baoxing song on 20.10.18.
//

#include "deNovoGenomeVariantCalling.h"


void deNovoGenomeVariantCalling(const std::string & refGffFilePath, const std::string & refFastaFilePath,
                                const std::string & targetGffFilePath, const std::string & targetFastaFilePath,
                                const size_t & minIntron, const size_t & minGene, std::map<std::string, std::string>& parameters,
                                const size_t & widownWidth, const std::string & outPutFilePath) {
    std::ofstream ofile;
    ofile.open(outPutFilePath);

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);

    std::map <std::string, std::vector<Gene>> refGeneMap;
    readGffFileWithEveryThing(refGffFilePath, refGeneMap);
    std::map <std::string, Fasta> refSequences;
    readFastaFile(refFastaFilePath, refSequences);
    updateGeneInformation(refGeneMap, nucleotideCodeSubstitutionMatrix, minIntron, refSequences);
    std::map <std::string, std::set<int32_t>> refToRemove;
    removeDuplication(refGeneMap, minGene, refToRemove);
    for( std::map<std::string, std::vector<Gene> >::iterator it0=refGeneMap.begin(); it0!=refGeneMap.end(); ++it0 ){ //remove non-coding gene and loss-of-function gene
        for( int32_t it1=0; it1<it0->second.size(); ++it1 ){
            if( it0->second[it1].getTranscripts()[0].getCdsVector().size() == 0 || (it0->second[it1].getTranscripts()[0].getIfOrfShift()) ){
                refToRemove[it0->first].insert(it1);
            }
        }
    }
    std::cout << "line 31" << std::endl;
    std::map <std::string, std::vector<int32_t>> refToRemoveVector;
    for( std::map <std::string, std::set<int32_t>>::iterator it0=refToRemove.begin(); it0!=refToRemove.end(); ++it0 ) { //remove non-coding gene and loss-of-function gene
        refToRemoveVector[it0->first]=std::vector<int32_t>();
        for( std::set<int32_t>::iterator it1=it0->second.begin(); it1!=it0->second.end(); ++it1 ){
            refToRemoveVector[it0->first].push_back(*it1);
        }
        std::sort(refToRemoveVector[it0->first].begin(), refToRemoveVector[it0->first].end(), [](int32_t a, int32_t  b) {
            return a < b;
        });
    }
    std::cout << "line 41" << std::endl;
    for( std::map<std::string, std::vector<Gene> >::iterator it0=refGeneMap.begin(); it0!=refGeneMap.end(); ++it0 ){
        for( int j=refToRemoveVector[it0->first].size()-1; j>=0; --j ){
            refGeneMap[it0->first].erase(refGeneMap[it0->first].begin()+refToRemoveVector[it0->first][j]);
        }
    }
    std::cout << "line 48" << std::endl;
    std::map <std::string, std::vector<Gene>> targetGeneMap;
    readGffFileWithEveryThing(targetGffFilePath, targetGeneMap);
    std::map <std::string, Fasta> targetSequences;
    readFastaFile(targetFastaFilePath, targetSequences);
    std::cout << "line 53" << std::endl;
    updateGeneInformation(targetGeneMap, nucleotideCodeSubstitutionMatrix, minIntron, targetSequences);
    std::map <std::string, std::set<int32_t>> targetToRemove;
    std::cout << "line 56" << std::endl;
    removeDuplication(targetGeneMap, minGene, targetToRemove);
    std::cout << "line 58" << std::endl;
    for( std::map<std::string, std::vector<Gene> >::iterator it0=targetGeneMap.begin(); it0!=targetGeneMap.end(); ++it0 ){ //remove non-coding gene and loss-of-function gene
        for( int32_t it1=0; it1<it0->second.size(); ++it1 ){
            if( it0->second[it1].getTranscripts()[0].getCdsVector().size() == 0 || (it0->second[it1].getTranscripts()[0].getIfOrfShift()) ){
                targetToRemove[it0->first].insert(it1);
            }
        }
    }
    std::cout << "line 63" << std::endl;
    std::map <std::string, std::vector<int32_t>> targetToRemoveVector;
    for( std::map <std::string, std::set<int32_t>>::iterator it0=targetToRemove.begin(); it0!=targetToRemove.end(); ++it0 ) { //remove non-coding gene and loss-of-function gene
        targetToRemoveVector[it0->first]=std::vector<int32_t>();
        for( std::set<int32_t>::iterator it1=it0->second.begin(); it1!=it0->second.end(); ++it1 ){
            targetToRemoveVector[it0->first].push_back(*it1);
        }
        std::sort(targetToRemoveVector[it0->first].begin(), targetToRemoveVector[it0->first].end(), [](int32_t a, int32_t  b) {
            return a < b;
        });
    }
    std::cout << "line 74" << std::endl;
    for( std::map<std::string, std::vector<Gene> >::iterator it0=targetGeneMap.begin(); it0!=targetGeneMap.end(); ++it0 ){
        for( int j=targetToRemoveVector[it0->first].size()-1; j>=0; --j ){
            targetGeneMap[it0->first].erase(targetGeneMap[it0->first].begin()+targetToRemoveVector[it0->first][j]);
        }
    }
    std::cout << "line 80" << std::endl;
    for( std::map<std::string, std::vector<Gene> >::iterator it0=targetGeneMap.begin(); it0!=targetGeneMap.end(); ++it0 ){
        int32_t maxScore = - pow(2, 31);
        std::string maxOne;
        std::vector<std::string> alignment_q;
        std::vector<std::string> alignment_d;
        for( std::map<std::string, std::vector<Gene>>::iterator itr=refGeneMap.begin(); itr!=refGeneMap.end(); ++itr ){
            std::vector<std::string> _alignment_q;
            std::vector<std::string> _alignment_d;
            _alignment_q.clear();
            _alignment_d.clear();
            int32_t score = globalAlignment(it0->second, itr->second, _alignment_q, _alignment_d);
            std::cout << "it0: " << it0->first << " itr: " << itr->first << " score: " << score << std::endl;
            if( score > maxScore ){
                maxScore=score;
                maxOne = itr->first;
                alignment_q=_alignment_q;
                alignment_d=_alignment_d;
            }
        }

        std::cout << " it0: " << it0->first << std::endl;
        std::cout << " maxOne: " << maxOne << std::endl;

        size_t startRef = 1;
        size_t startQuery = 1;
        size_t endRef;
        size_t endQuery;
        std::stringstream refAlign;
        std::stringstream queryAlign;
        int32_t refIndex=0;
        int32_t queryIndex=0;
        for(int32_t i=0; i<alignment_q.size(); ++i){
            if( alignment_q[i].compare(alignment_d[i]) == 0 ){
                assert(alignment_q[i].compare(it0->second[queryIndex].getName())==0);
                assert(alignment_d[i].compare(refGeneMap[maxOne][refIndex].getName())==0);
                if(refGeneMap[maxOne][refIndex].getTranscripts()[0].getName().compare( targetGeneMap[it0->first][queryIndex].getTranscripts()[0].getName()) ==0
                    && refGeneMap[maxOne][refIndex].getTranscripts()[0].getPStart()-1 > startRef
                    && targetGeneMap[it0->first][queryIndex].getTranscripts()[0].getPStart()-1 > startQuery ){
                    {
                        endRef = refGeneMap[maxOne][refIndex].getTranscripts()[0].getPStart()-1;
                        endQuery = targetGeneMap[it0->first][queryIndex].getTranscripts()[0].getPStart()-1;
                        std::string refSeq = getSubsequence( refSequences, maxOne, startRef, endRef);
                        std::string querySeq = getSubsequence( targetSequences, it0->first, startQuery, endQuery);
                        /*
                        NeedlemanWunsch_simd_Fast needlemanWunsch_simd_Fast(querySeq, refSeq, match_score, mis_match_score,
                                open_gap_penalty, extend_gap_penalty, nucleotideCodeSubstitutionMatrix);
                        __m256i* vProfile = needlemanWunsch_simd_Fast.query_profile_avx2_byte(nucleotideCodeSubstitutionMatrix);
                        needlemanWunsch_simd_Fast.ssw_avx2(vProfile, nucleotideCodeSubstitutionMatrix);
                        needlemanWunsch_simd_Fast.get_optimal_alignment(nucleotideCodeSubstitutionMatrix);
                        refAlign<<needlemanWunsch_simd_Fast.get_alignment_d();
                        queryAlign<<needlemanWunsch_simd_Fast.get_alignment_q();
                         */
                        std::string _alignment_q;
                        std::string _alignment_d;
                        alignSlidingWindow( querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, parameters, nucleotideCodeSubstitutionMatrix );
                        refAlign<<_alignment_d;
                        queryAlign<<_alignment_q;

                        std::string temp = _alignment_d;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(refSeq)==0);
                        temp = _alignment_q;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(querySeq)==0);
                    }
                    {
                        startRef = refGeneMap[maxOne][refIndex].getTranscripts()[0].getPStart();
                        startQuery= targetGeneMap[it0->first][queryIndex].getTranscripts()[0].getPStart();
                        endRef = refGeneMap[maxOne][refIndex].getTranscripts()[0].getPEnd();
                        endQuery = targetGeneMap[it0->first][queryIndex].getTranscripts()[0].getPEnd();
                        std::string refSeq = getSubsequence( refSequences, maxOne, startRef, endRef);
                        std::string querySeq = getSubsequence( targetSequences, it0->first, startQuery, endQuery);
                        /*
                        NeedlemanWunsch_simd_Fast needlemanWunsch_simd_Fast(querySeq, refSeq, match_score, mis_match_score,
                                                                            open_gap_penalty, extend_gap_penalty, nucleotideCodeSubstitutionMatrix);
                        __m256i* vProfile = needlemanWunsch_simd_Fast.query_profile_avx2_byte(nucleotideCodeSubstitutionMatrix);
                        needlemanWunsch_simd_Fast.ssw_avx2(vProfile, nucleotideCodeSubstitutionMatrix);
                        needlemanWunsch_simd_Fast.get_optimal_alignment(nucleotideCodeSubstitutionMatrix);
                        refAlign<<needlemanWunsch_simd_Fast.get_alignment_d();
                        queryAlign<<needlemanWunsch_simd_Fast.get_alignment_q();
                         */
                        std::string _alignment_q;
                        std::string _alignment_d;
                        alignSlidingWindow( querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, parameters, nucleotideCodeSubstitutionMatrix );
                        refAlign<<_alignment_d;
                        queryAlign<<_alignment_q;

                        std::string temp = _alignment_d;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(refSeq)==0);
                        temp = _alignment_q;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(querySeq)==0);

                    }
                    startRef = refGeneMap[maxOne][refIndex].getTranscripts()[0].getPEnd()+1;
                    startQuery= targetGeneMap[it0->first][queryIndex].getTranscripts()[0].getPEnd()+1;
                }
            }
            if( alignment_q[i].compare("-")!=0 ){
                ++queryIndex;
            }
            if( alignment_d[i].compare("-")!=0 ){
                ++refIndex;
            }
        }
        std::cout << "line 186" << std::endl;
        endRef = refSequences[maxOne].getSequence().length();
        endQuery = targetSequences[it0->first].getSequence().length();
        std::string refSeq = getSubsequence( refSequences, maxOne, startRef, endRef);
        std::string querySeq = getSubsequence( targetSequences, it0->first, startQuery, endQuery);
        /*
        NeedlemanWunsch_simd_Fast needlemanWunsch_simd_Fast(querySeq, refSeq, match_score, mis_match_score,
                                                            open_gap_penalty, extend_gap_penalty, nucleotideCodeSubstitutionMatrix);
        __m256i* vProfile = needlemanWunsch_simd_Fast.query_profile_avx2_byte(nucleotideCodeSubstitutionMatrix);
        needlemanWunsch_simd_Fast.ssw_avx2(vProfile, nucleotideCodeSubstitutionMatrix);
        needlemanWunsch_simd_Fast.get_optimal_alignment(nucleotideCodeSubstitutionMatrix);
        refAlign<<needlemanWunsch_simd_Fast.get_alignment_d();
        queryAlign<<needlemanWunsch_simd_Fast.get_alignment_q();
        */
        std::string _alignment_q;
        std::string _alignment_d;
        alignSlidingWindow( querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, parameters, nucleotideCodeSubstitutionMatrix );
        refAlign<<_alignment_d;
        queryAlign<<_alignment_q;

        std::string temp = _alignment_d;
        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
        assert(temp.compare(refSeq)==0);
        temp = _alignment_q;
        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
        assert(temp.compare(querySeq)==0);

        //std::cout << refAlign.str() << std::endl;
        //std::cout << queryAlign.str() << std::endl;

        temp = refAlign.str();
        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
        assert(temp.compare(refSequences[maxOne].getSequence())==0);
        temp = queryAlign.str();
        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
        assert(temp.compare(targetSequences[it0->first].getSequence())==0);

        std::cout << std::endl;
        std::string queryAlignSeq = queryAlign.str();
        std::string refAlignSeq = refAlign.str();
        std::cout << "line 227 " << std::endl;
        FirstLastList sdiRecords;
        int refLetterNumber = 0;
        for (int ai = 0; ai < refAlignSeq.length(); ai++) {
            //std::cout << "ai: " << ai << std::endl;
            if (refAlignSeq[ai] != '-') {
                ++refLetterNumber;
            }
//            if (refSeq[ai] != resultSeq[ai]) {
            if (refAlignSeq[ai] != queryAlignSeq[ai]) {
                if (queryAlignSeq[ai] == '-') {
//                            std::cout << "1754" << std::endl;
                    //int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber - 1;
                    std::string ori(1, refAlignSeq[ai]);
                    std::string result = "-";
                    Variant mapSingleRecord = Variant(it0->first, refLetterNumber, ori, result);
                    Data *data = new Data(mapSingleRecord);
                    sdiRecords.insertLast(data);
//                            std::cout << "1761" << std::endl;
                } else if (refAlignSeq[ai] == '-') {
//                            std::cout << "1763" << std::endl;
                    if (sdiRecords.getLast() == NULL) {
//                                std::cout << "1755" << std::endl;
                        int position = refLetterNumber+1;
                        std::string ori = "-";
                        std::string result(1, queryAlignSeq[ai]);
                        Variant mapSingleRecord = Variant(it0->first, position, ori, result);
                        Data *data = new Data(mapSingleRecord);
                        sdiRecords.insertLast(data);
//                                std::cout << "1772" << std::endl;
                    } else {
//                                std::cout << "1774" << std::endl;
                        if (NULL != sdiRecords.getLast() &&
                            sdiRecords.getLast()->getMapSingleRecord().getPosition() ==
                            ( refLetterNumber+1)
                            && sdiRecords.getLast()->getMapSingleRecord().getChanginglength() > 0 &&
                            sdiRecords.getLast()->getMapSingleRecord().getReference().compare("-") == 0) {

                            int position = refLetterNumber+1;
                            std::string ori = "-";
                            std::string result =
                                    sdiRecords.getLast()->getMapSingleRecord().getAlternative() + queryAlignSeq[ai];
                            Variant mapSingleRecord = Variant(it0->first, position, ori, result);
                            Data *data = new Data(mapSingleRecord);
                            sdiRecords.deleteLast();
                            sdiRecords.insertLast(data);
                        } else {
                            int position = refLetterNumber+1;
                            std::string ori = "-";
                            std::string result(1, queryAlignSeq[ai]);
                            Variant mapSingleRecord = Variant(it0->first, position, ori, result);
                            Data *data = new Data(mapSingleRecord);
                            sdiRecords.insertLast(data);
                        }
                    }
                } else {
                    int position = refLetterNumber;
                    std::string ori(1, refAlignSeq[ai]);
                    std::string result(1, queryAlignSeq[ai]);
                    Variant mapSingleRecord = Variant(it0->first, position, ori, result);
                    Data *data = new Data(mapSingleRecord);
                    sdiRecords.insertLast(data);
                }
            }
        }

        for( int runingCound = 0; runingCound<2; ++runingCound){
            std::cout << "round " << runingCound << std::endl;
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
                        (prevOne->getMapSingleRecord().getPosition()+abs(prevOne->getMapSingleRecord().getChanginglength()))==currOne->getMapSingleRecord().getPosition()
                        && prevOne->getMapSingleRecord().getAlternative().compare("-")==0 && currOne->getMapSingleRecord().getAlternative().compare("-")==0) { // merge to deletions
                        //std::cout << "797 delete prev" << std::endl;
                        int position = prevOne->getMapSingleRecord().getPosition();
                        std::string ori = prevOne->getMapSingleRecord().getReference() +
                                          currOne->getMapSingleRecord().getReference();
                        std::string result = "-";
                        Variant mapSingleRecord2(it0->first, position, ori, result);
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
                    } else if (currOne->getMapSingleRecord().getChanginglength()==0 && 0 == currOne->getMapSingleRecord().getReference().compare( currOne->getMapSingleRecord().getAlternative())){ // nonsense records
                        //delete current one
                        //std::cout << "820 delete prev" << std::endl;
                        prevOne->setNext( currOne->getNext() );
                        currOne->getNext()->setPrev(prevOne);
                        delete(currOne);
                        currOne = prevOne->getNext();
                        //std::cout << "825 delete prev" << std::endl;
                    } else if ( currOne->getMapSingleRecord().getChanginglength()<0 && prevOne->getMapSingleRecord().getChanginglength()>0
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

        std::cout << it0->first << " link data structure end" << std::endl;
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
                    Variant mapSingleRecord2(it0->first, position, ori, result);
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

        for( std::vector<Variant>::iterator itVariant = sdiRecordsThisOne.begin();
             itVariant!=sdiRecordsThisOne.end(); ++itVariant  ){
            ofile << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" << itVariant->getChanginglength() << "\t" <<
                  itVariant->getReference() << "\t" << itVariant->getAlternative() << std::endl;
        }
    }
    ofile.close();
}
