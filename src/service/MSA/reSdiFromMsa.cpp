//
// Created by baoxing on 10/10/17.
//

#include "reSdiFromMsa.h"

std::mutex gmutexReSdiFromMsa;

void msaFileReadPversion(std::string fileName, std::string& thisFoler, std::map<std::string, std::string>& sdiFiles, std::vector<MsaFileRecord>& msaFileRecords, std::atomic_int& number_of_runing_threads){
    try {
        std::regex reg("^(\\d+)_(\\d+)\\.mafft");

        std::smatch match;
        regex_search(fileName, match, reg);
        if( !match.empty() ) {
            MsaFileRecord msaFileRecord(stoi(match[1]), stoi(match[2]));

            std::string msaFileLocation = thisFoler + "/" + fileName;
            msaFileRead(msaFileRecord, msaFileLocation, sdiFiles);

            gmutexReSdiFromMsa.lock();
            msaFileRecords.push_back(msaFileRecord);
            gmutexReSdiFromMsa.unlock();
        }
    } catch ( ... ){
        --number_of_runing_threads;
        return;
    }
    --number_of_runing_threads;
}

void trimWindow(std::vector<TwoSeqOfMsaResult>& twoSeqOfMsaResults, std::string& accessionName, std::vector<MsaFileRecord>& msaFileRecords,
                std::map<std::string, Fasta>& referenceGenome, std::map<std::string, Fasta>& targetGenome, std::string & chrName,
                int& thisTargetChromosomeLength){
    //delete the extend sequence only keep the wanted region begin
    for( int i=0; i<msaFileRecords.size(); ++i ){
        MsaFileRecord msaFileRecord = msaFileRecords[i];
        int start = msaFileRecord.getStart();
        int end = msaFileRecord.getEnd();
//                --end;
        int targetStart = 0;
        int targetEnd = 0;

        MsaSingleRecord refMsaSingleRecord = msaFileRecord.getMsaSingleRecordRecords()["ref"];
        MsaSingleRecord targetMsaSingleRecord = msaFileRecord.getMsaSingleRecordRecords()[accessionName];
        int msaRefStart = refMsaSingleRecord.getStart();
        int msaTargetStart = targetMsaSingleRecord.getStart();

        //delete the extend sequence only keep the wanted region
        int refLetterNumber = 0;
        int targetLetterNumber = 0;
        std::stringstream refSeq;
        std::stringstream targetSeq;


        //for debug begin
        std::string thisResultSeqS = targetMsaSingleRecord.getSequence();
        thisResultSeqS.erase(std::remove(thisResultSeqS.begin(), thisResultSeqS.end(), '-'), thisResultSeqS.end());
//        std::replace( thisResultSeqS.begin(), thisResultSeqS.end(), "-", '\0');
        //if( targetMsaSingleRecord.getStart() < targetMsaSingleRecord.getEnd() ){

        //}else{
        if( thisResultSeqS.size() != (targetMsaSingleRecord.getEnd()-targetMsaSingleRecord.getStart()+1) ){
            std::cout << "there are problems with the MSA files, please check" << std::endl;
            std::cout << accessionName << " " << thisResultSeqS << " " << thisResultSeqS.size() << " " << targetMsaSingleRecord.getStart() << " " << targetMsaSingleRecord.getEnd() << std::endl <<
            refMsaSingleRecord.getSequence() << " " << refMsaSingleRecord.getStart() << " " << refMsaSingleRecord.getEnd() << std::endl <<
            msaFileRecord.getStart() << " " << msaFileRecord.getEnd() << std::endl;
        }
//        std::cout << accessionName << " " << chrName << " " << thisResultSeqS << " " << thisResultSeqS.size() << " " << targetMsaSingleRecord.getStart() << " " << targetMsaSingleRecord.getEnd() << std::endl;
        assert( thisResultSeqS.size() == (targetMsaSingleRecord.getEnd()-targetMsaSingleRecord.getStart()+1) );
        //}

        std::string thisRefSeqS = refMsaSingleRecord.getSequence();
        thisRefSeqS.erase(std::remove(thisRefSeqS.begin(), thisRefSeqS.end(), '-'), thisRefSeqS.end());

        int thisSubStart1 = refMsaSingleRecord.getStart();
        int thisSubEnd1 = refMsaSingleRecord.getEnd();
        int thisSubStart2 = targetMsaSingleRecord.getStart();
        int thisSubEnd2 = targetMsaSingleRecord.getEnd();

        std::string correctRefSeqs = getSubsequence(referenceGenome, chrName, thisSubStart1, thisSubEnd1);
        std::string correctResultSeqs = getSubsequence(targetGenome, chrName, thisSubStart2, thisSubEnd2);
        if( targetMsaSingleRecord.getStart() > targetMsaSingleRecord.getEnd()){
            correctResultSeqs="";
        }
//                std::cout << "1696 ref: " << thisRefSeqS << " correct: " << correctRefSeqs << " " << refMsaSingleRecord.getStart() << " " << refMsaSingleRecord.getEnd() << std::endl;
//                std::cout << "1697 result:" << thisResultSeqS << "correct: " << correctResultSeqs << " " << targetMsaSingleRecord.getStart() << " " << targetMsaSingleRecord.getEnd() << std::endl;
        if( thisRefSeqS.compare(correctRefSeqs) != 0 ){
            std::cout << "the reference sequence in the MSA files is not the same with the reference genome sequence file" << std::endl;
        }
        if( thisResultSeqS.compare(correctResultSeqs) != 0 ){
            std::cout << "the target genome sequence in the MSA files is not correct, please check" << std::endl;
        }
        assert(thisRefSeqS.compare(correctRefSeqs) == 0);
        assert(thisResultSeqS.compare(correctResultSeqs) == 0);
        //for debug end

        for( int ai=0; ai<refMsaSingleRecord.getSequence().length(); ++ai){
            if( refMsaSingleRecord.getSequence()[ai] != '-' ){
                ++refLetterNumber;
            }
            if( targetMsaSingleRecord.getSequence()[ai] != '-' ){
                ++targetLetterNumber;
            }
            if( start == (msaRefStart+refLetterNumber-1) && refMsaSingleRecord.getSequence()[ai] != '-' ){
                targetStart = msaTargetStart+targetLetterNumber-1;
                if( targetMsaSingleRecord.getSequence()[ai] == '-' ){
                    ++targetStart;
                }
            }
            if( (msaRefStart+refLetterNumber-1) == end && refMsaSingleRecord.getSequence()[ai] != '-'){
                targetEnd = msaTargetStart + targetLetterNumber -1;
            }
            if( start <= (msaRefStart+refLetterNumber-1) && (msaRefStart+refLetterNumber-1)<end ){
                refSeq << refMsaSingleRecord.getSequence()[ai];
                targetSeq << targetMsaSingleRecord.getSequence()[ai];
            }else if ( (msaRefStart+refLetterNumber-1)==end &&  refMsaSingleRecord.getSequence()[ai] != '-' ){
                refSeq << refMsaSingleRecord.getSequence()[ai];
                targetSeq << targetMsaSingleRecord.getSequence()[ai];
            }
        }
        if( ( i == msaFileRecords.size()-1 ) && targetEnd ==0){
            targetEnd = thisTargetChromosomeLength;
            std::cout << "1724 never never run here" << std::endl;
            assert(0);
        }
        std::string refSequence = refSeq.str();
        transform(refSequence.begin(), refSequence.end(), refSequence.begin(),::toupper);
        std::string targetSequence = targetSeq.str();
        transform(targetSequence.begin(), targetSequence.end(), targetSequence.begin(),::toupper);
        if( targetStart==0 && targetEnd==0  ){
            targetEnd = -1;
        }
        TwoSeqOfMsaResult twoSeqOfMsaResult( start, end, refSequence, targetStart, targetEnd, targetSequence );
        twoSeqOfMsaResults.push_back(twoSeqOfMsaResult);
//                std::cout << start << " " << end << std::endl;
//                std::cout << targetStart << " " << targetEnd << std::endl;
        assert(start <= end);
        std::string thisResultSeq = targetSequence;
        thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());
//                std::cout << targetStart << " " << targetEnd << " " << thisResultSeq << std::endl;
        assert( thisResultSeq.size() == (targetEnd-targetStart+1) );


        std::string thisResultSeq2 = twoSeqOfMsaResult.getResultSeq();
        thisResultSeq2.erase(std::remove(thisResultSeq2.begin(), thisResultSeq2.end(), '-'), thisResultSeq2.end());
        assert( thisResultSeq2.size() == (twoSeqOfMsaResult.getResultEnd()-twoSeqOfMsaResult.getResultStart()+1) );
        std::string thisRefSeq = twoSeqOfMsaResult.getRefSeq();
        thisResultSeq = targetSequence;
        thisRefSeq.erase(std::remove(thisRefSeq.begin(), thisRefSeq.end(), '-'), thisRefSeq.end());
        thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());

        int thisSubStart1_1 = twoSeqOfMsaResult.getRefStart();
        int thisSubEnd1_1 = twoSeqOfMsaResult.getRefEnd();

        int thisSubStart2_2 = twoSeqOfMsaResult.getResultStart();
        int thisSubEnd2_2 = twoSeqOfMsaResult.getResultEnd();

        std::string correctRefSeq = getSubsequence(referenceGenome, chrName, thisSubStart1_1, thisSubEnd1_1);
        std::string correctResultSeq = getSubsequence(targetGenome, chrName, thisSubStart2_2, thisSubEnd2_2);
        if( twoSeqOfMsaResult.getResultStart() > twoSeqOfMsaResult.getResultEnd() ){
            correctResultSeq="";
        }
//                std::cout << "1754 ref: " << thisRefSeq << " correct: " << correctRefSeq << " " << twoSeqOfMsaResult.getRefStart() << " " << twoSeqOfMsaResult.getRefEnd() << std::endl;
//                std::cout << "1755 result:" << thisResultSeq << "correct: " << correctResultSeq << " " << twoSeqOfMsaResult.getResultStart() << " " << twoSeqOfMsaResult.getResultEnd() << std::endl;
        assert(thisRefSeq.compare(correctRefSeq) == 0);
        assert(thisResultSeq.compare(correctResultSeq) == 0);
    }
}


void newSdiFileForOneAccession(std::string accessionName, std::map<std::string, std::string>& sdiFiles,
                               std::string & vcfFix, std::map<std::string, Fasta>& referenceGenome,
                               std::map<std::string, std::string>& parameters, std::string & chrName,
                               std::vector<MsaFileRecord>& msaFileRecords, const std::string & outputFolder, std::atomic_int& number_of_runing_threads){

    std::map<std::string, std::vector<Variant> > variantsMaps;
    readSdiFile(sdiFiles[accessionName], variantsMaps, vcfFix, referenceGenome);
    std::map<std::string, Fasta> targetGenome;
    getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome, parameters);
    std::string targetGenomeThisChrSequence = targetGenome[chrName].getSequence();
    int thisTargetChromosomeLength = targetGenomeThisChrSequence.length();

    FirstLastList sdiRecords;

    std::vector<TwoSeqOfMsaResult> twoSeqOfMsaResults;
    trimWindow(twoSeqOfMsaResults, accessionName, msaFileRecords, referenceGenome, targetGenome, chrName, thisTargetChromosomeLength);
    std::sort(twoSeqOfMsaResults.begin(), twoSeqOfMsaResults.end());

    std::cout << accessionName << " delete the extend sequence  done" << std::endl;

    for ( int iIndex = 0; iIndex < twoSeqOfMsaResults.size(); ++iIndex) {
        std::string thisResultSeq2 = twoSeqOfMsaResults[iIndex].getResultSeq();
        thisResultSeq2.erase(std::remove(thisResultSeq2.begin(), thisResultSeq2.end(), '-'), thisResultSeq2.end());
        assert( thisResultSeq2.size() == (twoSeqOfMsaResults[iIndex].getResultEnd()-twoSeqOfMsaResults[iIndex].getResultStart()+1) );
        std::string thisRefSeq = twoSeqOfMsaResults[iIndex].getRefSeq();
        std::string thisResultSeq = twoSeqOfMsaResults[iIndex].getResultSeq();
        thisRefSeq.erase(std::remove(thisRefSeq.begin(), thisRefSeq.end(), '-'), thisRefSeq.end());
        thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());
        int thisSubStart_1 = twoSeqOfMsaResults[iIndex].getRefStart();
        int thisSubEnd_1 = twoSeqOfMsaResults[iIndex].getRefEnd();

        int thisSubStart_2 = twoSeqOfMsaResults[iIndex].getResultStart();
        int thisSubEnd_2 = twoSeqOfMsaResults[iIndex].getResultEnd();
        std::string correctRefSeq = getSubsequence(referenceGenome, chrName, thisSubStart_1, thisSubEnd_1);
        std::string correctResultSeq = getSubsequence(targetGenome, chrName, thisSubStart_2, thisSubEnd_2);
        if( twoSeqOfMsaResults[iIndex].getResultStart() > twoSeqOfMsaResults[iIndex].getResultEnd() ){
            correctResultSeq="";
        }
//                std::cout << "1737ref: " << thisRefSeq << " correct: " << correctRefSeq << " " << twoSeqOfMsaResults[iIndex].getRefStart() << " " << twoSeqOfMsaResults[iIndex].getRefEnd() << std::endl;
//                std::cout << "1738 result:" << thisResultSeq << "correct: " << correctResultSeq << " " << twoSeqOfMsaResults[iIndex].getResultStart() << " " << twoSeqOfMsaResults[iIndex].getResultEnd() << std::endl;
        assert(thisRefSeq.compare(correctRefSeq) == 0);
        assert(thisResultSeq.compare(correctResultSeq) == 0);
    }

    // merge overlapped neighbor window begin
    for( int j=1; j<twoSeqOfMsaResults.size(); j++ ){
        //if( twoSeqOfMsaResults.size() > 1 ){
        if( j<1 ){
            j=1;
        }
        std::string thisResultSeq = twoSeqOfMsaResults[j-1].getResultSeq();
        thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());
        assert( thisResultSeq.size() == (twoSeqOfMsaResults[j-1].getResultEnd()-twoSeqOfMsaResults[j-1].getResultStart()+1) );

        thisResultSeq = twoSeqOfMsaResults[j].getResultSeq();
        thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());
        //std::cout << accessionName << " " << chrName << " " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << " " << thisResultSeq << std::endl;
        assert( thisResultSeq.size() == (twoSeqOfMsaResults[j].getResultEnd()-twoSeqOfMsaResults[j].getResultStart()+1) );

        if( twoSeqOfMsaResults[j].getResultEnd() < twoSeqOfMsaResults[j].getResultStart() ){
            //std::cout <<"1732 to be complete , This should never appear. If you find it, please contact me." << std::endl;
//                        std::cout << "j: " << j <<std::endl;
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::string thisResultSeq = twoSeqOfMsaResults[j].getResultSeq();
//                        std::cout << "j: " << thisResultSeq <<std::endl;

            if ( (twoSeqOfMsaResults[j-1].getRefEnd() >= twoSeqOfMsaResults[j].getRefStart()) ){
                int overLapLength = twoSeqOfMsaResults[j-1].getRefEnd() - twoSeqOfMsaResults[j].getRefStart();
                overLapLength++;
                std::string newRef = twoSeqOfMsaResults[j].getRefSeq();
                int corrected = 0;
                for( int t=0; t<newRef.length(); t++ ){
                    if(corrected<overLapLength && newRef[t]!='-' ){
                        newRef=newRef.substr(0, t)+"-"+newRef.substr(t+1, newRef.length()-t-1);
                        corrected++;
                    }
                }
                int largerRefEnd;
                if( twoSeqOfMsaResults[j].getRefEnd() > twoSeqOfMsaResults[j-1].getRefEnd() ){
                    largerRefEnd = twoSeqOfMsaResults[j].getRefEnd();
                }else{
                    largerRefEnd = twoSeqOfMsaResults[j-1].getRefEnd();
                }
//
//                            std::cout << "j: " << j <<std::endl;
//                            std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                            std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;

                twoSeqOfMsaResults[j-1].setRefEnd(largerRefEnd);
                twoSeqOfMsaResults[j-1].setRefSeq(twoSeqOfMsaResults[j-1].getRefSeq()+newRef);

                twoSeqOfMsaResults[j-1].setResultSeq(twoSeqOfMsaResults[j-1].getResultSeq()+twoSeqOfMsaResults[j].getResultSeq());
//                            std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                            std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
                twoSeqOfMsaResults.erase(twoSeqOfMsaResults.begin()+j);
//                            std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                            std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                            std::cout <<"1762 to be complete , This should never appear. If you find it, please contact me." << std::endl;
                --j;
                --j;
            } else if ( twoSeqOfMsaResults[j-1].getResultEnd() < twoSeqOfMsaResults[j-1].getResultStart()){
                twoSeqOfMsaResults[j].setResultStart(twoSeqOfMsaResults[j-1].getResultStart());
                twoSeqOfMsaResults[j].setResultEnd(twoSeqOfMsaResults[j-1].getResultStart()-1);
            } else {
                twoSeqOfMsaResults[j].setResultStart(twoSeqOfMsaResults[j-1].getResultEnd());
                twoSeqOfMsaResults[j].setResultEnd(twoSeqOfMsaResults[j-1].getResultEnd()-1);
            }
        } else if ( (twoSeqOfMsaResults[j-1].getResultEnd() >= twoSeqOfMsaResults[j].getResultStart()) && (twoSeqOfMsaResults[j-1].getRefEnd() >= twoSeqOfMsaResults[j].getRefStart())
                    && twoSeqOfMsaResults[j-1].getResultStart() <= twoSeqOfMsaResults[j-1].getResultEnd() ){
            int overLapLength = twoSeqOfMsaResults[j-1].getResultEnd() - twoSeqOfMsaResults[j].getResultStart();
            ++overLapLength;
            std::string newSequence = twoSeqOfMsaResults[j].getResultSeq();
            int corrected = 0;
            for( int t=0; t<newSequence.length(); t++ ){
                if(corrected<overLapLength && newSequence[t]!='-' ){
                    newSequence=newSequence.substr(0, t)+"-"+newSequence.substr(t+1, newSequence.length()-t-1);
                    corrected++;
                }
            }
            int largerEnd;
            if( twoSeqOfMsaResults[j].getResultEnd() > twoSeqOfMsaResults[j-1].getResultEnd() ){
                largerEnd = twoSeqOfMsaResults[j].getResultEnd();
            } else {
                largerEnd = twoSeqOfMsaResults[j-1].getResultEnd();
            }

            overLapLength = twoSeqOfMsaResults[j-1].getRefEnd() - twoSeqOfMsaResults[j].getRefStart();
            overLapLength++;
            std::string newRef = twoSeqOfMsaResults[j].getRefSeq();
            corrected = 0;
            for( int t=0; t<newRef.length(); t++ ){
                if(corrected<overLapLength && newRef[t]!='-' ){
                    newRef=newRef.substr(0, t)+"-"+newRef.substr(t+1, newRef.length()-t-1);
                    corrected++;
                }
            }
            int largerRefEnd;
            if( twoSeqOfMsaResults[j].getRefEnd() > twoSeqOfMsaResults[j-1].getRefEnd() ){
                largerRefEnd = twoSeqOfMsaResults[j].getRefEnd();
            }else{
                largerRefEnd = twoSeqOfMsaResults[j-1].getRefEnd();
            }
//                        std::cout << "j: " << j <<std::endl;
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;

            twoSeqOfMsaResults[j-1].setRefEnd(largerRefEnd);
            twoSeqOfMsaResults[j-1].setRefSeq(twoSeqOfMsaResults[j-1].getRefSeq()+newRef);

            twoSeqOfMsaResults[j-1].setResultEnd(largerEnd);
            twoSeqOfMsaResults[j-1].setResultSeq(twoSeqOfMsaResults[j-1].getResultSeq()+newSequence);
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
            twoSeqOfMsaResults.erase(twoSeqOfMsaResults.begin()+j);
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout <<"1777 to be complete , This should never appear. If you find it, please contact me." << std::endl;
            --j;
            --j;
        }else if ( (twoSeqOfMsaResults[j-1].getResultStart() >= twoSeqOfMsaResults[j].getResultStart()) && (twoSeqOfMsaResults[j-1].getRefEnd() >= twoSeqOfMsaResults[j].getRefStart())
                   && twoSeqOfMsaResults[j-1].getResultStart() > twoSeqOfMsaResults[j-1].getResultEnd() ){
            int overLapLength = twoSeqOfMsaResults[j-1].getRefEnd() - twoSeqOfMsaResults[j].getRefStart();
            overLapLength++;
            std::string newRef = twoSeqOfMsaResults[j].getRefSeq();
            int corrected = 0;
            for( int t=0; t<newRef.length(); t++ ){
                if(corrected<overLapLength && newRef[t]!='-' ){
                    newRef=newRef.substr(0, t)+"-"+newRef.substr(t+1, newRef.length()-t-1);
                    corrected++;
                }
            }
            int largerRefEnd;
            if( twoSeqOfMsaResults[j].getRefEnd() > twoSeqOfMsaResults[j-1].getRefEnd() ){
                largerRefEnd = twoSeqOfMsaResults[j].getRefEnd();
            }else{
                largerRefEnd = twoSeqOfMsaResults[j-1].getRefEnd();
            }

            overLapLength = twoSeqOfMsaResults[j-1].getResultStart() - twoSeqOfMsaResults[j].getResultStart();
            ++overLapLength;
            std::string newSequence = twoSeqOfMsaResults[j].getResultSeq();
            corrected = 0;
            for( int t=0; t<newSequence.length(); t++ ){
                if(corrected<overLapLength && newSequence[t]!='-' ){
                    newSequence=newSequence.substr(0, t)+"-"+newSequence.substr(t+1, newSequence.length()-t-1);
                    corrected++;
                }
            }
            int end = twoSeqOfMsaResults[j].getResultEnd(); // the current one is not empty
            int start =twoSeqOfMsaResults[j].getResultStart() + overLapLength;


            twoSeqOfMsaResults[j-1].setRefEnd(largerRefEnd);
            twoSeqOfMsaResults[j-1].setRefSeq(twoSeqOfMsaResults[j-1].getRefSeq()+newRef);

            if( start < end ){
                twoSeqOfMsaResults[j-1].setResultEnd(end);
                twoSeqOfMsaResults[j-1].setResultStart(start);
                twoSeqOfMsaResults[j-1].setResultSeq(twoSeqOfMsaResults[j-1].getResultSeq()+newSequence);
            } else {
                twoSeqOfMsaResults[j-1].setResultSeq(twoSeqOfMsaResults[j-1].getResultSeq()+twoSeqOfMsaResults[j].getResultSeq());
            }

//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
            twoSeqOfMsaResults.erase(twoSeqOfMsaResults.begin()+j);

//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout <<"line 1952" << std::endl;
            --j;
            --j;

        } else if( (twoSeqOfMsaResults[j-1].getResultEnd() >= twoSeqOfMsaResults[j].getResultStart()) && twoSeqOfMsaResults[j-1].getResultStart() <= twoSeqOfMsaResults[j-1].getResultEnd() ) {

            int overLapLength =  twoSeqOfMsaResults[j - 1].getResultEnd() - twoSeqOfMsaResults[j].getResultStart();
            ++overLapLength;
            std::string newSequence = twoSeqOfMsaResults[j].getResultSeq();
            int corrected = 0;
            for (int t = 0; t < newSequence.length(); ++t) {
                if (corrected < overLapLength && newSequence[t] != '-') {
                    newSequence = newSequence.substr(0, t) + "-" + newSequence.substr(t + 1, newSequence.length() - t - 1);
                    corrected++;
                }
            }
//
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
            twoSeqOfMsaResults[j].setResultSeq(newSequence);
            twoSeqOfMsaResults[j].setResultStart(twoSeqOfMsaResults[j].getResultStart()+overLapLength);
            if( twoSeqOfMsaResults[j].getResultStart() > twoSeqOfMsaResults[j].getResultEnd() ){
                twoSeqOfMsaResults[j].setResultStart( twoSeqOfMsaResults[j-1].getResultEnd() );
                twoSeqOfMsaResults[j].setResultEnd(twoSeqOfMsaResults[j].getResultStart() - 1);
            }
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout <<"line 1985" << std::endl;
        } else if( (twoSeqOfMsaResults[j-1].getResultStart()>=twoSeqOfMsaResults[j].getResultStart()) && twoSeqOfMsaResults[j-1].getResultStart()>twoSeqOfMsaResults[j-1].getResultEnd() ){

            int overLapLength =  twoSeqOfMsaResults[j - 1].getResultStart() - twoSeqOfMsaResults[j].getResultStart();

            ++overLapLength;
            std::string newSequence = twoSeqOfMsaResults[j].getResultSeq();
            int corrected = 0;
            for (int t = 0; t < newSequence.length(); ++t) {
                if (corrected < overLapLength && newSequence[t] != '-') {
                    newSequence = newSequence.substr(0, t) + "-" + newSequence.substr(t + 1, newSequence.length() - t - 1);
                    corrected++;
                }
            }
//
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
            twoSeqOfMsaResults[j].setResultSeq(newSequence);
            twoSeqOfMsaResults[j].setResultStart(twoSeqOfMsaResults[j].getResultStart()+overLapLength);
            if( twoSeqOfMsaResults[j].getResultStart() > twoSeqOfMsaResults[j].getResultEnd() ){
                twoSeqOfMsaResults[j].setResultStart( twoSeqOfMsaResults[j-1].getResultStart() );
                twoSeqOfMsaResults[j].setResultEnd(twoSeqOfMsaResults[j].getResultStart() - 1);
            }

//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout <<"merged 2008" << std::endl;

        }else if ( (twoSeqOfMsaResults[j-1].getRefEnd() >= twoSeqOfMsaResults[j].getRefStart()) ){
            int overLapLength = twoSeqOfMsaResults[j-1].getRefEnd() - twoSeqOfMsaResults[j].getRefStart();
            overLapLength++;
            std::string newRef = twoSeqOfMsaResults[j].getRefSeq();
            int corrected = 0;
            for( int t=0; t<newRef.length(); t++ ){
                if(corrected<overLapLength && newRef[t]!='-'){
                    newRef=newRef.substr(0, t)+"-"+newRef.substr(t+1, newRef.length()-t-1);
                    corrected++;
                }
            }
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
            twoSeqOfMsaResults[j].setRefSeq(newRef);
            twoSeqOfMsaResults[j].setRefStart(twoSeqOfMsaResults[j].getRefStart()+overLapLength);
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout <<"1824 to be complete , This should never appear. If you find it, please contact me." << std::endl;
        }else if( twoSeqOfMsaResults[j].getResultStart() < twoSeqOfMsaResults[j-1].getResultStart() ){
//                        std::cout << "j: " << j <<std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << " " << twoSeqOfMsaResults[j-1].getRefStart() <<" "<<  twoSeqOfMsaResults[j-1].getRefEnd() << std::endl;
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << " " << twoSeqOfMsaResults[j].getRefStart() << " " << twoSeqOfMsaResults[j].getRefEnd() << std::endl;
            std::cout <<"1987 to be complete , This should never appear. If you find it, please contact me." << std::endl; // this should never run
        } else {
//                        std::cout << "j: " << j <<std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout <<"1781 to be complete. Check it there any problem." << std::endl;
        }
        //}
    } // merge overlapped neighbor window end
    std::sort(twoSeqOfMsaResults.begin(), twoSeqOfMsaResults.end());

    for ( int iIndex = 0; iIndex < twoSeqOfMsaResults.size(); ++iIndex) {

        std::string thisResultSeq2 = twoSeqOfMsaResults[iIndex].getResultSeq();
        thisResultSeq2.erase(std::remove(thisResultSeq2.begin(), thisResultSeq2.end(), '-'), thisResultSeq2.end());
        assert( thisResultSeq2.size() == (twoSeqOfMsaResults[iIndex].getResultEnd()-twoSeqOfMsaResults[iIndex].getResultStart()+1) );

        std::string thisRefSeq = twoSeqOfMsaResults[iIndex].getRefSeq();
        std::string thisResultSeq = twoSeqOfMsaResults[iIndex].getResultSeq();

        thisRefSeq.erase(std::remove(thisRefSeq.begin(), thisRefSeq.end(), '-'), thisRefSeq.end());
        thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());
        int thisSubStart1 = twoSeqOfMsaResults[iIndex].getRefStart();
        int thisSubEnd1 = twoSeqOfMsaResults[iIndex].getRefEnd();

        int thisSubStart2 = twoSeqOfMsaResults[iIndex].getResultStart();
        int thisSubEnd2 = twoSeqOfMsaResults[iIndex].getResultEnd();

        std::string correctRefSeq = getSubsequence(referenceGenome, chrName, thisSubStart1, thisSubEnd1);
        std::string correctResultSeq = getSubsequence(targetGenome, chrName, thisSubStart2, thisSubEnd2);

        if( twoSeqOfMsaResults[iIndex].getResultStart() > twoSeqOfMsaResults[iIndex].getResultEnd() ){
            correctResultSeq="";
        }
//        std::cout << "1986 ref: " << thisRefSeq << " correct: " << correctRefSeq << " " << twoSeqOfMsaResults[iIndex].getRefStart() << " " << twoSeqOfMsaResults[iIndex].getRefEnd() << std::endl;
//        std::cout << "1987 result:" << thisResultSeq << "correct: " << correctResultSeq << " " << twoSeqOfMsaResults[iIndex].getResultStart() << " " << twoSeqOfMsaResults[iIndex].getResultEnd() << std::endl;
        assert(thisRefSeq.compare(correctRefSeq) == 0);
        assert(thisResultSeq.compare(correctResultSeq) == 0);
    }

    std::cout << accessionName << " begin: insert the SDI record before the first window" << std::endl;



    if( twoSeqOfMsaResults[0].getRefStart()>1){//begin: insert the SDI record before the first window
        std::cout << accessionName << " " << chrName << " start " << twoSeqOfMsaResults[0].getRefStart() << std::endl;
        int thisSubStart1 = 1;
        int thisSubEnd1 = twoSeqOfMsaResults[0].getRefStart()-1;
        std::string ori = getSubsequence(referenceGenome, chrName, thisSubStart1, thisSubEnd1);
        std::string result;
        if(twoSeqOfMsaResults[0].getResultStart()>1){
            int thisSubStart = 1;
            int thisSubEnd = twoSeqOfMsaResults[0].getResultStart()-1;

            result = getSubsequence( targetGenome, chrName, thisSubStart, thisSubEnd);
        }else{
            result = "-";
        }
        int position = 1;

        if( ori.compare(result) != 0 ){
            Variant mapSingleRecord = Variant(chrName, position, ori, result);
            Data* data = new Data(mapSingleRecord);
            sdiRecords.insertLast(data);
        }
    }else if(twoSeqOfMsaResults[0].getRefStart()==1 && twoSeqOfMsaResults[0].getResultStart()>1){
        int thisSubStart = 1;
        int thisSubEnd = twoSeqOfMsaResults[0].getResultStart()-1;

        std::string result = getSubsequence(targetGenome, chrName, thisSubStart, thisSubEnd);
        int position = 1;
        std::string ori = "-";
        Variant mapSingleRecord = Variant(chrName, position, ori, result);
        Data* data = new Data(mapSingleRecord);
        sdiRecords.insertLast(data);
    } //end: insert the first SDI record before the first windows


    std::cout << accessionName << " begin: insert the windows region SDI record" << std::endl;

    //begin: insert the windows region SDI record
    int i=0;
    while(i<(twoSeqOfMsaResults.size()-1)) {
        std::string refSeq = twoSeqOfMsaResults[i].getRefSeq();
        std::string resultSeq = twoSeqOfMsaResults[i].getResultSeq();
        //std::cout << "ref:       " << refSeq << std::endl;
        //std::cout << "resultSeq: " << resultSeq << std::endl;
        int refLetterNumber = 0;
        for (int ai = 0; ai < refSeq.length(); ai++) {
            if (refSeq[ai] != '-') {
                ++refLetterNumber;
            }
//            if (refSeq[ai] != resultSeq[ai]) {
                if (refSeq[ai] != resultSeq[ai]) {
                    if (resultSeq[ai] == '-') {
//                            std::cout << "1754" << std::endl;
                        int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber - 1;
                        std::string ori(1, refSeq[ai]);
                        std::string result = "-";
                        Variant mapSingleRecord = Variant(chrName, position, ori, result);
                        Data* data = new Data(mapSingleRecord);
                        sdiRecords.insertLast(data);
//                            std::cout << "1761" << std::endl;
                    } else if (refSeq[ai] == '-') {
//                            std::cout << "1763" << std::endl;
                        if (sdiRecords.getLast() == NULL) {
//                                std::cout << "1755" << std::endl;
                            int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                            std::string ori = "-";
                            std::string result(1, resultSeq[ai]);
                            Variant mapSingleRecord = Variant(chrName, position, ori, result);
                            Data* data = new Data(mapSingleRecord);
                            sdiRecords.insertLast(data);
//                                std::cout << "1772" << std::endl;
                        } else {
//                                std::cout << "1774" << std::endl;
                            if (NULL!=sdiRecords.getLast() && sdiRecords.getLast()->getMapSingleRecord().getPosition()==(twoSeqOfMsaResults[i].getRefStart() + refLetterNumber)
                                && sdiRecords.getLast()->getMapSingleRecord().getChanginglength() > 0 &&
                                sdiRecords.getLast()->getMapSingleRecord().getReference().compare("-") == 0) {

                                int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                                std::string ori = "-";
                                std::string result = sdiRecords.getLast()->getMapSingleRecord().getAlternative() + resultSeq[ai];
                                Variant mapSingleRecord = Variant(chrName, position, ori, result);
                                Data* data = new Data(mapSingleRecord);
                                sdiRecords.deleteLast();
                                sdiRecords.insertLast(data);
                            } else {
                                int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                                std::string ori = "-";
                                std::string result(1, resultSeq[ai]);
                                Variant mapSingleRecord = Variant(chrName, position, ori, result);
                                Data* data = new Data(mapSingleRecord);
                                sdiRecords.insertLast(data);
                            }
//                                std::cout << "1793" << std::endl;
                        }
                    } else {
//                            std::cout << "1796" << std::endl;
                        int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber - 1;
                        std::string ori(1, refSeq[ai]);
                        std::string result(1, resultSeq[ai]);
                        Variant mapSingleRecord = Variant(chrName, position, ori, result);
                        Data* data = new Data(mapSingleRecord);
                        sdiRecords.insertLast(data);
//                            std::cout << "1803" << std::endl;
                    }
                }
//            }
        }
        {
            //begin: insert the inter window region SDI record
            Variant mapSingleRecord;
            std::string oriSeq;
//                    std::cout << "1811" << std::endl;
            if ((twoSeqOfMsaResults[i].getRefEnd()) < (twoSeqOfMsaResults[i + 1].getRefStart() - 1)) {
                int thisSubStart = twoSeqOfMsaResults[i].getRefEnd() + 1;
                int thisSubEnd = twoSeqOfMsaResults[i + 1].getRefStart() - 1;
                oriSeq = getSubsequence(referenceGenome, chrName, thisSubStart, thisSubEnd);
            } else {
                oriSeq = "-";
            }
            std::string resultSeq;
//                    int resultSeqLength = 0;
            if ( (twoSeqOfMsaResults[i].getResultEnd() >= (twoSeqOfMsaResults[i].getResultStart()))  && (twoSeqOfMsaResults[i].getResultEnd() < (twoSeqOfMsaResults[i + 1].getResultStart() - 1)) ) {
                int thisSubStart = twoSeqOfMsaResults[i].getResultEnd() + 1;
                int thisSubEnd = twoSeqOfMsaResults[i + 1].getResultStart() - 1;
                resultSeq = getSubsequence(targetGenome, chrName, thisSubStart,thisSubEnd);
//                        resultSeqLength = resultSeq.length();
            } else if( (twoSeqOfMsaResults[i].getResultEnd() < (twoSeqOfMsaResults[i].getResultStart())) &&  (twoSeqOfMsaResults[i].getResultStart() < (twoSeqOfMsaResults[i + 1].getResultStart() - 1)) ){
                int thisSubStart = twoSeqOfMsaResults[i].getResultStart() + 1;
                int thisSubEnd = twoSeqOfMsaResults[i + 1].getResultStart() - 1;
                resultSeq = getSubsequence(targetGenome, chrName, thisSubStart, thisSubEnd);
            }// if the one is totally deleted
            else {
                resultSeq = "-";
            }
            if ( (oriSeq.compare(resultSeq)!=0) && (oriSeq.compare("-")==0) ){
                int position = twoSeqOfMsaResults[i].getRefEnd() + 1;
                mapSingleRecord = Variant(chrName, position, oriSeq, resultSeq);
//                        std::cout << "2140 " << chrName << " " << position << " " << oriSeq << " " << resultSeq << std::endl;
                Data* data = new Data (mapSingleRecord);
                sdiRecords.insertLast(data);
            }else if (oriSeq.compare(resultSeq) != 0) {
//                        std::cout << "1828" << std::endl;
                int position = twoSeqOfMsaResults[i].getRefEnd() + 1;
                mapSingleRecord = Variant(chrName, position, oriSeq, resultSeq);
//                        std::cout << "1981 " << chrName << " " << position << " " << oriSeq << " " << resultSeq << std::endl;
                Data* data = new Data (mapSingleRecord);
                sdiRecords.insertLast(data);
            }
            i++;
        }
        refLetterNumber =0;
    }
    //end: insert the windows region SDI record
    std::cout << "begin: insert the last windows region SDI record" << std::endl;


    //begin: insert the last windows region SDI record
    if(i == (twoSeqOfMsaResults.size()-1)){
        //System.out.println("i="+i+" twoSeqOfMsaResults.size()-1:"+(twoSeqOfMsaResults.size()-1));
        int refLetterNumber = 0;
        std::string refSeq = twoSeqOfMsaResults[i].getRefSeq();
        std::string resultSeq = twoSeqOfMsaResults[i].getResultSeq();
        for(size_t ai=0; ai<refSeq.length(); ai++){
            if (refSeq[ai] != '-') {
                ++refLetterNumber;
            }
            if (refSeq[ai] != resultSeq[ai]) {
                if (resultSeq[ai] == '-') {
//                            std::cout << "1754" << std::endl;
                    int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber - 1;
                    std::string ori(1, refSeq[ai]);
                    std::string result = "-";
                    Variant mapSingleRecord = Variant(chrName, position, ori, result);
                    Data* data = new Data(mapSingleRecord);
                    sdiRecords.insertLast(data);
//                            std::cout << "1761" << std::endl;
                } else if (refSeq[ai] == '-') {
//                            std::cout << "1763" << std::endl;
                    if (sdiRecords.getLast() == NULL) {
//                                std::cout << "1755" << std::endl;
                        int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                        std::string ori = "-";
                        std::string result(1, resultSeq[ai]);
                        Variant mapSingleRecord = Variant(chrName, position, ori, result);
                        Data* data = new Data(mapSingleRecord);
                        sdiRecords.insertLast(data);
//                                std::cout << "1772" << std::endl;
                    } else {
//                                std::cout << "1774" << std::endl;
                        if (NULL!=sdiRecords.getLast() && sdiRecords.getLast()->getMapSingleRecord().getPosition()==(twoSeqOfMsaResults[i].getRefStart() + refLetterNumber)
                            && sdiRecords.getLast()->getMapSingleRecord().getChanginglength() > 0 &&
                            sdiRecords.getLast()->getMapSingleRecord().getReference().compare("-") == 0) {

                            int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                            std::string ori = "-";
                            std::string result = sdiRecords.getLast()->getMapSingleRecord().getAlternative() + resultSeq[ai];
                            Variant mapSingleRecord = Variant(chrName, position, ori, result);
                            Data* data = new Data(mapSingleRecord);
                            sdiRecords.deleteLast();
                            sdiRecords.insertLast(data);
                        } else {
                            int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                            std::string ori = "-";
                            std::string result(1, resultSeq[ai]);
                            Variant mapSingleRecord = Variant(chrName, position, ori, result);
                            Data* data = new Data(mapSingleRecord);
                            sdiRecords.insertLast(data);
                        }
//                                std::cout << "1793" << std::endl;
                    }
                } else {
//                            std::cout << "1796" << std::endl;
                    int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber - 1;
                    std::string ori(1, refSeq[ai]);
                    std::string result(1, resultSeq[ai]);
                    Variant mapSingleRecord = Variant(chrName, position, ori, result);
                    Data* data = new Data(mapSingleRecord);
                    sdiRecords.insertLast(data);
//                            std::cout << "1803" << std::endl;
                }
            }
        }
    }
    //end: insert the last windows region SDI record
    std::cout << "insert the SDI record after last window" << std::endl;
    //begin: insert the SDI record after last window
    int endIndex = twoSeqOfMsaResults.size()-1;
    std::string oriSeq;
    if((twoSeqOfMsaResults[endIndex].getRefEnd()) < referenceGenome[chrName].getSequence().length()){
        int thisSubStart = twoSeqOfMsaResults[endIndex].getRefEnd()+1;
        int thisSubEnd = referenceGenome[chrName].getSequence().length();
        oriSeq = getSubsequence(referenceGenome, chrName, thisSubStart, thisSubEnd);
    }else{
        oriSeq="-";
    }
    if(twoSeqOfMsaResults[endIndex].getResultEnd() < (targetGenome[chrName].getSequence().length()-1)){
        int thisSubStart = twoSeqOfMsaResults[endIndex].getResultEnd() + 1;
        int thisSubEnd = targetGenome[chrName].getSequence().length();
        std::string resultSeq = getSubsequence(targetGenome, chrName, thisSubStart, thisSubEnd);
        if(oriSeq.compare(resultSeq) != 0){
            int position = twoSeqOfMsaResults[endIndex].getRefEnd()+1;
//                    std::cout << "1917" << std::endl;
            Variant mapSingleRecord = Variant(chrName, position, oriSeq, resultSeq);
            Data* data = new Data(mapSingleRecord);
            sdiRecords.insertLast(data);
        }
    }else if(twoSeqOfMsaResults[endIndex].getResultEnd() == (targetGenome[chrName].getSequence().length()-1)){
        int position = twoSeqOfMsaResults[endIndex].getRefEnd()+1;
        std::string result( 1, targetGenome[chrName].getSequence()[targetGenome[chrName].getSequence().length()-1] );
//                std::cout << "1925" << std::endl;
        if(oriSeq.compare(result) != 0) {
            Variant mapSingleRecord = Variant(chrName, position, oriSeq, result);
            Data *data = new Data(mapSingleRecord);
            sdiRecords.insertLast(data);
        }
    }else if(twoSeqOfMsaResults[endIndex].getResultEnd() == (targetGenome[chrName].getSequence().length())  && oriSeq.compare("-")!=0 ){
        int position = twoSeqOfMsaResults[endIndex].getRefEnd()+1;
        std::string result = "-";
//                std::cout << "1932" << std::endl;
        Variant mapSingleRecord = Variant(chrName, position, oriSeq, result);
        Data* data = new Data(mapSingleRecord);
        sdiRecords.insertLast(data);
    }
//    else if(oriSeq.length()>0){
//        std::cout << "should never run here 2333" << std::endl;
//    }
    //end: insert the SDI record after last window

    twoSeqOfMsaResults.clear();// for RAM saving
    //begin: merge link data structure
    std::cout << accessionName << " link data structure start" << std::endl;



//    --number_of_runing_threads;
//    return;


//            if( (sdiRecords[chrName].getFirst() != NULL)
//                && (sdiRecords[chrName].getFirst()->getNext() != NULL ) ){
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
                    Variant mapSingleRecord2(chrName, position, ori, result);
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

    std::cout << accessionName << " link data structure end" << std::endl;
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
                Variant mapSingleRecord2(chrName, position, ori, result);
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

    std::ofstream ofile;
    ofile.open(outputFolder + "/" + accessionName + ".sdi", std::ofstream::app | std::ofstream::out);
    for( std::vector<Variant>::iterator itVariant = sdiRecordsThisOne.begin();
        itVariant!=sdiRecordsThisOne.end(); ++itVariant  ){
        ofile << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" << itVariant->getChanginglength() << "\t" <<
                  itVariant->getReference() << "\t" << itVariant->getAlternative() << std::endl;
    }
    ofile.close();
    --number_of_runing_threads;
}


void constructSdiFromMsa(std::vector<std::string>& chromosomes, std::string& folder, std::string& outputFolder, std::string & referenceGenomeFilePath,
                         std::map<std::string, std::string>& sdiFiles, std::string & vcfFix, std::map<std::string, std::string>& parameters, int & maxThread){
    std::map<std::string, Fasta> referenceGenome;
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    readFastaFile(referenceGenomeFilePath, referenceGenome);

    int len = folder.length();
    if (folder[len-1] != '/'){
        folder = folder + "/";
    }

    //prepare output file begin
    createdFolder(outputFolder);
    for( std::map<std::string, std::string>::iterator itName=sdiFiles.begin(); itName!=sdiFiles.end(); ++itName ) {
        std::ofstream ofile;
        ofile.open(outputFolder + "/" + itName->first + ".sdi");
        ofile.close();
    }// prepare output file end

    for( size_t chromosomei=0; chromosomei<chromosomes.size(); ++chromosomei ){
        // read file begin
        std::string chrName = chromosomes[chromosomei];

        createdFolder(outputFolder + "/" + chrName);

        std::vector<std::string> files;
        std::string thisFoler = folder + chrName;
        getMafftResultListFromFolder(thisFoler, files);

        std::vector<MsaFileRecord> msaFileRecords;

        std::atomic_int number_of_runing_threads(0);
        for(std::vector<std::string>::iterator itf=files.begin();
            itf!=files.end(); ++ itf ){
            bool isThisThreadUnrun = true;
            while (isThisThreadUnrun) {
                if (number_of_runing_threads < maxThread) {
                    std::thread t(msaFileReadPversion, (*itf), std::ref(thisFoler), std::ref(sdiFiles), std::ref(msaFileRecords), std::ref(number_of_runing_threads));
                    ++number_of_runing_threads;
                    t.detach();
                    isThisThreadUnrun = false;
                    break;
                }else {
                    usleep(1);
                }
            }
        }
        while (number_of_runing_threads > 0) {// wait for all the thread
            usleep(1);
        }
        std::sort(msaFileRecords.begin(), msaFileRecords.end());
        // read files end
        std::cout << "get file list done" << std::endl;


        std::atomic_int number_of_runing_threads1(0);
        for( std::map<std::string, std::string>::iterator itName=sdiFiles.begin(); itName!=sdiFiles.end(); ++itName ){
            bool isThisThreadUnrun = true;
            while (isThisThreadUnrun) {
                if (number_of_runing_threads1 < maxThread) {
                    std::thread t(newSdiFileForOneAccession, itName->first, std::ref(sdiFiles), std::ref(vcfFix), std::ref(referenceGenome), std::ref(parameters), std::ref(chrName), std::ref(msaFileRecords), std::ref(outputFolder), std::ref(number_of_runing_threads1));
                    ++number_of_runing_threads1;
                    t.detach();
                    isThisThreadUnrun = false;
                    break;
                }else {
                    usleep(10);
                }
            }
        }
        while (number_of_runing_threads1 > 0) {// wait for all the thread
            usleep(100);
        }
    }
}
