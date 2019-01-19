//
// Created by baoxing on 10/10/17.
//

#include "transcriptRealignmentAndExonerate.h"


std::mutex gmutexTranscriptRealignmentAndExonerate;


bool transcriptRealignment( Transcript& targetTranscript, int& startTarget, int & endTarget, Transcript& referenceTranscript,
                                    NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                    std::map<std::string, Fasta>& targetGenome, std::string& refGenomeSequence,
                                    std::map<std::string, Fasta>& referenceGenome, std::string& chromosomeName,
                                    std::map<std::string, Transcript>& targetTranscriptsHashMap, int & lengthThread,
                                    std::map<std::string, std::string>& parameters, int & minIntron ){

    int startCodonPosition=1;
    int stopCodonPosition=refGenomeSequence.length()-2;
    std::vector<SpliceSitePosition> spliceSitePositions;

    targetTranscript.setSource("REALIGNMENT");
    int targetPosition=0;
    int referencePosition=0;
    //prepare for the new GenomeBasicFeature After re-alignment begin
    std::vector<int> targetCdsStarts;
    std::vector<int> targetCdsEnds;
    STRAND thisStrand = referenceTranscript.getStrand();
    std::string dna_b = getSubsequence(targetGenome, chromosomeName, startTarget, endTarget, thisStrand);

    if( dna_b.length()<=lengthThread && refGenomeSequence.length()<=lengthThread ){ //  if the sequence is too long, don't try to align it
        //  if the sequence is too long, don't try to align it
        if( referenceTranscript.getStrand() == POSITIVE ){
            if( referenceTranscript.getCdsVector().size()>1 ){
                for( size_t i=1; i<referenceTranscript.getCdsVector().size(); ++i ){
                    SpliceSitePosition spliceSitePosition(
                            referenceTranscript.getCdsVector()[i-1].getEnd()-referenceTranscript.getPStart()+2,
                            referenceTranscript.getCdsVector()[i].getStart()-referenceTranscript.getPStart());
                    spliceSitePositions.push_back(spliceSitePosition);
                }
            }
            alignNeedlemanForTranscript_simd_avx2int32 nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, spliceSitePositions, parameters, nucleotideCodeSubstitutionMatrix);
            /*if(referenceTranscript.getCdsVector().size()>1){
                nw.print_results();
            }*/
            for( size_t tp = 0; tp<nw.getAlignment_q().length(); ++tp ){
                if( nw.getAlignment_q()[tp] != '-' ){
                    ++targetPosition;
                }
                if( nw.getAlignment_d()[tp] != '-' ){
                    ++referencePosition;
                    for( std::vector<GenomeBasicFeature>::iterator it4=referenceTranscript.getCdsVector().begin();
                         it4!=referenceTranscript.getCdsVector().end(); ++it4){
                        if( referencePosition+referenceTranscript.getPStart()-1 == (*it4).getStart() ){
                            targetCdsStarts.push_back(startTarget + targetPosition-1); //(*i3) is the target transcript, refGenomeSequence.length() is the extend length
                        }
                        if( referencePosition+referenceTranscript.getPStart()-1 == (*it4).getEnd() ){ //todo recheck here and compare the code with the TransferGffWithNucmerResult
                            targetCdsEnds.push_back(startTarget + targetPosition-1);
                        }
                    }
                }
            }
        } else {
            if( referenceTranscript.getCdsVector().size()>1 ){
                for( size_t i=referenceTranscript.getCdsVector().size()-1; i>0; --i ){
                    SpliceSitePosition spliceSitePosition(
                            referenceTranscript.getPEnd()-referenceTranscript.getCdsVector()[i].getStart()+2,
                            referenceTranscript.getPEnd()-referenceTranscript.getCdsVector()[i-1].getEnd());
                    spliceSitePositions.push_back(spliceSitePosition);
                }
            }
            NeedlemanWunschForTranscript nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, spliceSitePositions, parameters, nucleotideCodeSubstitutionMatrix);
            for( size_t tp = 0; tp<nw.getAlignment_d().length(); ++tp ){
                if( nw.getAlignment_q()[tp] != '-' ){
                    ++targetPosition;
                }
                if( nw.getAlignment_d()[tp] != '-' ){
                    ++referencePosition;
                    for( std::vector<GenomeBasicFeature>::iterator it4=referenceTranscript.getCdsVector().begin();
                         it4!=referenceTranscript.getCdsVector().end(); it4++){
                        if( - referencePosition+referenceTranscript.getPEnd()+1 == (*it4).getStart() ){
                            targetCdsStarts.push_back(endTarget - targetPosition +1);
                        }
                        if( - referencePosition+referenceTranscript.getPEnd()+1 == (*it4).getEnd() ){
                            targetCdsEnds.push_back(endTarget - targetPosition +1);
                        }
                    }
                }
            }
        }
        for(size_t i5=0; i5<targetCdsStarts.size(); i5++ ){
            GenomeBasicFeature cds(targetCdsStarts[i5], targetCdsEnds[i5]);
            //std::cout << "CDS " << targetCdsStarts[i5] << " " << targetCdsEnds[i5] << std::endl;
            targetTranscript.addCds(cds);
        }
        TranscriptUpdateCdsInformation(targetTranscript, targetGenome);
        checkOrfState( targetTranscript, targetGenome, nucleotideCodeSubstitutionMatrix, minIntron);

        if( targetTranscript.getIfOrfShift() ){
            return false;
        }else{
            return true;
        }
    }
    return false;
}



void transcriptRealignmentAndExonerate( Transcript tartgetTranscript, Transcript referenceTranscript,
                                        NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                        std::map<std::string, Fasta>& targetGenome,
                                        std::map<std::string, Fasta>& referenceGenome, std::string chromosomeName, std::map<std::string,
                                        Transcript>& targetTranscriptsHashMap, std::atomic_int & number_of_runing_threads, std::string & prefixUuid, int & lengthThread,
                                        std::map<std::string, std::string>& parameters, int& minIntron ){
    std::cout << "realigning " << referenceTranscript.getName() << " begin" << std::endl;
    std::string refGenomeSequence = referenceTranscript.getGeneomeSequence();

    Transcript targetTranscript(tartgetTranscript.getName(), chromosomeName, tartgetTranscript.getStrand());
    int startTarget = tartgetTranscript.getPStart()-refGenomeSequence.length();
    int endTarget = tartgetTranscript.getPEnd()+refGenomeSequence.length();

    if( endTarget > targetGenome[chromosomeName].getSequence().length() ){
        endTarget = targetGenome[chromosomeName].getSequence().length();
    }

    if( transcriptRealignment( targetTranscript, startTarget, endTarget, referenceTranscript,
                                nucleotideCodeSubstitutionMatrix,
                                targetGenome, refGenomeSequence,
                                referenceGenome, chromosomeName,
                                targetTranscriptsHashMap, lengthThread,
                                parameters, minIntron ) ){

        if( targetTranscript.getIfOrfShift() ){
            std::string cdsSequence=referenceTranscript.getCdsSequence();
            if(cdsSequence.length() >0 && cdsSequence.length() < lengthThread*4){
                std::string transcriptName = referenceTranscript.getName();
                std::string chrName = referenceTranscript.getChromeSomeName();
                std::string targetSequence = getSubsequence(targetGenome, chrName,
                                                            startTarget, endTarget);

                runExonerateEst(transcriptName, cdsSequence, targetSequence,
                                nucleotideCodeSubstitutionMatrix, targetTranscriptsHashMap, startTarget, endTarget,
                                referenceTranscript.getStrand(), chromosomeName, prefixUuid, targetGenome, parameters, minIntron);
            }
        }else{
            gmutexTranscriptRealignmentAndExonerate.lock();
            targetTranscriptsHashMap[targetTranscript.getName()]=targetTranscript; // update the map data structure with modified transcript
            gmutexTranscriptRealignmentAndExonerate.unlock();
        }
    } else {
        std::string cdsSequence=referenceTranscript.getCdsSequence();
        if(cdsSequence.length() >0 && cdsSequence.length() < lengthThread*4 ){
            std::string transcriptName = referenceTranscript.getName();
            std::string chrName = referenceTranscript.getChromeSomeName();
            std::string targetSequence = getSubsequence(targetGenome, chrName,
                                                        startTarget, endTarget);

            runExonerateEst(transcriptName, cdsSequence, targetSequence,
                            nucleotideCodeSubstitutionMatrix, targetTranscriptsHashMap, startTarget, endTarget,
                            referenceTranscript.getStrand(), chromosomeName, prefixUuid, targetGenome, parameters, minIntron);

        }
    }
    --number_of_runing_threads;
    std::cout << "realigning " << referenceTranscript.getName() << " done" << std::endl;
}


//
//    int startCodonPosition=1;
//    int stopCodonPosition=refGenomeSequence.length()-2;
//    std::vector<SpliceSitePosition> spliceSitePositions;
//
//    targetTranscript.setSource("REALIGNMENT");
//    int targetPosition=0;
//    int referencePosition=0;
//    //prepare for the new GenomeBasicFeature After re-alignment begin
//    std::vector<int> targetCdsStarts;
//    std::vector<int> targetCdsEnds;
//    STRAND thisStrand =  referenceTranscript.getStrand();
//    std::string dna_b = getSubsequence(targetGenome, chromosomeName, startTarget, endTarget, thisStrand);
//
//    if( dna_b.length()<=lengthThread && refGenomeSequence.length()<=lengthThread ){ //  if the sequence is too long, don't try to align it
//        if( referenceTranscript.getStrand() == POSITIVE ){
//            if( referenceTranscript.getCdsVector().size()>1 ){
//                for( size_t i=1; i<referenceTranscript.getCdsVector().size(); ++i ){
//                    SpliceSitePosition spliceSitePosition(
//                            referenceTranscript.getCdsVector()[i-1].getEnd()-referenceTranscript.getStart()+2,
//                            referenceTranscript.getCdsVector()[i].getStart()-referenceTranscript.getStart());
//                    spliceSitePositions.push_back(spliceSitePosition);
//                }
//            }
//            NeedlemanWunschForTranscript nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, spliceSitePositions, parameters);
//            for( size_t tp = 0; tp<nw.getAlignment_a().length(); ++tp ){
//                if( nw.getAlignment_b()[tp] != '-' ){
//                    ++targetPosition;
//                }
//                if( nw.getAlignment_a()[tp] != '-' ){
//                    ++referencePosition;
//                    for( std::vector<GenomeBasicFeature>::iterator it4=referenceTranscript.getCdsVector().begin();
//                         it4!=referenceTranscript.getCdsVector().end(); ++it4){
//                        if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getStart() ){
//                            //if( nw.getAlignment_b()[tp] == '-' ){
//                            targetCdsStarts.push_back(startTarget + targetPosition-1); //(*i3) is the target transcript, refGenomeSequence.length() is the extend length
////                            }else{
////                                targetCdsStarts.push_back(startTarget + targetPosition -1);
////                            }
//                        }
//                        if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getEnd() ){
//                            //if( nw.getAlignment_b()[tp] == '-' ){
//                            targetCdsEnds.push_back(startTarget + targetPosition-1);
////                            }else{
////                                targetCdsEnds.push_back(startTarget + targetPosition -1);
////                            }
//                        }
//                    }
//                }
//            }
//
//
//        } else {
//            if( referenceTranscript.getCdsVector().size()>1 ){
//                for( size_t i=referenceTranscript.getCdsVector().size()-1; i>0; --i ){
//                    SpliceSitePosition spliceSitePosition(
//                            referenceTranscript.getEnd()-referenceTranscript.getCdsVector()[i].getStart()+2,
//                            referenceTranscript.getEnd()-referenceTranscript.getCdsVector()[i-1].getEnd()-1);
//                    spliceSitePositions.push_back(spliceSitePosition);
//                }
//            }
//            NeedlemanWunschForTranscript nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, spliceSitePositions, parameters);
//            for( size_t tp = 0; tp<nw.getAlignment_a().length(); ++tp ){
//                if( nw.getAlignment_b()[tp] != '-' ){
//                    ++targetPosition;
//                }
//                if( nw.getAlignment_a()[tp] != '-' ){
//                    ++referencePosition;
//                    for( std::vector<GenomeBasicFeature>::iterator it4=referenceTranscript.getCdsVector().begin();
//                         it4!=referenceTranscript.getCdsVector().end(); it4++){
//                        //if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getStart() ){
//                        if( - referencePosition+referenceTranscript.getEnd()+1 == (*it4).getStart() ){
//                            //if( nw.getAlignment_b()[tp] == '-' ){
//                            targetCdsStarts.push_back(endTarget - targetPosition +1);
////                            }else{
////                                targetCdsStarts.push_back(endTarget - targetPosition +1);
////                            }
//                        }
//                        //if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getEnd() ){
//                        if( - referencePosition+referenceTranscript.getEnd()+1 == (*it4).getEnd() ){
//                            //if( nw.getAlignment_b()[tp] == '-' ){
//                            targetCdsEnds.push_back(endTarget - targetPosition +1);
////                            }else{
////                                targetCdsEnds.push_back(endTarget - targetPosition +1);
////                            }
//                        }
//                    }
//                }
//            }
//        }
//
//
//        for(size_t i5=0; i5<targetCdsStarts.size(); i5++ ){
//            GenomeBasicFeature cds(targetCdsStarts[i5], targetCdsEnds[i5]);
//            //std::cout << "CDS " << targetCdsStarts[i5] << " " << targetCdsEnds[i5] << std::endl;
//            targetTranscript.addCds(cds);
//        }
//        TranscriptUpdateInformation(targetTranscript, targetGenome);
////        targetTranscript.updateInfor(targetGenome);
////        checkOrfState( targetTranscript, referenceTranscript,
////                       targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix);
//        checkOrfState( targetTranscript, targetGenome,  nucleotideCodeSubstitutionMatrix, minIntron);

