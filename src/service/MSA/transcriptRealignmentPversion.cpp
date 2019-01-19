//
// Created by baoxing on 10/10/17.
//

#include "transcriptRealignmentPversion.h"

std::mutex gmutexTranscriptRealignmentPversion;

void transcriptRealignmentPversion( Transcript& tartgetTranscript, Transcript& referenceTranscript,
                                    NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                    std::map<std::string, Fasta>& targetGenome, std::atomic_int & number_of_runing_threads,
                                    std::map<std::string, Fasta>& referenceGenome, std::string chromosomeName, std::map<std::string,
                                    Transcript>& targetTranscriptsHashMap, int & lengthThread,
                                    std::map<std::string, std::string>& parameters, int & minIntron ){
    std::string refGenomeSequence = referenceTranscript.getGeneomeSequence();
    Transcript targetTranscript(tartgetTranscript.getName(), chromosomeName, tartgetTranscript.getStrand());
    int startTarget = tartgetTranscript.getPStart()-refGenomeSequence.length();
    int endTarget = tartgetTranscript.getPEnd()+refGenomeSequence.length();

    if( endTarget > targetGenome[chromosomeName].getSequence().length() ){
        endTarget = targetGenome[chromosomeName].getSequence().length();
    }

    if( transcriptRealignment( targetTranscript, startTarget,  endTarget, referenceTranscript,
                               nucleotideCodeSubstitutionMatrix,
                           targetGenome, refGenomeSequence,
                           referenceGenome, chromosomeName,
                           targetTranscriptsHashMap, lengthThread,
                           parameters, minIntron )){

        if( targetTranscript.getIfOrfShift() ){

        }else{
            gmutexTranscriptRealignmentPversion.lock();
            targetTranscriptsHashMap[targetTranscript.getName()]=targetTranscript; // update the map data structure with modified transcript
            gmutexTranscriptRealignmentPversion.unlock();
        }
    }
    --number_of_runing_threads;
}
