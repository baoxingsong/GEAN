//
// Created by baoxing on 10/10/17.
//

#include "GetPseudoGenomeSequence.h"


int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, const std::string& sdiFile,
                            std::map<std::string, Fasta>& targetSequences, const std::string& vcfFix, std::map<std::string, std::string>& parameters){
    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix, referenceSequences);
    return getPseudoGenomeSequence(referenceSequences, sdiMaps, targetSequences, parameters);
}
int getPseudoGenomeSequence(const std::string& referenceGenomeFastaFile,
                            const std::string& sdiFile, std::map<std::string, Fasta>& targetSequences, const std::string& vcfFix, std::map<std::string, std::string>& parameters){
    std::map<std::string, Fasta> referenceSequences;
    readFastaFile(referenceGenomeFastaFile, referenceSequences);
    return getPseudoGenomeSequence(referenceSequences, sdiFile, targetSequences, vcfFix, parameters);
}
int getPseudoGenomeSequence(const std::string& referenceGenomeFastaFile,
                            const std::string& sdiFile, const std::string& outputFile, const std::string& vcfFix, std::map<std::string, std::string>& parameters){
    std::map<std::string, Fasta> targetSequences;
    int resultCode = getPseudoGenomeSequence(referenceGenomeFastaFile, sdiFile, targetSequences, vcfFix, parameters);
    std::ofstream ofile;
    ofile.open(outputFile);
    for( std::map<std::string, Fasta>::iterator i=targetSequences.begin() ; i!= targetSequences.end(); i ++ ){
        std::string name = i->first;
        std::string sequence = i->second.getSequence();
        int lineWidth=60;
        writeFasta(ofile, name, sequence,  lineWidth);
    }
    ofile.close();
    return resultCode;
}


int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences, std::map<std::string, std::string>& parameters) {
    std::set<std::string> chromosomes;
    for (std::map<std::string, Fasta>::iterator iterFasta = referenceSequences.begin();
         iterFasta != referenceSequences.end(); iterFasta++) {
        chromosomes.insert(iterFasta->first);
    }
    return getPseudoGenomeSequence(referenceSequences, variantsMap, targetSequences, chromosomes, parameters);
}

int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences, std::string& chromosome, std::map<std::string, std::string>& parameters){
    std::set<std::string> chromosomes;
    chromosomes.insert(chromosome);
    return getPseudoGenomeSequence(referenceSequences, variantsMap, targetSequences, chromosomes, parameters);
}

int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences, std::set<std::string>& chromosomes, std::map<std::string, std::string>& parameters){
    for( std::map<std::string, Fasta>::iterator iterFasta=referenceSequences.begin(); iterFasta!=referenceSequences.end(); iterFasta++ ) {
        std::string chromosome = iterFasta->first;
        if( chromosomes.find(chromosome)!=chromosomes.end() ){
            std::string referenceSequence = iterFasta->second.getSequence();
            int totalSize = referenceSequence.size();
            std::string reserveString;
            reserveString.reserve(totalSize * 1.5);
            std::stringstream sequencestream(reserveString);
            int currentPosition = 1;
            for (std::vector<Variant>::iterator iterSdi = variantsMap[chromosome].begin(); iterSdi != variantsMap[chromosome].end(); iterSdi++) {
                int position = (*iterSdi).getPosition();
                std::string ref = (*iterSdi).getReference();
                std::string alter = (*iterSdi).getAlternative();
                int changingLength = (*iterSdi).getChanginglength();
                if (currentPosition < position) {
                    sequencestream << referenceSequence.substr(currentPosition - 1, position - currentPosition);
                } else if( currentPosition == position+1 && changingLength>0 && ref.compare("-")==0 && (*(iterSdi-1)).getChanginglength()==0 && (*(iterSdi-1)).getReference().size()==1 ){
                    currentPosition = position+1;
                    sequencestream.seekp(-1,sequencestream.cur);
                    sequencestream << alter << (*(iterSdi-1)).getAlternative();
                    continue;
                } else if( currentPosition == position+1 && changingLength>0 && ref.compare("-")==0 && (*(iterSdi-1)).getChanginglength() < 0 && (*(iterSdi-1)).getAlternative().compare("-")==0 ){
                    currentPosition = position+1;
                    sequencestream.seekp(-1,sequencestream.cur);
                    sequencestream << alter << (*(iterSdi-1)).getAlternative();
                    continue;
                } else if (currentPosition > position) {
                    std::cerr << "the sdi file is not well sorted, it should be sorted with coordinate " <<std::endl
                              << chromosome << ": currentPosition:"<< currentPosition << " position:" << position << std::endl;
                    exit(1);
                }
                if (changingLength == 0) {
                    currentPosition = position + ref.size();
                    sequencestream << alter;
                }
                if (changingLength > 0) {
                    if ( ref.compare("-")==0 ) {
                        sequencestream << alter;
                        currentPosition = position;
                    } else {
                        sequencestream << alter;
                        currentPosition = position + ref.size();
                    }
                }
                if (changingLength < 0) {
                    if (alter.compare("-")==0 ) {
                        currentPosition = position - changingLength;
                    } else {
                        sequencestream << alter;
                        currentPosition = position + ref.size();
                    }
                }
            }

            sequencestream << referenceSequence.substr(currentPosition - 1, totalSize - currentPosition + 1);
            std::string targetSequence = sequencestream.str();
            Fasta fasta(chromosome, targetSequence);
            targetSequences[chromosome]=fasta;
        }
    }
    return 0;
}
