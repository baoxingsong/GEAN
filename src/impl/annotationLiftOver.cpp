//
// Created by baoxing on 10/10/17.
//

#include "annotationLiftOver.h"


void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                        std::string chromosome, std::map<std::string, std::string>& parameters, const int& minIntron){
    std::set<std::string> chromosomes;
    chromosomes.insert(chromosome);
    annotationLiftOver( refTranscriptHashSet,
                        targetTranscriptHashMap,
                        variantsMap,
                        targetGenome, referenceGenome,
                        chromosomes, parameters, minIntron);
}

void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                        std::map<std::string, std::string>& parameters, const int& minIntron){
    std::set<std::string> chromosomes;
    for(std::map<std::string,std::vector<Transcript>>::iterator it1=refTranscriptHashSet.begin(); it1!=refTranscriptHashSet.end(); it1++){
        std::string chromosome = it1->first;
        chromosomes.insert(chromosome);
    }
    annotationLiftOver(refTranscriptHashSet, targetTranscriptHashMap, variantsMap, targetGenome, referenceGenome, chromosomes, parameters, minIntron);
}

void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                        std::set<std::string>& chromosomes, std::map<std::string, std::string>& parameters, const int& minIntron){
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    for(std::map<std::string,std::vector<Transcript> >::iterator it1=refTranscriptHashSet.begin(); it1!=refTranscriptHashSet.end(); ++it1){
        std::string chromosome = it1->first;
        if( chromosomes.find(chromosome)!=chromosomes.end() ){
            for( std::vector<Transcript>::iterator it2=it1->second.begin(); it2!=it1->second.end(); ++it2){
                //Transcript referenceTranscript = (*it2);
                Transcript targetTranscript(it2->getName(), chromosome, it2->getStrand());
                for( std::vector<GenomeBasicFeature>::iterator it3=it2->getCdsVector().begin(); it3!=it2->getCdsVector().end(); ++it3 ){
                    int liftStart = getChangedFromBasement(chromosome, it3->getStart(), variantsMap);
                    int liftEnd = getChangedFromBasement(chromosome, it3->getEnd(), variantsMap);
                    GenomeBasicFeature cds(liftStart, liftEnd);
                    targetTranscript.addCds(cds);
                }
                targetTranscript.setSource("LIFTOVER");
                TranscriptUpdateCdsInformation(targetTranscript, targetGenome);
                checkOrfState(targetTranscript, targetGenome, nucleotideCodeSubstitutionMatrix, minIntron);
                targetTranscriptHashMap[targetTranscript.getName()]=targetTranscript;
            }
        }
    }
}

void annotationLiftOverPversion(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                        std::map<std::string, std::string>& parameters, int& minIntron, const int & maxThread){

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    for(std::map<std::string,std::vector<Transcript> >::iterator it1=refTranscriptHashSet.begin(); it1!=refTranscriptHashSet.end(); it1++){
        std::string chromosome = it1->first;
        for( std::vector<Transcript>::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++){
            Transcript referenceTranscript = (*it2);
            Transcript targetTranscript((*it2).getName(), chromosome, (*it2).getStrand());
            for( std::vector<GenomeBasicFeature>::iterator it3=referenceTranscript.getCdsVector().begin(); it3!=referenceTranscript.getCdsVector().end();it3++ ){
                int liftStart = getChangedFromBasement(chromosome, (*it3).getStart(), variantsMap);
                int liftEnd = getChangedFromBasement(chromosome, (*it3).getEnd(), variantsMap);
                GenomeBasicFeature cds(liftStart, liftEnd);
                targetTranscript.addCds(cds);
            }
            targetTranscript.setSource("LIFTOVER");
            //TranscriptUpdateInformation(targetTranscript, targetGenome);
            targetTranscriptHashMap[targetTranscript.getName()]=targetTranscript;
        }
    }

    std::cout << "liftover done, begin to check ORF state" << std::endl;

//    struct timespec ts, ts1;
//    ts.tv_nsec = 150000;    // 15ms
//    ts.tv_sec = 1;

    std::atomic_int number_of_runing_threads2(0);
    for( std::map<std::string, Transcript>::iterator it1=targetTranscriptHashMap.begin();
            it1!=targetTranscriptHashMap.end(); ++it1){
        bool isThisThreadUnrun = true;
        while (isThisThreadUnrun) {
            if (number_of_runing_threads2 < maxThread) {
                //std::thread t(checkOrfPversion, std::ref((it1->second)[index]), std::ref(referenceGenome), std::ref(nucleotideCodeSubstitutionMatrix), std::ref(minIntron), std::ref(number_of_runing_threads0));

                std::thread t2(checkOrfPversion, std::ref(it1->second), std::ref(targetGenome), std::ref(nucleotideCodeSubstitutionMatrix), std::ref(minIntron), std::ref(number_of_runing_threads2) );

                ++number_of_runing_threads2;
                t2.detach();
                isThisThreadUnrun = false;
                break;
            } else {
                usleep(10);
            }
        }
    }
    while (number_of_runing_threads2 > 0) {// wait for all the thread
        usleep(100);
    }
}
