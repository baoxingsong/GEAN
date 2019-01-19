//
// Created by baoxing on 10/10/17.
//

#include "reAnnotationAndExonerateAndNovo.h"

void reAnnotationAndExonerateAndNovo( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::string novoGffFilePath, std::string sdiFile,
                                      std::map<std::string, std::vector<Gene> >& genes, std::map<std::string, Transcript > & targetTranscriptsHashMap, int & maxThread, std::string & outputGffFile, int & lengthThread, std::string & vcfFix,
                                      std::map<std::string, std::string>& parameters, int& minIntron, bool & remove_reference_orf_shift){

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    std::map<std::string, Fasta> referenceGenome;
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    readFastaFile(referenceGenomeFilePath, referenceGenome);
    std::cout << "reference genome sequence reading done" << std::endl;
    std::string regex = get_parameters("cdsParentRegex", parameters);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);
    CheckAndUpdateTranscriptsEnds( referenceTranscriptHashSet, referenceGenome, nucleotideCodeSubstitutionMatrix, minIntron);

    std::cout << "reference genome annotation reading done" << std::endl;
    std::set<std::string> toRemoveChromosomes;

    // clean data begin
    for( std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
            it!=referenceTranscriptHashSet.end(); ++it){
        if( referenceGenome.find(it->first)==referenceGenome.end() ){
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
    if( remove_reference_orf_shift ){
        std::cout << "remove ORF shift allele begin" << std::endl;
        // remove ORF shift allele from reference dataset begin
        for ( std::map<std::string, std::vector<Transcript> >::iterator it1=referenceTranscriptHashSet.begin();
              it1!=referenceTranscriptHashSet.end(); ++it1){

            //check ORF begin
            std::atomic_int number_of_runing_threads0(0);
            for (int index=0; index < (it1->second.size()); ++index){
                bool isThisThreadUnrun = true;
                while (isThisThreadUnrun) {
                    if (number_of_runing_threads0 < maxThread) {
                        std::thread t(checkOrfPversion, std::ref((it1->second)[index]), std::ref(referenceGenome), std::ref(nucleotideCodeSubstitutionMatrix), std::ref(minIntron), std::ref(number_of_runing_threads0));
                        ++number_of_runing_threads0;
                        t.detach();
                        isThisThreadUnrun = false;
                        break;
                    }else {
                        usleep(10);
                    }
                }
            }
            while (number_of_runing_threads0 > 0) {// wait for all the thread
                usleep(10);
            }//check ORF end

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
    }

    std::map<std::string, std::vector<Variant> > variantsMaps;
    readSdiFile (sdiFile, variantsMaps, vcfFix, referenceGenome);
    std::map<std::string, Fasta> targetGenome;
    getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome, parameters);
    //std::map<std::string, Transcript> targetTranscriptsHashMap;// transcriptName Transcript
    annotationLiftOver(referenceTranscriptHashSet, targetTranscriptsHashMap, variantsMaps, targetGenome, referenceGenome, parameters, minIntron);
    std::cout << "annotationLiftOver done. maxThread: " << maxThread << std::endl;

    std::vector<std::thread> threads;

    std::string prefixUuid = outputGffFile;
    std::regex reg(get_parameters("temp_file_regex", parameters));
    std::smatch match;
    if( regex_search(outputGffFile, match, reg) ){
        prefixUuid=match[3];
    }

    std::string tempFolder = createdTempFolder(parameters); // create tempFolder
    std::atomic_int number_of_runing_threads(0);
    for(std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
            it!=referenceTranscriptHashSet.end(); ++it){
        for( std::vector<Transcript>::iterator it2=it->second.begin(); it2!=it->second.end();++it2 ){
            Transcript * it3 = & targetTranscriptsHashMap[it2->getName()];
            if( (*it3).getIfOrfShift() ){
                TranscriptUpdateCdsInformation((*it2), referenceGenome);
                bool isThisThreadUnrun = true;
                while (isThisThreadUnrun) {
                    if (number_of_runing_threads < maxThread) {
                        std::thread t(transcriptRealignmentAndExonerate, *it3, *it2, std::ref(nucleotideCodeSubstitutionMatrix),  std::ref(targetGenome),
                                      std::ref(referenceGenome),
                                      it->first, std::ref(targetTranscriptsHashMap), std::ref(number_of_runing_threads), std::ref(prefixUuid), std::ref(lengthThread), std::ref(parameters), std::ref(minIntron));
                        ++number_of_runing_threads;
                        t.detach();
                        isThisThreadUnrun = false;
                        break;
                    } else {
                        usleep(1000);
                    }
                }
            }
        }
    }
    while( number_of_runing_threads >0 ){// wait for all the thread
        usleep(1000);
    }
    std::cerr << "realignment done" << std::endl;

    //replace old orf-shift annotation with novo annotation begin
    std::map<std::string, std::vector<Transcript> > novoTranscriptHashSet;
    std::string novoRegex = get_parameters("novo_cdsParentRegex", parameters);
    readGffFile (novoGffFilePath, novoTranscriptHashSet, novoRegex);
    CheckAndUpdateTranscriptsEnds( novoTranscriptHashSet, referenceGenome, nucleotideCodeSubstitutionMatrix, minIntron);
    std::map<std::string, Transcript> adaptedNovoTranscriptAnnotation1;
    std::map<std::string, Transcript> adaptedNovoTranscriptAnnotation2;
    for( std::map<std::string, std::vector<Transcript> >::iterator it=novoTranscriptHashSet.begin(); it!=novoTranscriptHashSet.end(); ++it){
        if( referenceGenome.find(it->first) != referenceGenome.end() ){
            for(std::vector<Transcript>::iterator it2=novoTranscriptHashSet[it->first].begin();
                it2!=novoTranscriptHashSet[it->first].end(); ++it2){
                bool ifHasOverLap=false;
                for( std::map<std::string, Transcript>::iterator it3=targetTranscriptsHashMap.begin();
                     it3!=targetTranscriptsHashMap.end(); ++it3){
                    if( it2->ifOverLap(it3->second) ){
                        ifHasOverLap=true;
                        if( it3->second.getIfOrfShift() ){
                            TranscriptUpdateCdsInformation((*it2), targetGenome);
                            checkOrfState( (*it2), targetGenome, nucleotideCodeSubstitutionMatrix, minIntron);
                            if( it2->getIfOrfShift() ){
                                std::cout << "There are de novo gene structure with ORF incomplete\n" << it2->getName() << " " << it2->getMetaInformation()  << std::endl;
                            }else{
                                it2->setSource("DENOVOOVERLAP");
                                it2->setName(it3->second.getName());
                                adaptedNovoTranscriptAnnotation1[it2->getName()]=(*it2);
                            }
                        }
                    }
                }
                if( !ifHasOverLap ){
                    TranscriptUpdateCdsInformation((*it2), targetGenome);
                    //it2->updateInfor(targetGenome);
                    checkOrfState( (*it2), targetGenome, nucleotideCodeSubstitutionMatrix, minIntron);
                    if( it2->getIfOrfShift() ){
                        std::cout << "There are de novo gene structure with ORF incomplete\n" << it2->getName() << " " << it2->getMetaInformation()  << std::endl;
                    }else {
                        it2->setSource("DENOVO");
                        adaptedNovoTranscriptAnnotation2[it2->getName()] = (*it2);
                    }
                }
            }
        }
    }//replace old orf-shift annotation with novo annotation end

    std::cerr << "de novo annotation included" << std::endl;
    // remove those ORF-shift annotation begin
    std::vector<std::string> transcriptToRemove;
    for( std::map<std::string, Transcript>::iterator it=targetTranscriptsHashMap.begin();
         it!=targetTranscriptsHashMap.end(); ++it){
        if( it->second.getIfOrfShift() ){
            transcriptToRemove.push_back(it->first);
        }else{
            std::cerr << "lift over " <<it->first << " is good"  << std::endl;
        }
    }

    for( std::vector<std::string>::iterator it=transcriptToRemove.begin();
         it!=transcriptToRemove.end(); ++it){
        if( targetTranscriptsHashMap.find(*it) != targetTranscriptsHashMap.end() ){
            targetTranscriptsHashMap.erase(*it);
        }
    }

    for(  std::map<std::string, Transcript>::iterator it=adaptedNovoTranscriptAnnotation1.begin();
          it!=adaptedNovoTranscriptAnnotation1.end(); ++it){
        targetTranscriptsHashMap[it->first]=it->second;
    }
    std::string transcript_to_gene_regex_reference_gff = get_parameters("transcript_to_gene_regex_reference_gff", parameters);
    std::string transcript_to_gene_regex_novo_gff = get_parameters("transcript_to_gene_regex_novo_gff", parameters);
    // remove those ORF-shift annotation end
    std::cerr << "ORF shift transcripts removed" << std::endl;
    std::map<std::string, std::string > transcript_to_gene_map;
    get_transcript_to_gene_map_from_gff (referenceGffFilePath, transcript_to_gene_map, transcript_to_gene_regex_reference_gff);
    TranscriptsTogenes(transcript_to_gene_map, genes, targetTranscriptsHashMap);

    transcript_to_gene_map.empty();
    get_transcript_to_gene_map_from_gff (novoGffFilePath, transcript_to_gene_map, transcript_to_gene_regex_novo_gff);
    std::map<std::string, std::vector<Gene>> genes2;
    TranscriptsTogenes(transcript_to_gene_map, genes2, adaptedNovoTranscriptAnnotation2);

    for(  std::map<std::string, Transcript>::iterator it=adaptedNovoTranscriptAnnotation2.begin();
          it!=adaptedNovoTranscriptAnnotation2.end(); ++it){
        targetTranscriptsHashMap[it->first]=it->second;
    }

    for( std::map<std::string, std::vector<Gene>>::iterator it=genes2.begin();
            it!=genes2.end(); ++it){
        if( genes.find(it->first) == genes.end() ){
            genes[it->first]=std::vector<Gene>();
        }
        for( std::vector<Gene>::iterator it2=genes2[it->first].begin();
             it2!=genes2[it->first].end(); ++it2 ){
            genes[it->first].push_back(*it2);
        }
    }
    std::cerr << "transcript structure to gene structure done" << std::endl;
}
