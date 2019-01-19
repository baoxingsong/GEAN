//
// Created by baoxing on 10/10/17.
//

#include "TranscriptsTogenes.h"


void TranscriptsTogenes(std::map<std::string, std::string >& transcript_to_gene_map, std::map<std::string, std::vector<Gene> >& genes,
                        std::map<std::string, Transcript>& targetTranscriptsHashMap){
    std::map< std::string, std::map<std::string, Gene > > tempGenes;
    for( std::map<std::string, Transcript>::iterator it=targetTranscriptsHashMap.begin();
         it!=targetTranscriptsHashMap.end(); ++it){
        std::string transcriptId = it->first;
        std::string chromosome = it->second.getChromeSomeName();
        if(tempGenes.find(chromosome) == tempGenes.end() ){
            tempGenes[chromosome]=std::map<std::string, Gene >();
        }
        if( transcript_to_gene_map.find(transcriptId)!=transcript_to_gene_map.end() ){
            std::string geneName = transcript_to_gene_map[transcriptId];
            if( tempGenes[chromosome].find(geneName) == tempGenes[chromosome].end() ){
                STRAND strand = it->second.getStrand();
                Gene gene(geneName, strand);
                tempGenes[chromosome].insert( std::pair<std::string, Gene>(geneName, gene) );
            }
            tempGenes[chromosome][geneName].addTranscript(it->second);
            std::string source = it->second.getSource();
            tempGenes[chromosome][geneName].setSource(source);
        }else{
            std::cerr << "could not map of transcript ID " << transcriptId << " to gene ID" << std::endl;
            exit(1);
        }
    }
    for( std::map< std::string, std::map<std::string, Gene > >::iterator it=tempGenes.begin();
            it!=tempGenes.end(); ++it ){
        if( genes.find(it->first) == genes.end()  ){
            genes[it->first]=std::vector<Gene>();
        }
        for( std::map<std::string, Gene >::iterator it2=it->second.begin();
             it2!=it->second.end(); ++it2){

            int start = INT_MAX;
            int end=0;
            for( std::string t : it2->second.getTranscriptVector() ){
                if( targetTranscriptsHashMap[t].getPEnd() > end ){
                    end = targetTranscriptsHashMap[t].getPEnd();
                }
                if(  targetTranscriptsHashMap[t].getPStart() < start ){
                    start = targetTranscriptsHashMap[t].getPStart();
                }
            }
            it2->second.setStart(start);
            it2->second.setEnd(end);
            genes[it->first].push_back(it2->second);
        }
    }
}
