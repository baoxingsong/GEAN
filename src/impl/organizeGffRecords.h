//
// Created by Baoxing song on 20.10.18.
//

#ifndef ZSDP_ORGANIZEGFFRECORDS_H
#define ZSDP_ORGANIZEGFFRECORDS_H

#include "../model/model.h"
#include "../util/util.h"
#include "./TranscriptUpdateInformation.h"
#include "./checkOrfState.h"

bool exactlySame( Gene & g1, Gene & g2);
void mergeToFirstOne( Gene & g1, Gene & g2);
void updateGeneInformation(std::map<std::string, std::vector<Gene> > & geneMap, NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix,
                           const size_t & minIntron,  std::map<std::string, Fasta> & querySequences);
void removeDuplication(std::map<std::string, std::vector<Gene> > & geneMap, const size_t & minGene, std::map<std::string, std::set<int32_t>> & toRemove);
#endif //ZSDP_ORGANIZEGFFRECORDS_H
