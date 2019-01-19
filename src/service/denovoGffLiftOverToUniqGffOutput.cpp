//
// Created by Baoxing song on 19.10.18.
//

#include "denovoGffLiftOverToUniqGffOutput.h"


void denovoGffLiftOverToUniqGffOutput(const std::string & gffFilePath, const std::string & outPutFilePath,
        const std::string & queryFastaFilePath, const size_t & minIntron, const size_t & minGene,
        std::map<std::string, std::string>& parameters){
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);

    std::ofstream ofile;
    ofile.open(outPutFilePath);
    std::map<std::string, std::vector<Gene> > geneMap;
    readGffFileWithEveryThing ( gffFilePath, geneMap);

    std::map<std::string, Fasta> querySequences;
    readFastaFile(queryFastaFilePath, querySequences);

    updateGeneInformation(geneMap, nucleotideCodeSubstitutionMatrix, minIntron, querySequences);

    std::map<std::string, std::set<int32_t>> toRemove;
    removeDuplication(geneMap, minGene, toRemove);

    for( std::map<std::string, std::vector<Gene> >::iterator it0=geneMap.begin(); it0!=geneMap.end(); ++it0 ){
        size_t index=0;
        for( std::vector<Gene>::iterator it1=it0->second.begin(); it1!=it0->second.end(); ++it1 ){
            if( toRemove[it0->first].find(index) == toRemove[it0->first].end() ){
                std::string st = "+";
                if( NEGATIVE ==  it1->getTranscripts()[0].getStrand() ){
                    st="-";
                }
                ofile << it1->getTranscripts()[0].getChromeSomeName() << "\t"+it1->getTranscripts()[0].getSource()+"\tgene\t" << it1->getStart() << "\t" <<
                      it1->getEnd() << "\t.\t"<< st <<"\t.\t"<< it1->getLastColumnInformation() << std::endl;
                for( std::vector<Transcript>::iterator it=it1->getTranscripts().begin();
                     it!=it1->getTranscripts().end(); ++it ){
                    ofile << it->getChromeSomeName() << "\t"+it->getSource()+"\t"<< it->getType() <<"\t" << it->getStart() << "\t" <<
                          it->getEnd() << "\t" << it->getScore() << "\t"<< st <<"\t.\t"<< it->getLastColumnInformation() << std::endl;
                    std::vector<GenomeBasicFeature> GenomeBasicFeatures;
                    if( it->getFivePrimerUtr().size()>0){
                        for (std::vector<GenomeBasicFeature>::iterator it4 = it->getFivePrimerUtr().begin();
                             it4 != it->getFivePrimerUtr().end(); ++it4) {
                            GenomeBasicFeatures.push_back(*it4);
                        }
                    }
                    if( it->getCdsVector().size() > 0 ) {
                        for (std::vector<GenomeBasicFeature>::iterator it4 = it->getCdsVector().begin();
                             it4 != it->getCdsVector().end(); ++it4) {
                            GenomeBasicFeatures.push_back(*it4);
                        }
                    }
                    if( it->getExonVector().size() > 0 ) {
                        for (std::vector<GenomeBasicFeature>::iterator it4 = it->getExonVector().begin();
                             it4 != it->getExonVector().end(); ++it4) {
                            GenomeBasicFeatures.push_back(*it4);
                        }
                    }
                    if( it->getThreePrimerUtr().size()>0){
                        for (std::vector<GenomeBasicFeature>::iterator it4 = it->getThreePrimerUtr().begin();
                             it4 != it->getThreePrimerUtr().end(); ++it4) {
                            GenomeBasicFeatures.push_back(*it4);
                        }
                    }
                    if( GenomeBasicFeatures.size() > 0 ) {
                        std::sort(GenomeBasicFeatures.begin(), GenomeBasicFeatures.end(), [](GenomeBasicFeature a, GenomeBasicFeature b) {
                            return a.getStart() < b.getStart();
                        });
                        for (std::vector<GenomeBasicFeature>::iterator it4 = GenomeBasicFeatures.begin();
                             it4 != GenomeBasicFeatures.end(); ++it4) {
                            ofile << it->getChromeSomeName() << "\t" + it->getSource() + "\t"+(*it4).getType()+"\t" << (*it4).getStart() << "\t"
                                  << (*it4).getEnd() << "\t" << it->getScore() << "\t" << st << "\t" << it4->getCodonFrame() << "\t"
                                  << it4->getLastColumnInformation() << std::endl;
                        }
                    }
                    if( it->getCdsVector().size() > 0 ){
                        ofile << "#metainformation: " << it->getMetaInformation() << std::endl;
                        ofile << "#genome sequence: " << it->getGeneomeSequence() << std::endl;
                        ofile << "#CDS sequence: " << it->getCdsSequence() << std::endl;
                    }
                }
            }
            ++index;
        }
    }
    ofile.close();
}
