//
// Created by baoxing on 10/10/17.
//

#include "myReAnnotationAndExonerateAndNovo.h"


void myReAnnotationAndExonerateAndNovo( std::string& referenceGenomeFile, std::string& inputGffFile,std::string& novoGffFilePath,
                                        std::string& variantsFile, std::string& outputGffFile, int &maxThread, int& lengthThread, std::string & vcfFix, std::map<std::string, std::string>& parameters, int & minIntron, bool & remove_reference_orf_shift){

    std::map<std::string, std::vector<Gene> > genes;
    std::map<std::string, Transcript > targetTranscriptsHashMap;
    reAnnotationAndExonerateAndNovo( referenceGenomeFile, inputGffFile, novoGffFilePath, variantsFile, genes, targetTranscriptsHashMap, maxThread, outputGffFile, lengthThread, vcfFix, parameters, minIntron, remove_reference_orf_shift);
    hereOutPutLiftOrOrthologousResult(genes, targetTranscriptsHashMap, outputGffFile );
}


void hereOutPutLiftOrOrthologousResult(std::map<std::string, std::vector<Gene> >& genes,
        std::map<std::string, Transcript> & targetTranscriptsHashMap, std::string& outputGffFile ){

    std::ofstream ofile;
    ofile.open(outputGffFile);

    for( std::map<std::string, std::vector<Gene> >::iterator it=genes.begin(); it!=genes.end(); ++it ){
        std::sort(it->second.begin(), it->second.end(), compare_gene());
        for( std::vector<Gene>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){
            std::string st = "+";

            if( NEGATIVE == (*it2).getStrand() ){
                st="-";
            }
            std::string geneResource = (*it2).getSource();
            if( geneResource.length()<1 ){
                geneResource="LIFTOVER";
            }
            ofile << it->first << "\t"+geneResource+"\tgene\t" << (*it2).getStart() << "\t" << (*it2).getEnd()
                  <<"\t.\t"<< st <<"\t.\tID="<< (*it2).getName()<< std::endl;
            for( std::vector<std::string>::iterator it3 = (*it2).getTranscriptVector().begin();
                 it3!=(*it2).getTranscriptVector().end(); ++it3  ){

                std::string transcriptResource = targetTranscriptsHashMap[*it3].getSource();
                if( transcriptResource.length()<1 ){
                    transcriptResource="LIFTOVER";
                }

                ofile << it->first << "\t"+transcriptResource+"\tmRNA\t" << targetTranscriptsHashMap[*it3].getPStart() << "\t" <<
                              targetTranscriptsHashMap[*it3].getPEnd() << "\t.\t"<< st <<"\t.\tID="<<
                              targetTranscriptsHashMap[*it3].getName() << ";Parent=" << (*it2).getName() << std::endl;
                //int cdsId = 1;
                for( std::vector<GenomeBasicFeature>::iterator it4=targetTranscriptsHashMap[*it3].getCdsVector().begin();
                     it4!=targetTranscriptsHashMap[*it3].getCdsVector().end(); ++it4 ){
                    ofile << it->first << "\t"+transcriptResource+"\tCDS\t" << (*it4).getStart() << "\t" <<
                          (*it4).getEnd() << "\t.\t" <<st << "\t.\tParent=" << targetTranscriptsHashMap[*it3].getName()<< std::endl;
                    // ++cdsId;
                }
                ofile << "#metainformation: " << targetTranscriptsHashMap[*it3].getMetaInformation() << std::endl;
                ofile << "#genome sequence: " << targetTranscriptsHashMap[*it3].getGeneomeSequence() << std::endl;
                ofile << "#CDS sequence: " << targetTranscriptsHashMap[*it3].getCdsSequence() << std::endl;
                std::string cdsSequence = targetTranscriptsHashMap[*it3].getCdsSequence();
                //std::string proteinSequence = nA2AA(cdsSequence, nucleotideCodeSubstitutionMatrix );
                //ofile << "#protein sequence: " << proteinSequence << std::endl;
            }
        }
    }
    ofile.close();
}
