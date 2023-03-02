//
// Created by baoxing on 10/10/17.
//

#include "getSequencesFromGff.h"
#include "CheckAndUpdateTranscriptsEnds.h"

#include <stdexcept>

void getSequences(const std::string& gffFile, const std::string& genomeFile, const std::string& outputProteinSequences,
                  const std::string& outputCdsSequences, const std::string& outputGenomeSequences,
                  std::map<std::string, std::string>& parameters, const int & minIntron,
                  const bool & autoAdjust, const size_t & orfCount){
    const bool writeProtein = !outputProteinSequences.empty();
    const bool writeCds = !outputCdsSequences.empty();
    const bool writeGenome = !outputGenomeSequences.empty();
    if( !writeProtein && !writeCds && !writeGenome ){
        throw std::invalid_argument("at least one sequence output is required");
    }

    std::string regex = get_parameters("cdsParentRegex", parameters);
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);

    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
    readGffFile (gffFile, transcriptHashSet, regex);
    std::map<std::string, Fasta> genome;
    readFastaFile(genomeFile, genome);
    if( autoAdjust ){
        CheckAndUpdateTranscriptsEnds( transcriptHashSet, genome, nucleotideCodeSubstitutionMatrix, minIntron, orfCount);
    }

    std::ofstream oPfile;
    std::ofstream oCfile;
    std::ofstream oGfile;
    if( writeProtein ){
        oPfile.open(outputProteinSequences);
    }
    if( writeCds ){
        oCfile.open(outputCdsSequences);
    }
    if( writeGenome ){
        oGfile.open(outputGenomeSequences);
    }

    for( std::map<std::string, std::vector<Transcript> >::iterator it1=transcriptHashSet.begin();
            it1!=transcriptHashSet.end(); ++it1 ){
        if( genome.find(it1->first) != genome.end() ){
            for ( std::vector<Transcript>::iterator it2=it1->second.begin(); it2!=it1->second.end(); ++it2 ) {
                TranscriptUpdateCdsInformation((*it2), genome);
                checkOrfState( (*it2), genome, nucleotideCodeSubstitutionMatrix, minIntron);

                std::string cdsSequence = (*it2).getCdsSequence();

                if( writeProtein ){
                    oPfile << ">" << (*it2).getName() << std::endl;
                }
                if( writeCds ){
                    oCfile << ">" << (*it2).getName() << " metaInformation:" <<(*it2).getMetaInformation() << std::endl;
                    oCfile << (*it2).getCdsSequence() << std::endl;
                }
                if( writeGenome ){
                    oGfile << ">" << (*it2).getName() << std::endl;
                    oGfile << (*it2).getGeneomeSequence() << std::endl;
                }

                if( writeProtein && cdsSequence.length() >2 ){
                    std::string proteinSequence = nA2AA(cdsSequence, nucleotideCodeSubstitutionMatrix);
                    oPfile << proteinSequence << std::endl;
                }
            }
        }
    }
    if( writeProtein ){
        oPfile.close();
    }
    if( writeCds ){
        oCfile.close();
    }
    if( writeGenome ){
        oGfile.close();
    }
}
