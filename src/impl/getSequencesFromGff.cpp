//
// Created by baoxing on 10/10/17.
//

#include "getSequencesFromGff.h"
#include "CheckAndUpdateTranscriptsEnds.h"

void getSequences(const std::string& gffFile, const std::string& genomeFile, const std::string& outputProteinSequences,
                  const std::string& outputCdsSequences, const std::string& outputGenomeSequences, std::map<std::string, std::string>& parameters, const int & minIntron){
    std::string regex = get_parameters("cdsParentRegex", parameters);
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);

    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
    readGffFile (gffFile, transcriptHashSet, regex);
    std::map<std::string, Fasta> genome;
    readFastaFile(genomeFile, genome);
    CheckAndUpdateTranscriptsEnds( transcriptHashSet, genome, nucleotideCodeSubstitutionMatrix, minIntron);

    std::ofstream oPfile;
    oPfile.open(outputProteinSequences);
    std::ofstream oCfile;
    oCfile.open(outputCdsSequences);
    std::ofstream oGfile;
    oGfile.open(outputGenomeSequences);

    for( std::map<std::string, std::vector<Transcript> >::iterator it1=transcriptHashSet.begin();
            it1!=transcriptHashSet.end(); ++it1 ){
        if( genome.find(it1->first) != genome.end() ){
            for ( std::vector<Transcript>::iterator it2=it1->second.begin(); it2!=it1->second.end(); ++it2 ) {
                TranscriptUpdateCdsInformation((*it2), genome);
                checkOrfState( (*it2), genome, nucleotideCodeSubstitutionMatrix, minIntron);

                std::string cdsSequence = (*it2).getCdsSequence();

                oPfile << ">" << (*it2).getName() << std::endl;
                oCfile << ">" << (*it2).getName() << " metaInformation:" <<(*it2).getMetaInformation() << std::endl;
                oGfile << ">" << (*it2).getName() << std::endl;

                oCfile << (*it2).getCdsSequence() << std::endl;
                oGfile << (*it2).getGeneomeSequence() << std::endl;

                if( cdsSequence.length() >2 ){
                    std::string proteinSequence = nA2AA(cdsSequence, nucleotideCodeSubstitutionMatrix);
                    oPfile << proteinSequence << std::endl;
                }
            }
        }
    }
    oPfile.close();
    oCfile.close();
    oGfile.close();
}
