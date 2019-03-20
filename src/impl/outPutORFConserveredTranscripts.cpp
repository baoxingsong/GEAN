//
// Created by Baoxing Song on 2019-03-17.
//

#include "outPutORFConserveredTranscripts.h"

void outPutORFConserveredTranscripts( const std::string & genomeFile, const std::string & gffFile,
                                      const std::string & outputGffFile, const int & minIntron,
                                      std::map<std::string, std::string>& parameters ){
     std::map<std::string, Fasta> genome;
     readFastaFile(genomeFile, genome);
     NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
     std::map<std::string, std::vector<Gene> > genes;
     readGffFileWithEveryThing (gffFile, genes);
     std::ofstream ofile;
     ofile.open(outputGffFile);
     for( std::map<std::string, std::vector<Gene> >::iterator it=genes.begin(); it!=genes.end(); ++it ){
         if(genome.find(it->first)!=genome.end()){
             for( Gene gene : it->second ){
                 if( gene.getTranscripts().size() > 0 ) {
                     double thisScore = 0;

                     for (int index = 0; index < gene.getTranscripts().size(); ++index) {
                         Transcript transcript = gene.getTranscripts()[index];
                         TranscriptUpdateCdsInformation(transcript, genome);
                         checkOrfState(transcript, genome, nucleotideCodeSubstitutionMatrix, minIntron);
                         if (transcript.getIfOrfShift()) {
                             // this score is for gene, if there is one ORF conserved transcript, then give this gene a positive score
                         } else {
                             thisScore += 1.0;
                         }
                         gene.getTranscripts()[index] = transcript;
                     }
                     if( thisScore>0.0 ){
                         for (std::vector<Transcript>::iterator it = gene.getTranscripts().begin();
                          it != gene.getTranscripts().end(); ++it) {
                             if (it->getIfOrfShift()) {

                             } else {
                                 int thisTranscriptStart = std::numeric_limits<int>::max();
                                 int thisTranscriptEnd = 0;
                                 if (it->getExonVector().size() > 0) {
                                     for (std::vector<GenomeBasicFeature>::iterator it4 = it->getExonVector().begin();
                                          it4 != it->getExonVector().end(); ++it4) {
                                         if (thisTranscriptStart > it4->getStart()) {
                                             thisTranscriptStart = it4->getStart();
                                         }
                                         if (thisTranscriptEnd < it4->getEnd()) {
                                             thisTranscriptEnd = it4->getEnd();
                                         }
                                     }
                                 }
                                 if (it->getCdsVector().size() > 0) {
                                     for (std::vector<GenomeBasicFeature>::iterator it4 = it->getCdsVector().begin();
                                          it4 != it->getCdsVector().end(); ++it4) {
                                         if (thisTranscriptStart > it4->getStart()) {
                                             thisTranscriptStart = it4->getStart();
                                         }
                                         if (thisTranscriptEnd < it4->getStart()) {
                                             thisTranscriptEnd = it4->getEnd();
                                         }
                                     }
                                 }
                                 it->setStart(thisTranscriptStart);
                                 it->setEnd(thisTranscriptEnd);
                             }
                         }
                         int thisStart = INT32_MAX;
                         int thisEnd = 0;
                         for (std::vector<Transcript>::iterator it = gene.getTranscripts().begin();
                              it != gene.getTranscripts().end(); ++it) {
                             if (it->getIfOrfShift()) {

                             } else {
                                 if (thisStart > it->getStart()) {
                                     thisStart = it->getStart();
                                 }
                                 if (thisEnd < it->getEnd()) {
                                     thisEnd = it->getEnd();
                                 }
                             }
                         }
                         std::string st = "+";
                         if (NEGATIVE == gene.getTranscripts()[0].getStrand()) {
                             st = "-";
                         }
                         std::string transcriptResource = gene.getTranscripts()[0].getSource();
                         if (transcriptResource.length() < 1) {
                             transcriptResource = "LIFTOVER";
                         }
                         ofile << gene.getTranscripts()[0].getChromeSomeName() << "\t" + transcriptResource + "\tgene\t"
                               << thisStart << "\t" <<
                               thisEnd << "\t.\t" << st << "\t.\tID=" << gene.getName() << ";" << std::endl;
                         for (std::vector<Transcript>::iterator it = gene.getTranscripts().begin();
                              it != gene.getTranscripts().end(); ++it) {
                             if (it->getIfOrfShift()) {

                             } else {
                                 transcriptResource = it->getSource();
                                 if (transcriptResource.length() < 1) {
                                     transcriptResource = "LIFTOVER";
                                 }

                                 ofile << it->getChromeSomeName() << "\t" + transcriptResource + "\t" << it->getType()
                                       << "\t" << it->getStart() << "\t" <<
                                       it->getEnd() << "\t" << it->getScore() << "\t" << st << "\t.\tID=" << it->getName()
                                       << ";Parent=" << it->getName() << ";" << std::endl;
                                 std::vector<GenomeBasicFeature> GenomeBasicFeatures;
                                 if (it->getFivePrimerUtr().size() > 0) {
                                     for (std::vector<GenomeBasicFeature>::iterator it4 = it->getFivePrimerUtr().begin();
                                          it4 != it->getFivePrimerUtr().end(); ++it4) {
                                         GenomeBasicFeatures.push_back(*it4);
                                     }
                                 }
                                 if (it->getCdsVector().size() > 0) {
                                     for (std::vector<GenomeBasicFeature>::iterator it4 = it->getCdsVector().begin();
                                          it4 != it->getCdsVector().end(); ++it4) {
                                         GenomeBasicFeatures.push_back(*it4);
                                     }
                                 }
                                 if (it->getExonVector().size() > 0) {
                                     for (std::vector<GenomeBasicFeature>::iterator it4 = it->getExonVector().begin();
                                          it4 != it->getExonVector().end(); ++it4) {
                                         GenomeBasicFeatures.push_back(*it4);
                                     }
                                 }
                                 if (it->getThreePrimerUtr().size() > 0) {
                                     for (std::vector<GenomeBasicFeature>::iterator it4 = it->getThreePrimerUtr().begin();
                                          it4 != it->getThreePrimerUtr().end(); ++it4) {
                                         GenomeBasicFeatures.push_back(*it4);
                                     }
                                 }
                                 if (GenomeBasicFeatures.size() > 0) {
                                     std::sort(GenomeBasicFeatures.begin(), GenomeBasicFeatures.end(),
                                               [](GenomeBasicFeature a, GenomeBasicFeature b) {
                                                   return a.getStart() < b.getStart();
                                               });
                                     for (std::vector<GenomeBasicFeature>::iterator it4 = GenomeBasicFeatures.begin();
                                          it4 != GenomeBasicFeatures.end(); ++it4) {
                                         ofile << it->getChromeSomeName()
                                               << "\t" + transcriptResource + "\t" + (*it4).getType() + "\t"
                                               << (*it4).getStart() << "\t"
                                               << (*it4).getEnd() << "\t" << it->getScore() << "\t" << st << "\t"
                                               << it4->getCodonFrame() << "\tParent="
                                               << it->getName() << ";" << std::endl;
                                     }
                                 }
                                 if (it->getCdsVector().size() > 0) {
                                     ofile << "#metainformation: " << it->getMetaInformation() << std::endl;
                                     ofile << "#genome sequence: " << it->getGeneomeSequence() << std::endl;
                                     ofile << "#CDS sequence: " << it->getCdsSequence() << std::endl;
                                 }
                                 std::string cdsSequence = it->getCdsSequence();
                             }
                         }
                    }
                 }
             }
         }
     }
    ofile.close();
}
