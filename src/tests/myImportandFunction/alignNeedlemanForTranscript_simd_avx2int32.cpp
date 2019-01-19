//
// Created by song on 8/2/18.
//


#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include "../../myImportandFunction/alignNeedlemanForTranscript_simd_avx2int32.h"

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>


TEST(alignNeedlemanForTranscript_simd_avx2int32, c1){
    std::string parameterFile = "configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "./");
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    int minIntron = 5;

    std::string referenceGenomeFilePath = "/home/song/Dropbox/tair11/Col.fa";
    std::map<std::string, Fasta> referenceGenome;

    readFastaFile(referenceGenomeFilePath, referenceGenome);
    std::cout << "reference genome sequence reading done" << std::endl;

    std::string referenceGffFilePath = "/home/song/Dropbox/tair11/TAIR10_GFF3_genes_no_UM.gff";
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    std::string regex = get_parameters("cdsParentRegex", parameters);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);

    CheckAndUpdateTranscriptsEnds( referenceTranscriptHashSet, referenceGenome, nucleotideCodeSubstitutionMatrix, minIntron);
    std::cout << "gff file reading done" << std::endl;
    std::map<std::string, std::vector<Variant> > variantsMaps;
    std::string sdiFile = "/home/song/Dropbox/tair11/PA9996.sdi";
    readSdiFile (sdiFile, variantsMaps, "", referenceGenome);
    std::cout << "variant tables reading done" << std::endl;

    std::map<std::string, Fasta> targetGenome;
    getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome, parameters);
    std::cout << "get pseudoGenome sequence done" << std::endl;
    std::map<std::string, Transcript> targetTranscriptsHashMap;// transcriptName Transcript
    annotationLiftOver(referenceTranscriptHashSet, targetTranscriptsHashMap, variantsMaps, targetGenome, referenceGenome, parameters, minIntron);
    std::cout << "annotationLiftOver done" << std::endl;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    for(std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
        it!=referenceTranscriptHashSet.end(); ++it) {
        if( referenceGenome.find(it->first) != referenceGenome.end() ) {
            for (std::vector<Transcript>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                Transcript referenceTranscript = (*it2);
                Transcript targetTranscript = targetTranscriptsHashMap[it2->getName()];
                TranscriptUpdateCdsInformation(referenceTranscript, referenceGenome);
                std::string refGenomeSequence = referenceTranscript.getGeneomeSequence();

                int startCodonPosition = 1;
                int stopCodonPosition = refGenomeSequence.length() - 2;
                std::vector<SpliceSitePosition> spliceSitePositions;

                targetTranscript.setSource("REALIGNMENT");
                //prepare for the new GenomeBasicFeature After re-alignment begin
                std::vector<int> targetCdsStarts;
                std::vector<int> targetCdsEnds;
                STRAND thisStrand = referenceTranscript.getStrand();

                int startTarget = targetTranscript.getPStart() - refGenomeSequence.length();
                int endTarget = targetTranscript.getPEnd() + refGenomeSequence.length();

                std::string dna_b = getSubsequence(targetGenome, it->first, startTarget, endTarget, thisStrand);
                if (referenceTranscript.getStrand() == POSITIVE && referenceTranscript.getCdsVector().size() > 1) {
                    if (referenceTranscript.getCdsVector().size() > 1) {
                        for (size_t i = 1; i < referenceTranscript.getCdsVector().size(); ++i) {
                            SpliceSitePosition spliceSitePosition(
                                    referenceTranscript.getCdsVector()[i - 1].getEnd() -
                                    referenceTranscript.getPStart() + 2,
                                    referenceTranscript.getCdsVector()[i].getStart() - referenceTranscript.getPStart());
                            spliceSitePositions.push_back(spliceSitePosition);
                        }
                    }
                    if (dna_b.length() < 10000) {
                        alignNeedlemanForTranscript_simd_avx2int32 nw(refGenomeSequence, dna_b, startCodonPosition,
                                                                      stopCodonPosition, spliceSitePositions,
                                                                      parameters,
                                                                      nucleotideCodeSubstitutionMatrix);
//                        std::cout << "standard " << referenceTranscript.getName() << " begin" << std::endl;
                        NeedlemanWunschForTranscript nw2(refGenomeSequence, dna_b, startCodonPosition,
                                                         stopCodonPosition,
                                                         spliceSitePositions, parameters,
                                                         nucleotideCodeSubstitutionMatrix);
                        std::string rq = nw.getAlignment_q();
                        rq.erase(std::remove(rq.begin(), rq.end(), '-'), rq.end());
                        std::string rd = nw.getAlignment_d();
                        rd.erase(std::remove(rd.begin(), rd.end(), '-'), rd.end());

                        if (nw.getAlignment_d().compare(nw2.getAlignment_d()) != 0) {
                            std::cout << "error d" << it2->getName() << std::endl;
                            nw.print_results();
                            nw2.print_results();
                        } else if (nw.getAlignment_q().compare(nw2.getAlignment_q()) != 0) {
                            std::cout << "error q" << it2->getName() << std::endl;
                            nw.print_results();
                            nw2.print_results();
                        } else if (rq.compare(dna_b) != 0) {
                            std::cout << "error rq" << it2->getName() << std::endl;
                            nw.print_results();
                            nw2.print_results();
                        } else if (rd.compare(refGenomeSequence) != 0) {
                            std::cout << "error rd" << it2->getName() << std::endl;
                            nw.print_results();
                            nw2.print_results();
                        }else{
                            std::cout << "good" << it2->getName() << std::endl;
                            nw.print_results();
                            nw2.print_results();
                        }
                    }
                }
            }
        }
    }
    ASSERT_EQ(0, 0);
}


TEST(alignNeedlemanForTranscript_simd_avx2int32, c2){
    std::string parameterFile = "configure";
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, "./");
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    int minIntron = 5;

    std::string referenceGenomeFilePath = "/biodata/dep_tsiantis/grp_gan/song/rdINDELallHere/inputData/fullSdiFile/col_0.fa";
    std::map<std::string, Fasta> referenceGenome;

    readFastaFile(referenceGenomeFilePath, referenceGenome);
    std::cout << "reference genome sequence reading done" << std::endl;

    std::string referenceGffFilePath = "/biodata/dep_tsiantis/grp_gan/song/gennepredication/TAIR10annotationclassification/unmodifyed/TAIR10_GFF3_genes.gff";
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    std::string regex = get_parameters("cdsParentRegex", parameters);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);

    CheckAndUpdateTranscriptsEnds( referenceTranscriptHashSet, referenceGenome, nucleotideCodeSubstitutionMatrix, minIntron);
    std::cout << "gff file reading done" << std::endl;
    std::map<std::string, std::vector<Variant> > variantsMaps;
    std::string sdiFile = "/biodata/dep_tsiantis/grp_gan/song/rdINDELallHere/inputData/fullSdiFile/PA9830.sdi";
    readSdiFile (sdiFile, variantsMaps, "", referenceGenome);
    std::cout << "variant tables reading done" << std::endl;

    std::map<std::string, Fasta> targetGenome;
    getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome, parameters);
    std::cout << "get pseudoGenome sequence done" << std::endl;
    std::map<std::string, Transcript> targetTranscriptsHashMap;// transcriptName Transcript
    annotationLiftOver(referenceTranscriptHashSet, targetTranscriptsHashMap, variantsMaps, targetGenome, referenceGenome, parameters, minIntron);
    std::cout << "annotationLiftOver done" << std::endl;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    for(std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
        it!=referenceTranscriptHashSet.end(); ++it) {
        if( referenceGenome.find(it->first) != referenceGenome.end() ) {
            for (std::vector<Transcript>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                Transcript referenceTranscript = (*it2);
                Transcript targetTranscript = targetTranscriptsHashMap[it2->getName()];
                TranscriptUpdateCdsInformation(referenceTranscript, referenceGenome);
                std::string refGenomeSequence = referenceTranscript.getGeneomeSequence();

                int startCodonPosition = 1;
                int stopCodonPosition = refGenomeSequence.length() - 2;
                std::vector<SpliceSitePosition> spliceSitePositions;

                targetTranscript.setSource("REALIGNMENT");
                int targetPosition = 0;
                int referencePosition = 0;
                //prepare for the new GenomeBasicFeature After re-alignment begin
                std::vector<int> targetCdsStarts;
                std::vector<int> targetCdsEnds;
                STRAND thisStrand = referenceTranscript.getStrand();

                int startTarget = targetTranscript.getPStart() - refGenomeSequence.length();
                int endTarget = targetTranscript.getPEnd() + refGenomeSequence.length();

                std::string dna_b = getSubsequence(targetGenome, it->first, startTarget, endTarget, thisStrand);
                if (referenceTranscript.getStrand() == POSITIVE && referenceTranscript.getCdsVector().size() > 1) {
                    if (referenceTranscript.getCdsVector().size() > 1) {
                        for (size_t i = 1; i < referenceTranscript.getCdsVector().size(); ++i) {
                            SpliceSitePosition spliceSitePosition(
                                    referenceTranscript.getCdsVector()[i - 1].getEnd() -
                                    referenceTranscript.getPStart() + 2,
                                    referenceTranscript.getCdsVector()[i].getStart() - referenceTranscript.getPStart());
                            spliceSitePositions.push_back(spliceSitePosition);
                        }
                    }
                    if (dna_b.length() < 5000) {
                        alignNeedlemanForTranscript_simd_avx2int32 nw(refGenomeSequence, dna_b, startCodonPosition,
                                                                      stopCodonPosition, spliceSitePositions,
                                                                      parameters,
                                                                      nucleotideCodeSubstitutionMatrix);
                    }
                }
            }
        }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    for(std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
        it!=referenceTranscriptHashSet.end(); ++it) {
        if( referenceGenome.find(it->first) != referenceGenome.end() ) {
            for (std::vector<Transcript>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                Transcript referenceTranscript = (*it2);
                Transcript targetTranscript = targetTranscriptsHashMap[it2->getName()];
                TranscriptUpdateCdsInformation(referenceTranscript, referenceGenome);
                std::string refGenomeSequence = referenceTranscript.getGeneomeSequence();

                int startCodonPosition = 1;
                int stopCodonPosition = refGenomeSequence.length() - 2;
                std::vector<SpliceSitePosition> spliceSitePositions;

                targetTranscript.setSource("REALIGNMENT");
                int targetPosition = 0;
                int referencePosition = 0;
                //prepare for the new GenomeBasicFeature After re-alignment begin
                std::vector<int> targetCdsStarts;
                std::vector<int> targetCdsEnds;
                STRAND thisStrand = referenceTranscript.getStrand();

                int startTarget = targetTranscript.getPStart() - refGenomeSequence.length();
                int endTarget = targetTranscript.getPEnd() + refGenomeSequence.length();

                std::string dna_b = getSubsequence(targetGenome, it->first, startTarget, endTarget, thisStrand);
                if (referenceTranscript.getStrand() == POSITIVE && referenceTranscript.getCdsVector().size() > 1) {
                    if (referenceTranscript.getCdsVector().size() > 1) {
                        for (size_t i = 1; i < referenceTranscript.getCdsVector().size(); ++i) {
                            SpliceSitePosition spliceSitePosition(
                                    referenceTranscript.getCdsVector()[i - 1].getEnd() -
                                    referenceTranscript.getPStart() + 2,
                                    referenceTranscript.getCdsVector()[i].getStart() - referenceTranscript.getPStart());
                            spliceSitePositions.push_back(spliceSitePosition);
                        }
                    }
                    if (dna_b.length() < 5000) {
                        NeedlemanWunschForTranscript nw(refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition,
                                                        spliceSitePositions, parameters,
                                                        nucleotideCodeSubstitutionMatrix);
                    }
                }
            }
        }
    }
    std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_span1 = t2 - t1;
    std::chrono::duration<double, std::milli> time_span2 = t3 - t2;
    std::cout << "time costing " << time_span1.count() << std::endl;
    std::cout << "time costing " << time_span2.count()<< std::endl;
    ASSERT_EQ(0, 0);
}
