//
// Created by song on 8/4/18.
//

#include "TransferGffWithNucmerResult.h"

struct NewGffRecord{
    std::string geneName;
    STRAND strand;
    std::vector<Transcript> transcripts;
};

bool overlap2( const Transcript & transcript, const Range & range ){
    if( ( transcript.getPStart()<=range.getStart() && range.getStart()<=transcript.getPEnd()) ||
        ( transcript.getPStart()<=range.getEnd() && range.getEnd()<=transcript.getPEnd()) ||
        ( range.getStart()<=transcript.getPStart() && transcript.getPStart()<=range.getEnd())||
        ( range.getStart()<=transcript.getPEnd() && transcript.getPEnd()<=range.getEnd())  ){
        return true;
    }
    return false;
}

void getAlltheOverLappedTranscripts2( std::map<std::string, std::vector<Transcript> > & transcriptHashSet,
                                      const std::string & chr, const AlignmentMatch & alignmentMatch, std::vector<Transcript*> & overLappedTranscripts,
                                      size_t & startShitfDistance, size_t & endShiftDistance){
    startShitfDistance=0;
    endShiftDistance=0;
    int i=0, j=transcriptHashSet[chr].size()-1, newi, newj, k=10;
    bool iChanged=true, jChanged=true;
    int interVal = 3;
    while( k>0 && (iChanged || jChanged)  ){
        --k;
        iChanged=false;
        jChanged=false;
        if( transcriptHashSet[chr][i].getPEnd() < alignmentMatch.getDatabase().getStart() ){
            newi = i+(j-i)/interVal;
            if( transcriptHashSet[chr][newi].getPEnd() < alignmentMatch.getDatabase().getStart() ){
                iChanged=true;
                i=newi;
            }
        }
        if( transcriptHashSet[chr][j].getPStart() > alignmentMatch.getDatabase().getEnd()  ){
            newj = j-(j-i)/interVal;
            if( transcriptHashSet[chr][newj].getPStart() > alignmentMatch.getDatabase().getEnd() ){
                jChanged=true;
                j=newj;
            }
        }
        if( !iChanged && !jChanged ){
            interVal=2*interVal;
            iChanged=true;
        }
    }
    i-=10;
    j+=10;
    if( i<0 ){
        i=0;
    }
    if( j >= transcriptHashSet[chr].size() ){
        j=transcriptHashSet[chr].size()-1;
    }
    for( ; i<=j; ++i ){
        if( overlap2(transcriptHashSet[chr][i], alignmentMatch.getDatabase()) ){
            if( transcriptHashSet[chr][i].getPStart() < alignmentMatch.getDatabase().getStart() ){
                if( startShitfDistance < (transcriptHashSet[chr][i].getPEnd()-transcriptHashSet[chr][i].getPStart()) ){
                    startShitfDistance = (transcriptHashSet[chr][i].getPEnd()-transcriptHashSet[chr][i].getPStart());
                }
            }
            if( transcriptHashSet[chr][i].getPEnd() > alignmentMatch.getDatabase().getEnd() ){
                if( endShiftDistance < (transcriptHashSet[chr][i].getPEnd()-transcriptHashSet[chr][i].getPStart()) ){
                    endShiftDistance = (transcriptHashSet[chr][i].getPEnd()-transcriptHashSet[chr][i].getPStart());
                }
            }
            overLappedTranscripts.push_back(&transcriptHashSet[chr][i]);
        }
    }
}

void algnmentMatchUpdate2( size_t & startShitfDistance, size_t & endShiftDistance, AlignmentMatch & alignmentMatch,
                           std::map<std::string, Fasta> & databaseSequences, std::map<std::string, Fasta> & querySequences) {
    size_t queryStart;
    size_t queryEnd;
    size_t databasetStart;
    size_t databaseEnd;
    endShiftDistance = 2 * (endShiftDistance + 1);
    databaseEnd = alignmentMatch.getDatabaseEnd() + endShiftDistance;
    if (databaseEnd >
        databaseSequences[alignmentMatch.getDatabaseChr()].getSequence().length()) {
        databaseEnd = databaseSequences[alignmentMatch.getDatabaseChr()].getSequence().length();
    }
    alignmentMatch.setDatabaseEnd(databaseEnd);

    startShitfDistance = 2 * (startShitfDistance + 1);
    if (alignmentMatch.getDatabaseStart() <= startShitfDistance) {
        databasetStart = 1;
    } else {
        databasetStart = alignmentMatch.getDatabaseStart() - startShitfDistance;
    }
    alignmentMatch.setDatabaseStart(databasetStart);
    if (POSITIVE == alignmentMatch.getQueryStrand()){
        queryEnd =
                alignmentMatch.getQueryEnd() + endShiftDistance; // here it is related with strand
        if (alignmentMatch.getQueryStart() <= startShitfDistance) {
            queryStart = 1;
        } else {
            queryStart = alignmentMatch.getQueryStart() - startShitfDistance;
        }
        //update algnmentMatch end
    }else{
        if (alignmentMatch.getQueryStart() <= endShiftDistance) {
            queryStart = 1;
        } else {
            queryStart = alignmentMatch.getQueryStart()-endShiftDistance;
        }
        queryEnd =
                alignmentMatch.getQueryEnd() + startShitfDistance; // here it is related with strand
        //update algnmentMatch end
    }
    alignmentMatch.setQueryStart(queryStart);
    if (queryEnd > querySequences[alignmentMatch.getQueryChr()].getSequence().length()) {
        queryEnd = querySequences[alignmentMatch.getQueryChr()].getSequence().length();
    }
    alignmentMatch.setQueryEnd(queryEnd);
}

bool overlap( const Transcript & transcript, const Range & range ){
    if( range.getStart()<=transcript.getPStart() && transcript.getPEnd()<=range.getEnd() ){
        return true;
    }
    return false;
}

void getAlltheOverLappedTranscripts( std::map<std::string, std::vector<Transcript> > & transcriptHashSet,
                                     const std::string & chr, const AlignmentMatch & alignmentMatch, std::vector<Transcript*> & overLappedTranscripts,
                                     size_t & startShitfDistance, size_t & endShiftDistance){
    startShitfDistance=0;
    endShiftDistance=0;
    int i=0, j=transcriptHashSet[chr].size()-1, newi, newj, k=10;
    bool iChanged=true, jChanged=true;
    int interVal = 3;
    while( k>0 && (iChanged || jChanged)  ){
        --k;
        iChanged=false;
        jChanged=false;
        if( transcriptHashSet[chr][i].getPEnd() < alignmentMatch.getDatabase().getStart() ){
            newi = i+(j-i)/interVal;
            if( transcriptHashSet[chr][newi].getPEnd() < alignmentMatch.getDatabase().getStart() ){
                iChanged=true;
                i=newi;
            }
        }
        if( transcriptHashSet[chr][j].getPStart() > alignmentMatch.getDatabase().getEnd()  ){
            newj = j-(j-i)/interVal;
            if( transcriptHashSet[chr][newj].getPStart() > alignmentMatch.getDatabase().getEnd() ){
                jChanged=true;
                j=newj;
            }
        }
        if( !iChanged && !jChanged ){
            interVal=2*interVal;
            iChanged=true;
        }
    }
    for( ; i<=j; ++i ){
        if( overlap(transcriptHashSet[chr][i], alignmentMatch.getDatabase()) ){
            overLappedTranscripts.push_back(&transcriptHashSet[chr][i]);
        }
    }
}


void slidingWinAlnAndGeneRateAnnotation(AlignmentMatch & alignmentMatch,
                                        std::map<std::string, Fasta> & databaseSequences,
                                        std::map<std::string, Fasta> & querySequences,
                                        std::vector<Transcript*> & overLappedTranscripts,
                                        std::string & alignQuerySequence,
                                        std::string & alignDatabaseSequence, std::map<int, int> & importantPositions,
                                        std::map<std::string, std::string>& parameters,
                                        NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix){
    int newStart;
    int queryPosition = 0;
    std::string databaseSequence = getSubsequence(databaseSequences, alignmentMatch.getDatabaseChr(),
                                                  alignmentMatch.getDatabaseStart(),
                                                  alignmentMatch.getDatabaseEnd());
    std::string querySequence = getSubsequence(querySequences, alignmentMatch.getQueryChr(),
                                               alignmentMatch.getQueryStart(),
                                               alignmentMatch.getQueryEnd(), alignmentMatch.getQueryStrand());
    alignSlidingWindow(querySequence, databaseSequence, alignQuerySequence, alignDatabaseSequence, 50, parameters, nucleotideCodeSubstitutionMatrix);
    //generate annotation according to sequence alignment begin

    for (Transcript *referenceTranscript : overLappedTranscripts) {
        for (std::vector<GenomeBasicFeature>::iterator it4 = referenceTranscript->getCdsVector().begin();
             it4 != referenceTranscript->getCdsVector().end(); ++it4) {
            importantPositions[(*it4).getStart()]=0;
            importantPositions[(*it4).getEnd()]=0;
        }
    }
    int databasePosition = alignmentMatch.getDatabaseStart() - 1;
    for (size_t tp = 0; tp < alignQuerySequence.length(); ++tp) {
        if (alignQuerySequence[tp] != '-') {
            ++queryPosition;
        }
        if (alignDatabaseSequence[tp] != '-') {
            ++databasePosition;
            if( importantPositions.find(databasePosition) != importantPositions.end() ){
                if (POSITIVE == alignmentMatch.getQueryStrand()) {
                    newStart = alignmentMatch.getQueryStart() + queryPosition - 1;
                } else {
                    newStart = alignmentMatch.getQueryEnd() - queryPosition + 1;
                }
                if (newStart < 1) {
                    newStart = 1;
                } else if (newStart >
                           querySequences[alignmentMatch.getQueryChr()].getSequence().length()) {
                    newStart = querySequences[alignmentMatch.getQueryChr()].getSequence().length();
                }
                importantPositions[databasePosition] = newStart;
            }
        }
    }
}

void geneStruAlnAndGeneRateAnnotation(Transcript & newTranscript, Transcript * referenceTranscript, Transcript & newTranscript2,
                                      std::map<std::string, Fasta> & databaseSequences,
                                      std::map<std::string, Fasta> & querySequences, AlignmentMatch & alignmentMatch,
                                      std::map<std::string, std::string>& parameters,
                                      NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix){
    int startCodonPosition = 1;
    int i, newStart;
    size_t transcriptLength = newTranscript.getPEnd() - newTranscript.getPStart();
    size_t startTarget;
    size_t endTarget;
    if (newTranscript.getPStart() > transcriptLength) {
        startTarget = newTranscript.getPStart() - transcriptLength;
    } else {
        startTarget = 1;
    }
    if ((newTranscript.getPEnd() + transcriptLength) >
        querySequences[alignmentMatch.getQueryChr()].getSequence().length()) {
        endTarget = querySequences[alignmentMatch.getQueryChr()].getSequence().length();
    } else {
        endTarget = newTranscript.getPEnd() + transcriptLength;
    }
    std::vector<SpliceSitePosition> spliceSitePositions;
    std::string refGenomeSequence = getSubsequence(databaseSequences,
                                                   referenceTranscript->getChromeSomeName(),
                                                   referenceTranscript->getPStart(),
                                                   referenceTranscript->getPEnd(),
                                                   referenceTranscript->getStrand());
    int stopCodonPosition = refGenomeSequence.length() - 2;
    std::string dna_b;
    if( POSITIVE == alignmentMatch.getQueryStrand() ) {
        dna_b = getSubsequence(querySequences,
                               alignmentMatch.getQueryChr(),
                               startTarget, endTarget, referenceTranscript->getStrand());
    } else {
        STRAND thisStrand;
        if (referenceTranscript->getStrand() == POSITIVE) {
            thisStrand = NEGATIVE;
        } else {
            thisStrand = POSITIVE;
        }
        dna_b = getSubsequence(querySequences,
                               alignmentMatch.getQueryChr(),
                               startTarget, endTarget, thisStrand);
    }
    if (referenceTranscript->getStrand() == POSITIVE) {
        if (referenceTranscript->getCdsVector().size() > 1) {
            for (i = 1; i < referenceTranscript->getCdsVector().size(); ++i) {
                SpliceSitePosition spliceSitePosition(
                        referenceTranscript->getCdsVector()[i - 1].getEnd() -
                        referenceTranscript->getPStart() + 2,
                        referenceTranscript->getCdsVector()[i].getStart() -
                        referenceTranscript->getPStart());
                spliceSitePositions.push_back(spliceSitePosition);
            }
        }
    }else{
        if (referenceTranscript->getCdsVector().size() > 1) {
            for (i = 1; i < referenceTranscript->getCdsVector().size(); ++i) {
                SpliceSitePosition spliceSitePosition(
                        referenceTranscript->getPEnd()-referenceTranscript->getCdsVector()[i].getStart()+2,
                        referenceTranscript->getPEnd()-referenceTranscript->getCdsVector()[i-1].getEnd());
                spliceSitePositions.push_back(spliceSitePosition);
            }
        }
    }
    AlignTranscript nw(refGenomeSequence, dna_b,
                                                  startCodonPosition,
                                                  stopCodonPosition,
                                                  spliceSitePositions,
                                                  parameters,
                                                  nucleotideCodeSubstitutionMatrix);
    std::string alignQuerySequence = nw.getAlignment_q();
    std::string alignDatabaseSequence = nw.getAlignment_d();
    int queryPosition = 0;
    std::map<int, int> importantPositions;
    for (std::vector<GenomeBasicFeature>::iterator it4 = referenceTranscript->getCdsVector().begin();
         it4 != referenceTranscript->getCdsVector().end(); ++it4) {
        importantPositions[(*it4).getStart()]=0;
        importantPositions[(*it4).getEnd()]=0;
    }
    if (referenceTranscript->getStrand() == POSITIVE) {
        int databasePosition =referenceTranscript->getPStart()-1;
        for (size_t tp = 0; tp < alignQuerySequence.length(); ++tp) {
            if (alignQuerySequence[tp] != '-') {
                ++queryPosition;
            }
            if (alignDatabaseSequence[tp] != '-') {
                ++databasePosition;
                if (importantPositions.find(databasePosition ) !=
                    importantPositions.end()) {
                    if( POSITIVE == alignmentMatch.getQueryStrand() ) {
                        newStart = startTarget + queryPosition - 1;
                    }else{
                        newStart = endTarget - queryPosition + 1;
                    }
                    if (newStart < 1) {
                        newStart = 1;
                    } else if (newStart >
                               querySequences[alignmentMatch.getQueryChr()].getSequence().length()) {
                        newStart = querySequences[alignmentMatch.getQueryChr()].getSequence().length();
                    }
                    importantPositions[databasePosition] = newStart;
                }
            }
        }
    }else{
        int databasePosition = referenceTranscript->getPEnd() + 1;
        for (size_t tp = 0; tp < alignQuerySequence.length(); ++tp) {
            if (alignQuerySequence[tp] != '-') {
                ++queryPosition;
            }
            if (alignDatabaseSequence[tp] != '-') {
                --databasePosition;
                if (importantPositions.find(databasePosition) !=
                    importantPositions.end()) {
                    if( POSITIVE == alignmentMatch.getQueryStrand() ) {
                        newStart = endTarget - queryPosition + 1;
                    }else{
                        newStart = startTarget + queryPosition - 1;
                    }
                    if (newStart < 1) {
                        newStart = 1;
                    } else if (newStart > querySequences[alignmentMatch.getQueryChr()].getSequence().length()) {
                        newStart = querySequences[alignmentMatch.getQueryChr()].getSequence().length();
                    }
                    importantPositions[databasePosition] = newStart;
                }
            }
        }
    }
    for (size_t i5 = 0;
         i5 < referenceTranscript->getCdsVector().size(); ++i5) {
        GenomeBasicFeature cds(importantPositions[referenceTranscript->getCdsVector()[i5].getStart()],
                               importantPositions[referenceTranscript->getCdsVector()[i5].getEnd()]);
        newTranscript2.addCds(cds);
    }
}
void newTranscriptAddCds(Transcript & referenceTranscript, Transcript & newTranscript, std::map<int, int> & importantPositions,
                         std::map<std::string, Fasta> & querySequences, NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix, int minIntron){
    if( referenceTranscript.getCdsVector().size() >0 ){
        for (size_t i5=0; i5 < referenceTranscript.getCdsVector().size(); ++i5) {
            GenomeBasicFeature cds(importantPositions[referenceTranscript.getCdsVector()[i5].getStart()],
                                   importantPositions[referenceTranscript.getCdsVector()[i5].getEnd()]);
            cds.setLastColumnInformation(referenceTranscript.getCdsVector()[i5].getLastColumnInformation());
            cds.setType(referenceTranscript.getCdsVector()[i5].getType());
            cds.setCodonFrame(referenceTranscript.getCdsVector()[i5].getCodonFrame());
            newTranscript.addCds(cds);
        }
        TranscriptUpdateCdsInformation(newTranscript, querySequences);
        checkOrfState(newTranscript, querySequences, nucleotideCodeSubstitutionMatrix, minIntron);
    }else{
        newTranscript.setIfOrfShift(false);
        newTranscript.setMetaInformation("this is not an coding transcript");
        newTranscript.setCdsSequence("NA");
    }
}
void TransferGffWithNucmerResult( const std::string & gffFilePath, const std::string & databaseFastaFilePath,
                                  const std::string & queryFastaFilePath, const std::string & nucmerFilePath,
                                  std::map<std::string, std::string>& parameters, const std::string & outPutFilePath,
                                  const size_t & maxLengthForStructureAlignment  ){

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    size_t minIntron = 5;

    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
    readGffFile(gffFilePath, transcriptHashSet);

    std::map<std::string, Fasta> databaseSequences;
    readFastaFile(databaseFastaFilePath, databaseSequences);
    CheckAndUpdateTranscriptsEnds( transcriptHashSet, databaseSequences, nucleotideCodeSubstitutionMatrix, minIntron);
    std::map<std::string, Fasta> querySequences;
    readFastaFile(queryFastaFilePath, querySequences);

    std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMap;
    nucmerRead(nucmerFilePath, alignmentMatchsMap);

    size_t startShitfDistance=0;
    size_t endShiftDistance=0;

    std::vector<Transcript> newTranscripts;
    int i;

    for(std::map<std::string, std::vector<AlignmentMatch>>::iterator it1=alignmentMatchsMap.begin();
        it1!=alignmentMatchsMap.end(); ++it1){
//        std::cout << it1->first << std::endl;
        if( databaseSequences.find(it1->first) != databaseSequences.end() ){
            for( AlignmentMatch alignmentMatch : it1->second ){
//                std::cout << "\t" << alignmentMatch.getDatabaseStart() << "\t" << alignmentMatch.getDatabaseEnd() << std::endl;
                if( querySequences.find(alignmentMatch.getQueryChr())!=querySequences.end() ){
                    //update algnmentMatch begin
                    std::vector<Transcript*> overLappedTranscripts;
//                    std::cout << "line 351" << std::endl;
                    getAlltheOverLappedTranscripts( transcriptHashSet, it1->first, alignmentMatch, overLappedTranscripts,
                                                    startShitfDistance, endShiftDistance);
//                    std::cout << "line 354" << std::endl;
                    if( overLappedTranscripts.size()>0 ) {
//                        std::cout << "line 356" << std::endl;
//                        algnmentMatchUpdate( startShitfDistance, endShiftDistance, alignmentMatch,
//                                             databaseSequences, querySequences);
//                        std::cout << "line 359" << std::endl;
                        std::string alignQuerySequence = "";
                        std::string alignDatabaseSequence = "";
                        std::map<int, int> importantPositions;
                        slidingWinAlnAndGeneRateAnnotation(alignmentMatch, databaseSequences, querySequences,
                                                           overLappedTranscripts,
                                                           alignQuerySequence, alignDatabaseSequence, importantPositions, parameters, nucleotideCodeSubstitutionMatrix);
//                        std::cout << "line 366" << std::endl;
                        if( POSITIVE == alignmentMatch.getQueryStrand() ) {
                            for ( i=0; i<overLappedTranscripts.size(); ++i) {
                                Transcript *referenceTranscript = overLappedTranscripts[i];
                                Transcript newTranscript(referenceTranscript->getName(),
                                                         alignmentMatch.getQueryChr(),
                                                         referenceTranscript->getStrand());
//                                std::cout << "line 373" << std::endl;
                                newTranscriptAddCds((*referenceTranscript), newTranscript, importantPositions, querySequences, nucleotideCodeSubstitutionMatrix, minIntron);
//                                std::cout << "line 375" << std::endl;
                                if ( !newTranscript.getIfOrfShift() ||
                                     ((referenceTranscript->getPEnd()-referenceTranscript->getPStart()) >= maxLengthForStructureAlignment)) {
                                    newTranscripts.push_back(newTranscript);
                                } else {
                                    Transcript newTranscript2(referenceTranscript->getName(),
                                                              alignmentMatch.getQueryChr(),
                                                              referenceTranscript->getStrand());
//                                    std::cout << "line 383" << std::endl;
                                    geneStruAlnAndGeneRateAnnotation( newTranscript, referenceTranscript, newTranscript2,
                                                                      databaseSequences,
                                                                      querySequences, alignmentMatch, parameters, nucleotideCodeSubstitutionMatrix);
//                                    std::cout << "line 387" << std::endl;
                                    TranscriptUpdateCdsInformation(newTranscript2, querySequences);
//                                    std::cout << "line 389" << std::endl;
                                    checkOrfState(newTranscript2, querySequences,
                                                  nucleotideCodeSubstitutionMatrix, minIntron);
//                                    std::cout << "line 392" << std::endl;
                                    if (newTranscript2.getIfOrfShift()) {
                                        newTranscripts.push_back(newTranscript);
                                    }else{
                                        newTranscript2.setSource("REALIGNMENT");
                                        newTranscripts.push_back(newTranscript2);
                                    }
                                }
                            }
                        } else {
                            for (Transcript *referenceTranscript : overLappedTranscripts) {
                                STRAND thisStrand;
                                if (referenceTranscript->getStrand() == POSITIVE) {
                                    thisStrand = NEGATIVE;
                                } else {
                                    thisStrand = POSITIVE;
                                }
                                Transcript newTranscript(referenceTranscript->getName(),
                                                         alignmentMatch.getQueryChr(), thisStrand);
//                                std::cout << "line 411" << std::endl;
                                newTranscriptAddCds((*referenceTranscript), newTranscript, importantPositions, querySequences, nucleotideCodeSubstitutionMatrix, minIntron);
//                                std::cout << "line 413" << std::endl;
//                                std::cout << "line 415" << std::endl;
                                if ( !newTranscript.getIfOrfShift()) {
                                    newTranscripts.push_back(newTranscript);
                                } else if( (referenceTranscript->getPEnd()-referenceTranscript->getPStart()) <
                                           maxLengthForStructureAlignment ) {
                                    Transcript newTranscript2(referenceTranscript->getName(),
                                                              alignmentMatch.getQueryChr(),
                                                              thisStrand);
//                                    std::cout << "line 423" << std::endl;
                                    geneStruAlnAndGeneRateAnnotation( newTranscript, referenceTranscript, newTranscript2,
                                                                      databaseSequences,
                                                                      querySequences, alignmentMatch, parameters, nucleotideCodeSubstitutionMatrix);
                                    TranscriptUpdateCdsInformation(newTranscript2, querySequences); // there is something wrong with here
//                                    std::cout << "line 425" << std::endl;
                                    checkOrfState(newTranscript2, querySequences, nucleotideCodeSubstitutionMatrix, minIntron);
//                                    std::cout << "line 427" << std::endl;
                                    if (newTranscript2.getIfOrfShift()) {
                                        newTranscripts.push_back(newTranscript);
                                    } else {
                                        newTranscript2.setSource("REALIGNMENT");
                                        newTranscripts.push_back(newTranscript2);
                                    }
                                }
                            }//generate annotation according to sequence alignment end
                        }
                    }
                }
            }
        }
    }
    std::ofstream ofile;

    ofile.open( outPutFilePath);
    for( std::vector<Transcript>::iterator it=newTranscripts.begin(); it!=newTranscripts.end(); ++it ){
        std::string transcriptResource = it->getSource();
        if( transcriptResource.length()<1 ){
            transcriptResource="LIFTOVER";
        }
        std::string st = "+";
        if( NEGATIVE == it->getStrand() ){
            st="-";
        }
        ofile << it->getChromeSomeName() << "\t"+transcriptResource+"\tmRNA\t" << it->getPStart() << "\t" <<
              it->getPEnd() << "\t.\t"<< st <<"\t.\tID="<< it->getName() << std::endl;
        //int cdsId = 1;
        for( std::vector<GenomeBasicFeature>::iterator it4=it->getCdsVector().begin();
             it4!=it->getCdsVector().end(); ++it4 ){
            ofile << it->getChromeSomeName() << "\t"+transcriptResource+"\tCDS\t" << (*it4).getStart() << "\t" <<
                  (*it4).getEnd() << "\t.\t" <<st << "\t.\tParent=" << it->getName()<< std::endl;
        }
        ofile << "#metainformation: " << it->getMetaInformation() << std::endl;
        ofile << "#genome sequence: " << it->getGeneomeSequence() << std::endl;
        ofile << "#CDS sequence: " << it->getCdsSequence() << std::endl;
        std::string cdsSequence = it->getCdsSequence();
    }
    ofile.close();
}











//for gene
bool overlap2( const Gene & gene, const Range & range ){
    //std::cout << "gene.getStart() " << gene.getStart() << " gene.getEnd() " << gene.getEnd() << " range.getStart() " << range.getStart() << " range.getEnd() " << range.getEnd() << std::endl;
    if(   ( gene.getStart()<=range.getStart() && range.getStart()<=gene.getEnd()) ||
          ( gene.getStart()<=range.getEnd() && range.getEnd()<=gene.getEnd()) ||
          ( range.getStart()<=gene.getStart() && gene.getStart()<=range.getEnd())||
          ( range.getStart()<=gene.getEnd() && gene.getEnd()<=range.getEnd())   ){
        return true;
    }
    return false;
}
void getAlltheOverLappedGenes2( std::map<std::string, std::vector<std::string> > & geneNameMap,
                                std::map<std::string, Gene > & geneHashMap, const std::string & chr,
                                const AlignmentMatch & alignmentMatch, std::vector<Gene*> & overLappedGenes,
                                size_t & startShitfDistance, size_t & endShiftDistance){
    startShitfDistance=0;
    endShiftDistance=0;
    int newStartShitfDistance=0;
    int newEndShiftDistance=0;
    int i=0, j=geneNameMap[chr].size()-1, newi, newj, k=14;
    bool iChanged=true, jChanged=true;
    int interVal = 3;
    while( k>0 && (iChanged || jChanged)  ){
        --k;
        iChanged=false;
        jChanged=false;
        if( geneHashMap[geneNameMap[chr][i]].getEnd() < alignmentMatch.getDatabase().getStart() ){
            newi = i+(j-i)/interVal;
            if( geneHashMap[geneNameMap[chr][newi]].getEnd() < alignmentMatch.getDatabase().getStart() ){
                iChanged=true;
                i=newi;
            }
        }
        if( geneHashMap[geneNameMap[chr][j]].getStart() > alignmentMatch.getDatabase().getEnd()  ){
            newj = j-(j-i)/interVal;
            if( geneHashMap[geneNameMap[chr][newj]].getStart() > alignmentMatch.getDatabase().getEnd() ){
                jChanged=true;
                j=newj;
            }
        }
        if( !iChanged && !jChanged ){
            interVal=2*interVal;
            iChanged= true;
        }
    }
    --i;
    if( i<0 ){
        i=0;
    }
    ++j;
    if( j >= geneNameMap[chr].size() ){
        j = geneNameMap[chr].size()-1;
    }
   // std::cout << "i " << i << " j " << j << std::endl;
    for( ; i<=j; ++i ){
//        std::cout << "line 599 i " << i << " geneNameMap[chr][i] " << geneNameMap[chr][i] << std::endl;
        if( overlap2(geneHashMap[geneNameMap[chr][i]], alignmentMatch.getDatabase()) ){
  //          std::cout << "line 601" << std::endl;
            newStartShitfDistance = startShitfDistance;
            newEndShiftDistance = endShiftDistance;
            if( geneHashMap[geneNameMap[chr][i]].getStart() < alignmentMatch.getDatabase().getStart() ){
                if( startShitfDistance <= (geneHashMap[geneNameMap[chr][i]].getEnd()-geneHashMap[geneNameMap[chr][i]].getStart()) ){
//                    std::cout << "geneHashMap[geneNameMap[chr][i]].getStart(): " << geneHashMap[geneNameMap[chr][i]].getStart() <<
//                    " alignmentMatch.getDatabase().getStart(): " << alignmentMatch.getDatabase().getStart() <<
//                    " geneHashMap[geneNameMap[chr][i]].getEnd(): " << geneHashMap[geneNameMap[chr][i]].getEnd() << std::endl;
                    newStartShitfDistance = (geneHashMap[geneNameMap[chr][i]].getEnd()-geneHashMap[geneNameMap[chr][i]].getStart());
                }
            }
            if( geneHashMap[geneNameMap[chr][i]].getEnd() > alignmentMatch.getDatabase().getEnd() ){
                if( endShiftDistance < (geneHashMap[geneNameMap[chr][i]].getEnd()-geneHashMap[geneNameMap[chr][i]].getStart()) ){
                    //std::cout << "line 612" << std::endl;
                    newEndShiftDistance = (geneHashMap[geneNameMap[chr][i]].getEnd()-geneHashMap[geneNameMap[chr][i]].getStart());
                }
            }
//            std::cout << "newStartShitfDistance: " << newStartShitfDistance << " newEndShiftDistance: " << newEndShiftDistance << std::endl;
            if( newStartShitfDistance <= 100000 && newEndShiftDistance <= 100000 ){
                startShitfDistance=newStartShitfDistance;
                endShiftDistance=newEndShiftDistance;
                overLappedGenes.push_back(&geneHashMap[geneNameMap[chr][i]]);
            }
        }
    }
}

bool overlap( const Gene & gene, const Range & range ){
    if( range.getStart()<=gene.getStart() && gene.getEnd()<=range.getEnd() ){
        return true;
    }
    return false;
}

void getAlltheOverLappedGenes( std::map<std::string, std::vector<std::string> > & geneNameMap,
                               std::map<std::string, Gene > & geneHashMap, const std::string & chr,
                               const AlignmentMatch & alignmentMatch, std::vector<Gene*> & overLappedGenes,
                               size_t & startShitfDistance, size_t & endShiftDistance){
    startShitfDistance=0;
    endShiftDistance=0;
    int i=0, j=geneNameMap[chr].size()-1, newi, newj, k=10;
    bool iChanged=true, jChanged=true;
    int interVal = 3;
    while( k>0 && (iChanged || jChanged)  ){
        --k;
        iChanged=false;
        jChanged=false;
        if( geneHashMap[geneNameMap[chr][i]].getEnd() < alignmentMatch.getDatabase().getStart() ){
            newi = i+(j-i)/interVal;
            if( geneHashMap[geneNameMap[chr][newi]].getEnd() < alignmentMatch.getDatabase().getStart() ){
                iChanged=true;
                i=newi;
            }
        }
        if( geneHashMap[geneNameMap[chr][j]].getStart() > alignmentMatch.getDatabase().getEnd()  ){
            newj = j-(j-i)/interVal;
            if( geneHashMap[geneNameMap[chr][newj]].getStart() > alignmentMatch.getDatabase().getEnd() ){
                jChanged=true;
                j=newj;
            }
        }
        if( !iChanged && !jChanged ){
            interVal=2*interVal;
            iChanged= true;
        }
    }
    for( ; i<=j; ++i ){
        if( overlap(geneHashMap[geneNameMap[chr][i]], alignmentMatch.getDatabase()) ){
            overLappedGenes.push_back(&geneHashMap[geneNameMap[chr][i]]);
        }
    }
}

void slidingWinAlnAndGeneRateAnnotation(AlignmentMatch & alignmentMatch,
                                        std::map<std::string, Fasta> & databaseSequences,
                                        std::map<std::string, Fasta> & querySequences,
                                        std::vector<Gene*> & overLappedGenes, std::map<std::string, Transcript> & transcriptHashMap,
                                        std::map<int, int> & importantPositions,
                                        const int & slidingWindowSize, const int & startShitfDistance,
                                        const int & endShiftDistance, std::map<std::string, std::string>& parameters,
                                        NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix){


    int queryStart=alignmentMatch.getQueryStart();
    int queryEnd= alignmentMatch.getQueryEnd();
    int databaseStart=alignmentMatch.getDatabaseStart();
    int databaseEnd=alignmentMatch.getDatabaseEnd();

//    std::cout << "line 677" << std::endl;
    std::string alignQuerySequence="";
    std::string alignDatabaseSequence="";

    std::string alignQuerySequence1="";
    std::string alignDatabaseSequence1="";

    std::string alignQuerySequence2="";
    std::string alignDatabaseSequence2="";

    std::string querySequence0="";
    std::string databaseSequence0="";

    if(startShitfDistance>0){
        databaseStart = alignmentMatch.getDatabaseStart()-startShitfDistance;
        if( databaseStart < 1 ){
            databaseStart = 1;
        }
        if( alignmentMatch.getDatabaseStart()>1 ) {
            databaseSequence0 = getSubsequence(databaseSequences, alignmentMatch.getDatabaseChr(),
                                               databaseStart, alignmentMatch.getDatabaseStart() - 1);
        }

        if( POSITIVE == alignmentMatch.getQueryStrand() && alignmentMatch.getQueryStart()>1 ){
            queryStart = alignmentMatch.getQueryStart()-startShitfDistance;
            if( queryStart < 1 ){
                queryStart=1;
            }
            querySequence0 = getSubsequence(querySequences, alignmentMatch.getQueryChr(), queryStart,
                                            alignmentMatch.getQueryStart()-1);
        } else if ( NEGATIVE == alignmentMatch.getQueryStrand() && alignmentMatch.getQueryEnd() < querySequences[alignmentMatch.getQueryChr()].getSequence().length() ) {
            queryEnd = alignmentMatch.getQueryEnd()+startShitfDistance;
            if( queryEnd > querySequences[alignmentMatch.getQueryChr()].getSequence().length() ){
                queryEnd = querySequences[alignmentMatch.getQueryChr()].getSequence().length();
            }
            querySequence0 = getSubsequence(querySequences, alignmentMatch.getQueryChr(),
                                            alignmentMatch.getQueryEnd()+1,
                                            queryEnd, alignmentMatch.getQueryStrand());
        }
        int longerLength = querySequence0.length();
        if( longerLength < databaseSequence0.length()){
            longerLength = databaseSequence0.length();
        }
        alignSlidingWindow(querySequence0, databaseSequence0, alignQuerySequence, alignDatabaseSequence, longerLength, parameters, nucleotideCodeSubstitutionMatrix);
    }
    std::string querySequence1 = getSubsequence(querySequences, alignmentMatch.getQueryChr(),
                                               alignmentMatch.getQueryStart(),
                                               alignmentMatch.getQueryEnd(), alignmentMatch.getQueryStrand());
    std::string databaseSequence1 = getSubsequence(databaseSequences, alignmentMatch.getDatabaseChr(),
                                                   alignmentMatch.getDatabaseStart(),
                                                   alignmentMatch.getDatabaseEnd());

    std::string querySequence2="";
    std::string databaseSequence2="";
    if( endShiftDistance>0 ){
        databaseEnd = alignmentMatch.getDatabaseEnd()+endShiftDistance;
        if( databaseEnd > databaseSequences[alignmentMatch.getDatabaseChr()].getSequence().length() ){
            databaseEnd = databaseSequences[alignmentMatch.getDatabaseChr()].getSequence().length();
        }
        if( alignmentMatch.getDatabaseEnd() < databaseSequences[alignmentMatch.getDatabaseChr()].getSequence().length() ) {
            databaseSequence2 = getSubsequence(databaseSequences, alignmentMatch.getDatabaseChr(),
                                               alignmentMatch.getDatabaseEnd()+1, databaseEnd);
        }

        if( POSITIVE == alignmentMatch.getQueryStrand() && alignmentMatch.getQueryEnd() < querySequences[alignmentMatch.getQueryChr()].getSequence().length() ) {
            queryEnd = alignmentMatch.getQueryEnd() + endShiftDistance;
            if (queryEnd > querySequences[alignmentMatch.getQueryChr()].getSequence().length() ) {
                queryEnd = querySequences[alignmentMatch.getQueryChr()].getSequence().length();
            }
            querySequence2 = getSubsequence(querySequences, alignmentMatch.getQueryChr(),
                             alignmentMatch.getQueryEnd() + 1, queryEnd,
                             alignmentMatch.getQueryStrand());
        } else if ( NEGATIVE == alignmentMatch.getQueryStrand() && alignmentMatch.getQueryStart()>1){
            queryStart = alignmentMatch.getQueryStart() - endShiftDistance;
            if (queryStart < 1 ) {
                queryStart = 1;
            }
            querySequence2 = getSubsequence(querySequences, alignmentMatch.getQueryChr(),
                                            queryStart, alignmentMatch.getQueryStart() - 1,
                                            alignmentMatch.getQueryStrand());
        }
        int longerLength = querySequence2.length();
        if( longerLength < databaseSequence2.length()){
            longerLength = databaseSequence2.length();
        }
        alignSlidingWindow(querySequence2, databaseSequence2, alignQuerySequence2, alignDatabaseSequence2, longerLength, parameters, nucleotideCodeSubstitutionMatrix);
    }


    alignSlidingWindow(querySequence1, databaseSequence1, alignQuerySequence1, alignDatabaseSequence1, slidingWindowSize, parameters, nucleotideCodeSubstitutionMatrix);

    alignQuerySequence = alignQuerySequence + alignQuerySequence1;
    alignQuerySequence = alignQuerySequence + alignQuerySequence2;
    alignDatabaseSequence = alignDatabaseSequence + alignDatabaseSequence1;
    alignDatabaseSequence = alignDatabaseSequence + alignDatabaseSequence2;
    //generate annotation according to sequence alignment begin
//    std::cout << "line 661" << std::endl;
    for( Gene * gene : overLappedGenes ){
//        std::cout << gene->getName() << std::endl;
        for (int i=0; i< gene->getTranscriptVector().size(); ++i) {
//            std::cout << gene->getTranscriptVector()[i] << std::endl;
            if( transcriptHashMap[gene->getTranscriptVector()[i]].getCdsVector().size() > 0 ){
                for (int it4 =0; it4< transcriptHashMap[gene->getTranscriptVector()[i]].getCdsVector().size(); ++it4) {
                    GenomeBasicFeature * cds = & transcriptHashMap[gene->getTranscriptVector()[i]].getCdsVector()[it4];
                    importantPositions[(*cds).getStart()] = 0;
                    importantPositions[(*cds).getEnd()] = 0;
                }
            }
            if( transcriptHashMap[gene->getTranscriptVector()[i]].getExonVector().size()>0 ){
                for (int it4=0; it4< transcriptHashMap[gene->getTranscriptVector()[i]].getExonVector().size(); ++it4) {
                    GenomeBasicFeature *exon = & transcriptHashMap[gene->getTranscriptVector()[i]].getExonVector()[it4];
                    importantPositions[(*exon).getStart()] = 0;
                    importantPositions[(*exon).getEnd()] = 0;
                }
            }
            if( transcriptHashMap[gene->getTranscriptVector()[i]].getFivePrimerUtr().size()>0 ){
                for (int it4=0; it4< transcriptHashMap[gene->getTranscriptVector()[i]].getFivePrimerUtr().size(); ++it4) {
                    importantPositions[transcriptHashMap[gene->getTranscriptVector()[i]].getFivePrimerUtr()[it4].getStart()] = 0;
                    importantPositions[transcriptHashMap[gene->getTranscriptVector()[i]].getFivePrimerUtr()[it4].getEnd()] = 0;
                }
            }
            if( transcriptHashMap[gene->getTranscriptVector()[i]].getThreePrimerUtr().size()>0 ){
                for (int it4=0; it4< transcriptHashMap[gene->getTranscriptVector()[i]].getThreePrimerUtr().size(); ++it4) {
                    importantPositions[transcriptHashMap[gene->getTranscriptVector()[i]].getThreePrimerUtr()[it4].getStart()] = 0;
                    importantPositions[transcriptHashMap[gene->getTranscriptVector()[i]].getThreePrimerUtr()[it4].getEnd()] = 0;
                }
            }
        }
    }
//    std::cout << "line 682" << std::endl;
    int databasePosition = databaseStart - 1;
    int newStart;
    int queryPosition = 0;
    for (size_t tp = 0; tp < alignQuerySequence.length(); ++tp) {
        if (alignQuerySequence[tp] != '-') {
            ++queryPosition;
        }
        if (alignDatabaseSequence[tp] != '-') {
            ++databasePosition;
            if( importantPositions.find(databasePosition) != importantPositions.end() ){
                if (POSITIVE == alignmentMatch.getQueryStrand()) {
                    newStart = queryStart + queryPosition - 1;
                } else {
                    newStart = queryEnd - queryPosition + 1;
                }
                if (newStart < 1) {
                    newStart = 1;
                } else if (newStart >
                           querySequences[alignmentMatch.getQueryChr()].getSequence().length()) {
                    newStart = querySequences[alignmentMatch.getQueryChr()].getSequence().length();
                }
                importantPositions[databasePosition] = newStart;
            }
        }
    }
}

void geneGeneStruAlnAndGeneRateAnnotation(Transcript & newTranscript, Transcript * referenceTranscript, Transcript & newTranscript2,
                                          std::map<std::string, Fasta> & databaseSequences,
                                          std::map<std::string, Fasta> & querySequences, AlignmentMatch & alignmentMatch,
                                          std::map<std::string, std::string>& parameters,
                                          NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix,
                                          std::map<int, int> & importantPositions2){
    //std::cout << "ZSDP " << referenceTranscript->getName() << std::endl;
    int i;
    std::string refGenomeSequence = getSubsequence(databaseSequences,
                                                   referenceTranscript->getChromeSomeName(),
                                                   referenceTranscript->getPStart(),
                                                   referenceTranscript->getPEnd(),
                                                   referenceTranscript->getStrand());

    int startCodonPosition = 1;
    int stopCodonPosition = refGenomeSequence.length() - 2;

    size_t transcriptLength =
            newTranscript.getPEnd() - newTranscript.getPStart();
    size_t startTarget;
    size_t endTarget;
    if (newTranscript.getPStart() > transcriptLength) {
        startTarget = newTranscript.getPStart() - transcriptLength;
    } else {
        startTarget = 1;
    }
    if ((newTranscript.getPEnd() + transcriptLength) >
        querySequences[alignmentMatch.getQueryChr()].getSequence().length()) {
        endTarget = querySequences[alignmentMatch.getQueryChr()].getSequence().length();
    } else {
        endTarget = newTranscript.getPEnd() + transcriptLength;
    }

    STRAND thisStrand;
    if( POSITIVE == alignmentMatch.getQueryStrand() ) {
        thisStrand = referenceTranscript->getStrand();
    }else {
        if (referenceTranscript->getStrand() == POSITIVE) {
            thisStrand = NEGATIVE;
        } else {
            thisStrand = POSITIVE;
        }
    }
    std::string dna_b = getSubsequence(querySequences,
                                       alignmentMatch.getQueryChr(),
                                       startTarget, endTarget, thisStrand);
    std::vector<SpliceSitePosition> spliceSitePositions;
    if (referenceTranscript->getStrand() == POSITIVE) {
        if (referenceTranscript->getCdsVector().size() > 1) {
            for (i = 1; i < referenceTranscript->getCdsVector().size(); ++i) {
                SpliceSitePosition spliceSitePosition(
                        referenceTranscript->getCdsVector()[i - 1].getEnd() -
                        referenceTranscript->getPStart() + 2,
                        referenceTranscript->getCdsVector()[i].getStart() -
                        referenceTranscript->getPStart());
                spliceSitePositions.push_back(spliceSitePosition);
            }
        }
    } else {
        if (referenceTranscript->getCdsVector().size() > 1) {
            for (i = 1; i < referenceTranscript->getCdsVector().size(); ++i) {
                SpliceSitePosition spliceSitePosition(
                        referenceTranscript->getPEnd() -
                        referenceTranscript->getCdsVector()[i].getStart() + 2,
                        referenceTranscript->getPEnd() -
                        referenceTranscript->getCdsVector()[i - 1].getEnd());
                spliceSitePositions.push_back(spliceSitePosition);
            }
        }
    }
    AlignTranscript nw(refGenomeSequence, dna_b,
                                                  startCodonPosition,
                                                  stopCodonPosition,
                                                  spliceSitePositions,
                                                  parameters,
                                                  nucleotideCodeSubstitutionMatrix);
    std::string alignQuerySequence = nw.getAlignment_q();
    std::string alignDatabaseSequence = nw.getAlignment_d();
    int queryPosition = 0;
    int newStart;

    for (std::vector<GenomeBasicFeature>::iterator it4 = referenceTranscript->getCdsVector().begin();
         it4 != referenceTranscript->getCdsVector().end(); ++it4) {
        importantPositions2[(*it4).getStart()] = 0;
        importantPositions2[(*it4).getEnd()] = 0;
    }
    if (referenceTranscript->getStrand() == POSITIVE) {
        int databasePosition =referenceTranscript->getPStart()-1;
        for (size_t tp = 0; tp < alignQuerySequence.length(); ++tp) {
            if (alignQuerySequence[tp] != '-') {
                ++queryPosition;
            }
            if (alignDatabaseSequence[tp] != '-') {
                ++databasePosition;
                if (importantPositions2.find(databasePosition ) !=
                    importantPositions2.end()) {
                    if( POSITIVE == alignmentMatch.getQueryStrand() ) {
                        newStart = startTarget + queryPosition - 1;
                    }else{
                        newStart = endTarget - queryPosition + 1;
                    }
                    if (newStart < 1) {
                        newStart = 1;
                    } else if (newStart >
                               querySequences[alignmentMatch.getQueryChr()].getSequence().length()) {
                        newStart = querySequences[alignmentMatch.getQueryChr()].getSequence().length();
                    }
                    importantPositions2[databasePosition] = newStart;
                }
            }
        }
    }else{
        int databasePosition = referenceTranscript->getPEnd() + 1;
        for (size_t tp = 0; tp < alignQuerySequence.length(); ++tp) {
            if (alignQuerySequence[tp] != '-') {
                ++queryPosition;
            }
            if (alignDatabaseSequence[tp] != '-') {
                --databasePosition;
                if (importantPositions2.find(databasePosition) !=
                    importantPositions2.end()) {
                    if( POSITIVE == alignmentMatch.getQueryStrand() ) {
                        newStart = endTarget - queryPosition + 1;
                    }else{
                        newStart = startTarget + queryPosition - 1;
                    }
                    if (newStart < 1) {
                        newStart = 1;
                    } else if (newStart > querySequences[alignmentMatch.getQueryChr()].getSequence().length()) {
                        newStart = querySequences[alignmentMatch.getQueryChr()].getSequence().length();
                    }
                    importantPositions2[databasePosition] = newStart;
                }
            }
        }
    }
    for (size_t i5 = 0; i5 < referenceTranscript->getCdsVector().size(); ++i5) {
        GenomeBasicFeature cds(
                importantPositions2[referenceTranscript->getCdsVector()[i5].getStart()],
                importantPositions2[referenceTranscript->getCdsVector()[i5].getEnd()]);
        cds.setLastColumnInformation(referenceTranscript->getCdsVector()[i5].getLastColumnInformation());
        cds.setType(referenceTranscript->getCdsVector()[i5].getType());
        cds.setCodonFrame(referenceTranscript->getCdsVector()[i5].getCodonFrame());
        newTranscript2.addCds(cds);
    }
}
void updateNewTranscript( Transcript & newTranscript, Transcript * referenceTranscript, std::map<int, int> & importantPositions ){
    for (size_t i5 = 0;
         i5 < referenceTranscript->getExonVector().size(); ++i5) {
        GenomeBasicFeature exon(
                importantPositions[referenceTranscript->getExonVector()[i5].getStart()],
                importantPositions[referenceTranscript->getExonVector()[i5].getEnd()]);
        exon.setLastColumnInformation(
                referenceTranscript->getExonVector()[i5].getLastColumnInformation());
        exon.setType(referenceTranscript->getExonVector()[i5].getType());
        exon.setCodonFrame(referenceTranscript->getExonVector()[i5].getCodonFrame());
        newTranscript.addExon(exon);
    }
    newTranscript.setScore(referenceTranscript->getScore());
    newTranscript.setLastColumnInformation(referenceTranscript->getLastColumnInformation());
    if( referenceTranscript->getFivePrimerUtr().size()>0 ){
        for (int it4=0; it4< referenceTranscript->getFivePrimerUtr().size(); ++it4) {
            GenomeBasicFeature fivePrimerUtr(
                    importantPositions[referenceTranscript->getFivePrimerUtr()[it4].getStart()],
                    importantPositions[referenceTranscript->getFivePrimerUtr()[it4].getEnd()]);
            fivePrimerUtr.setType(referenceTranscript->getFivePrimerUtr()[it4].getType());
            fivePrimerUtr.setLastColumnInformation(referenceTranscript->getFivePrimerUtr()[it4].getLastColumnInformation());
            fivePrimerUtr.setCodonFrame(referenceTranscript->getFivePrimerUtr()[it4].getCodonFrame());
            newTranscript.addFivePrimerUtr(fivePrimerUtr);
        }
    }
    if( referenceTranscript->is_ifThreePrimerUtr() ){
        for (int it4=0; it4< referenceTranscript->getThreePrimerUtr().size(); ++it4) {
            GenomeBasicFeature threePrimerUtr(
                    importantPositions[referenceTranscript->getThreePrimerUtr()[it4].getStart()],
                    importantPositions[referenceTranscript->getThreePrimerUtr()[it4].getEnd()]);
            threePrimerUtr.setType(referenceTranscript->getThreePrimerUtr()[it4].getType());
            threePrimerUtr.setLastColumnInformation(referenceTranscript->getThreePrimerUtr()[it4].getLastColumnInformation());
            threePrimerUtr.setCodonFrame(referenceTranscript->getThreePrimerUtr()[it4].getCodonFrame());
            newTranscript.addThreePrimerUtr(threePrimerUtr);
        }
    }
}

void outputGffRecords(std::vector<NewGffRecord> &newGffRecords,
                      std::map<std::string, std::vector<std::string> > &geneNameMap,
                      std::map<std::string, Gene> &geneHashMap,
                      std::map<std::string, Transcript> &transcriptHashMap, std::ofstream &ofile){
    for( std::vector<NewGffRecord>::iterator it0=newGffRecords.begin(); it0!=newGffRecords.end(); ++it0 ){
        if( it0->transcripts.size() > 0 ){
            for(  std::vector<Transcript>::iterator it=it0->transcripts.begin();
                  it!=it0->transcripts.end(); ++it){
                int thisTranscriptStart=std::numeric_limits<int>::max();
                int thisTranscriptEnd=0;
                if( transcriptHashMap[it->getName()].getExonVector().size() > 0 ) {
                    for (std::vector<GenomeBasicFeature>::iterator it4 = it->getExonVector().begin();
                         it4 != it->getExonVector().end(); ++it4) {
                        if(thisTranscriptStart>it4->getStart()){
                            thisTranscriptStart=it4->getStart();
                        }
                        if(thisTranscriptEnd<it4->getStart()){
                            thisTranscriptEnd=it4->getEnd();
                        }
                    }
                }
                if( transcriptHashMap[it->getName()].getCdsVector().size() > 0 ) {
                    for (std::vector<GenomeBasicFeature>::iterator it4 = it->getCdsVector().begin();
                         it4 != it->getCdsVector().end(); ++it4) {
                        if(thisTranscriptStart>it4->getStart()){
                            thisTranscriptStart=it4->getStart();
                        }
                        if(thisTranscriptEnd<it4->getStart()){
                            thisTranscriptEnd=it4->getEnd();
                        }
                    }
                }
                it->setStart(thisTranscriptStart);
                it->setEnd(thisTranscriptEnd);
            }
            int thisStart = it0->transcripts[0].getStart();
            int thisEnd = it0->transcripts[0].getEnd();
            for(  std::vector<Transcript>::iterator it=it0->transcripts.begin();
                  it!=it0->transcripts.end(); ++it){
                if( thisStart> it->getStart() ){
                    thisStart = it->getStart();
                }
                if( thisEnd> it->getEnd() ){
                    thisEnd = it->getEnd();
                }
            }
            std::string st = "+";
            if( NEGATIVE ==  it0->transcripts[0].getStrand() ){
                st="-";
            }
            std::string transcriptResource = it0->transcripts[0].getSource();
            if( transcriptResource.length()<1 ){
                transcriptResource="LIFTOVER";
            }
            ofile << it0->transcripts[0].getChromeSomeName() << "\t"+transcriptResource+"\tgene\t" << thisStart << "\t" <<
                  thisEnd << "\t.\t"<< st <<"\t.\t"<< geneHashMap[it0->geneName].getLastColumnInformation() << std::endl;
            for( std::vector<Transcript>::iterator it=it0->transcripts.begin();
                  it!=it0->transcripts.end(); ++it ){
                transcriptResource = it->getSource();
                if( transcriptResource.length()<1 ){
                    transcriptResource="LIFTOVER";
                }
                ofile << it->getChromeSomeName() << "\t"+transcriptResource+"\t"<< transcriptHashMap[it->getName()].getType() <<"\t" << it->getStart() << "\t" <<
                      it->getEnd() << "\t" << it->getScore() << "\t"<< st <<"\t.\t"<< it->getLastColumnInformation() << std::endl;
                std::vector<GenomeBasicFeature> GenomeBasicFeatures;
                if( transcriptHashMap[it->getName()].getFivePrimerUtr().size()>0){
                    for (std::vector<GenomeBasicFeature>::iterator it4 = it->getFivePrimerUtr().begin();
                         it4 != it->getFivePrimerUtr().end(); ++it4) {
                        GenomeBasicFeatures.push_back(*it4);
                    }
                }
                if( transcriptHashMap[it->getName()].getCdsVector().size() > 0 ) {
                    for (std::vector<GenomeBasicFeature>::iterator it4 = it->getCdsVector().begin();
                         it4 != it->getCdsVector().end(); ++it4) {
                        GenomeBasicFeatures.push_back(*it4);
                    }
                }
                if( transcriptHashMap[it->getName()].getExonVector().size() > 0 ) {
                    for (std::vector<GenomeBasicFeature>::iterator it4 = it->getExonVector().begin();
                         it4 != it->getExonVector().end(); ++it4) {
                        GenomeBasicFeatures.push_back(*it4);
                    }
                }
                if( transcriptHashMap[it->getName()].getThreePrimerUtr().size()>0){
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
                        ofile << it->getChromeSomeName() << "\t" + transcriptResource + "\t"+(*it4).getType()+"\t" << (*it4).getStart() << "\t"
                              << (*it4).getEnd() << "\t" << it->getScore() << "\t" << st << "\t" << it4->getCodonFrame() << "\t"
                              << it4->getLastColumnInformation() << std::endl;
                    }
                }
                if( transcriptHashMap[it->getName()].getCdsVector().size() > 0 ){
                    ofile << "#metainformation: " << it->getMetaInformation() << std::endl;
                    ofile << "#genome sequence: " << it->getGeneomeSequence() << std::endl;
                    ofile << "#CDS sequence: " << it->getCdsSequence() << std::endl;
                }
                std::string cdsSequence = it->getCdsSequence();
            }
        }
    }
}

void splitCIGAR( std::string & cigarString, std::vector<std::string> & cigarElems) {
    std::regex reg("([0-9]+[MIDNSHPX=])");
    std::smatch match;
    while( regex_search(cigarString, match, reg) ){
        cigarElems.push_back( match[1] );
        cigarString = match.suffix().str();
    }
}

void TransferAllExonWithSpliceAlignmentResult( const std::string & gffFilePath, const std::string & databaseFastaFilePath,
                                      const std::string & queryFastaFilePath, const std::string & samFile,
                                      std::map<std::string, std::string>& parameters, const std::string & outPutFilePath,
                                      const size_t & minIntron,  const int & slidingWindowSize, const size_t & maxLengthForStructureAlignment){
    const bool slowMode=true;
    std::ifstream infile(samFile);
    if( ! infile.good()){
        std::cerr << "error in opening sam file " << samFile << std::endl;
        exit (1);
    }

    std::map<std::string, std::vector<std::string> > geneNameMap;
    std::map<std::string, Gene> geneHashMap;
    std::map<std::string, Transcript> transcriptHashMap;
    readGffFileWithEveryThing(gffFilePath, geneNameMap, geneHashMap, transcriptHashMap);

    //std::cout << " reading gffFile done " << std::endl;
    std::string line;
    char delim = '\t';
    std::vector<std::string> elems;
    std::string databaseChr;
    size_t databaseStart;
    size_t databaseEnd;
    std::string queryChr;
    size_t queryStart;
    size_t queryEnd;
    size_t windowSize = 0;
    std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMap;
    while (std::getline(infile, line)){ // no matter the transcript in on forward strand or reverse strand, it should do not matter
        if( line[0] != '@' ){ //ignore the header
            elems.clear();
            split(line, delim, elems);
            queryStart=stoi(elems[3]);
//            if(databaseStart>0 && transcriptHashMap.find(elems[0])!=transcriptHashMap.end() ){ // ignore those none mapping records
            if( transcriptHashMap.find(elems[0])!=transcriptHashMap.end() ){ // ignore those none mapping records
                //std::cout << "begain to analysis " << line << std::endl;
                databaseChr=transcriptHashMap[elems[0]].getChromeSomeName();
                queryChr=elems[2];

                databaseStart = transcriptHashMap[elems[0]].getPStart();
                //databaseEnd=databaseStart;
                queryEnd=queryStart;
                std::vector<std::string> cigarElems;
                splitCIGAR(elems[5], cigarElems);
//                std::cout << "cigar split done" << std::endl;
                for(int i=0;i<cigarElems.size();++i) {
                    std::string cVal = cigarElems[i];
                    char cLetter = cVal[cVal.length() - 1];
                    if( i == cigarElems.size()-1 && (cLetter == 'H' || cLetter == 'S' ) ){ // ignore the last soft/hard clipping
                        continue;
                    }
//                    std::cout << cVal << std::endl;
//                    std::cout << cLetter << std::endl;
                    int cLen = stoi(cVal.substr(0, cVal.length() - 1));
//                    std::cout << cLen << std::endl;
                    switch (cLetter) {
                        case 'H':
                            //databaseStart += cLen;
                            //databaseEnd = databaseStart;
                            break;
                        case 'S':
                            //databaseStart += cLen;
                            //databaseEnd = databaseStart;
                            break;
                        case 'M':
                            queryEnd += cLen;
                            //databaseEnd += cLen;
                            break;
                        case '=':
                            queryEnd += cLen;
                            //databaseEnd += cLen;
                            break;
                        case 'X':
                            queryEnd += cLen;
                            //databaseEnd += cLen;
                            break;
                        case 'I':
                            //databaseEnd += cLen;
                            break;
                        case 'D':
                            queryEnd += cLen;
                            break;
                        case 'N':
                            queryEnd += cLen;
                            break;
                        case 'P':
                            break;
                        default:
                            std::cout << "current we could not deal with cigar letter " << cLetter << std::endl;
                            break;
                    }
                }
                databaseEnd = transcriptHashMap[elems[0]].getPEnd();
               // std::cout << "cigar parsing done" << std::endl;
//                --databaseEnd;
                --queryEnd;
                if( alignmentMatchsMap.find(databaseChr) == alignmentMatchsMap.end() ){
                    alignmentMatchsMap[databaseChr] = std::vector<AlignmentMatch>();
                }
                int samFlag = stoi(elems[1]);
                if( (0 == samFlag % 32 && transcriptHashMap[elems[0]].getStrand()==POSITIVE)
                    || (0 != samFlag % 32 && transcriptHashMap[elems[0]].getStrand()==NEGATIVE) ){
                    AlignmentMatch alignmentMatch(queryChr, queryStart, queryEnd, POSITIVE, databaseChr, databaseStart, databaseEnd, windowSize);
                    alignmentMatchsMap[databaseChr].push_back(alignmentMatch);
//                    std::cout << "adding POSITIVE " << alignmentMatch.getDatabaseChr() << " " << alignmentMatch.getDatabaseStart() << " " << alignmentMatch.getDatabaseEnd()
//                              << " " << alignmentMatch.getQueryChr() << " " << alignmentMatch.getQueryStart() << " " << alignmentMatch.getQueryEnd() << std::endl;
                }else{
                    AlignmentMatch alignmentMatch(queryChr, queryStart, queryEnd, NEGATIVE, databaseChr, databaseStart, databaseEnd, windowSize);
                    alignmentMatchsMap[databaseChr].push_back(alignmentMatch);

//                    std::cout << "adding NEGATIVE " << alignmentMatch.getDatabaseChr() << " " << alignmentMatch.getDatabaseStart() << " " << alignmentMatch.getDatabaseEnd()
//                              << " " << alignmentMatch.getQueryChr() << " " << alignmentMatch.getQueryStart() << " " << alignmentMatch.getQueryEnd() << std::endl;
                }
            }else{
                std::cout << "could not analysis line: " << line << std::endl;
            }
        }
    }
    TransferAllExonWithNucmerResult( geneNameMap,
                                     geneHashMap, transcriptHashMap, databaseFastaFilePath, queryFastaFilePath, alignmentMatchsMap, parameters,
                                     outPutFilePath, minIntron, slowMode, slidingWindowSize, maxLengthForStructureAlignment);
}

void TransferAllExonWithNucmerResult( const std::string & gffFilePath, const std::string & databaseFastaFilePath,
                                      const std::string & queryFastaFilePath, const std::string & nucmerFilePath,
                                      std::map<std::string, std::string>& parameters, const std::string & outPutFilePath,
                                      const size_t & minIntron , const bool & slowMode, const int & slidingWindowSize, const size_t & maxLengthForStructureAlignment){
    std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMap;
    nucmerRead(nucmerFilePath, alignmentMatchsMap);

    std::map<std::string, std::vector<std::string> > geneNameMap;
    std::map<std::string, Gene> geneHashMap;
    std::map<std::string, Transcript> transcriptHashMap;
    readGffFileWithEveryThing(gffFilePath, geneNameMap, geneHashMap, transcriptHashMap);

    TransferAllExonWithNucmerResult( geneNameMap,
            geneHashMap, transcriptHashMap, databaseFastaFilePath, queryFastaFilePath, alignmentMatchsMap, parameters,
            outPutFilePath, minIntron, slowMode, slidingWindowSize, maxLengthForStructureAlignment);
}

void TransferAllExonWithNucmerResult(  std::map<std::string, std::vector<std::string> > & geneNameMap,
        std::map<std::string, Gene> & geneHashMap, std::map<std::string, Transcript> &transcriptHashMap,
        const std::string & databaseFastaFilePath,
        const std::string & queryFastaFilePath, std::map<std::string, std::vector<AlignmentMatch>> & alignmentMatchsMap,
        std::map<std::string, std::string>& parameters, const std::string & outPutFilePath,
        const size_t & minIntron , const bool & slowMode, const int & slidingWindowSize, const size_t & maxLengthForStructureAlignment){

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);

    std::map<std::string, Fasta> databaseSequences;
    readFastaFile(databaseFastaFilePath, databaseSequences);
    CheckAndUpdateTranscriptsEnds( transcriptHashMap, databaseSequences, nucleotideCodeSubstitutionMatrix, minIntron);

    std::map<std::string, Fasta> querySequences;
    readFastaFile(queryFastaFilePath, querySequences);

    size_t startShitfDistance=1;
    size_t endShiftDistance=1;

    int newStart;
    int newEnd;

    int i;

    std::ofstream ofile;
    ofile.open(outPutFilePath);

    for(std::map<std::string, std::vector<AlignmentMatch>>::iterator it1=alignmentMatchsMap.begin();
        it1!=alignmentMatchsMap.end(); ++it1){
        if( databaseSequences.find(it1->first) != databaseSequences.end() ){
            for( AlignmentMatch alignmentMatch : it1->second ){
//                std::cout << alignmentMatch.getDatabaseChr() << " " << alignmentMatch.getDatabaseStart() << " " << alignmentMatch.getDatabaseEnd()
//                          << " " << alignmentMatch.getQueryChr() << " " << alignmentMatch.getQueryStart() << " " << alignmentMatch.getQueryEnd() << std::endl;
                if( querySequences.find(alignmentMatch.getQueryChr())!=querySequences.end() ){
                    //update algnmentMatch begin
                    std::vector<Gene*> overLappedGenes;
                    if( slowMode ){
                        //std::cout << "line 1317" << std::endl;
                        getAlltheOverLappedGenes2( geneNameMap, geneHashMap, it1->first, alignmentMatch, overLappedGenes,
                                                   startShitfDistance, endShiftDistance);
                    }else{
                        getAlltheOverLappedGenes( geneNameMap, geneHashMap, it1->first, alignmentMatch, overLappedGenes,
                                                  startShitfDistance, endShiftDistance);
                        //std::cout << "line 1323" << std::endl;
                    }
                    //std::exit(0);
                    //std::cout << "line 982 overLappedGenes.size() " << overLappedGenes.size() << std::endl;
                    if( overLappedGenes.size()>0 ) {
                        std::vector<NewGffRecord> newGffRecords;
//                        std::cout << "line 984" << std::endl;
//                        algnmentMatchUpdate( startShitfDistance, endShiftDistance, alignmentMatch,
//                                             databaseSequences, querySequences);

                        std::map<int, int> importantPositions;
//                        std::cout << "line 990" << std::endl;
                        if( alignmentMatch.getWindowSize()>1 ){
                            slidingWinAlnAndGeneRateAnnotation( alignmentMatch, databaseSequences, querySequences,
                                                                overLappedGenes, transcriptHashMap, importantPositions,
                                                                alignmentMatch.getWindowSize(), startShitfDistance, endShiftDistance, parameters, nucleotideCodeSubstitutionMatrix);
                        }else{
                            slidingWinAlnAndGeneRateAnnotation( alignmentMatch, databaseSequences, querySequences,
                                                            overLappedGenes, transcriptHashMap, importantPositions,
                                                            slidingWindowSize, startShitfDistance, endShiftDistance, parameters, nucleotideCodeSubstitutionMatrix);
                        }
                        //std::cout << "genome sequence alignment done" << std::endl;
                        if( POSITIVE == alignmentMatch.getQueryStrand() ) {
                            for( Gene * gene : overLappedGenes ) {
                                NewGffRecord newGffRecord;
                                newGffRecord.geneName=gene->getName();
                                newGffRecord.strand=gene->getStrand();
//                                std::cout << gene->getName() << " " << gene->getTranscriptVector().size() << std::endl;
                                for (i=0; i<gene->getTranscriptVector().size();++i ) {
//                                    std::cout << transcriptHashMap[gene->getTranscriptVector()[i]].getName() << std::endl;
                                    std::string referenceTranscriptId = gene->getTranscriptVector()[i];
//                                    std::cout << referenceTranscriptId << " " << transcriptHashMap[referenceTranscriptId].getStrand() << std::endl;
                                    Transcript * referenceTranscript = & transcriptHashMap[referenceTranscriptId];
                                    Transcript newTranscript(referenceTranscriptId,
                                                             alignmentMatch.getQueryChr(),
                                                             referenceTranscript->getStrand());
//                                    std::cout << "line 1009" << std::endl;
                                    newTranscriptAddCds(transcriptHashMap[referenceTranscriptId], newTranscript, importantPositions, querySequences, nucleotideCodeSubstitutionMatrix, minIntron);
//                                    std::cout << "line 1011" << std::endl;
//                                    if( transcriptHashMap[referenceTranscriptId].getCdsVector().size()>0 ){
//                                        newTranscript.setIfOrfShift(true); //todo for debuging
//                                    }
                                    if (!newTranscript.getIfOrfShift() || (    (referenceTranscript->getPEnd() - referenceTranscript->getPStart()) >
                                                                               maxLengthForStructureAlignment    )) {
                                        updateNewTranscript( newTranscript, referenceTranscript, importantPositions ); // for exon and UTR
                                        newGffRecord.transcripts.push_back(newTranscript);
//                                        std::cout << "line 1060" << std::endl;
                                    } else {
                                        //std::cout << "realign " << referenceTranscript->getName() << std::endl;
                                        Transcript newTranscript2(referenceTranscript->getName(),
                                                                  alignmentMatch.getQueryChr(),
                                                                  referenceTranscript->getStrand());
//                                        std::cout << "line 980" << std::endl;
                                        std::map<int, int> importantPositions2;
                                        geneGeneStruAlnAndGeneRateAnnotation(newTranscript, referenceTranscript, newTranscript2,
                                                                             databaseSequences, querySequences, alignmentMatch,
                                                                             parameters, nucleotideCodeSubstitutionMatrix, importantPositions2);
                                        TranscriptUpdateCdsInformation(newTranscript2, querySequences);
//                                        std::cout << "line 987" << std::endl;
                                        checkOrfState(newTranscript2, querySequences,
                                                      nucleotideCodeSubstitutionMatrix, minIntron);
//                                        std::cout << "line 990" << std::endl;
                                        if (newTranscript2.getIfOrfShift() ) {
                                            updateNewTranscript( newTranscript, referenceTranscript, importantPositions );
                                            newGffRecord.transcripts.push_back(newTranscript);
                                        } else {
//                                            std::cout << "line 995" << std::endl;
                                            newTranscript2.setScore(referenceTranscript->getScore());
                                            //update exon begin
                                            for (size_t i5 = 0;
                                                 i5 < referenceTranscript->getExonVector().size(); ++i5) {
                                                if( importantPositions2.find(referenceTranscript->getExonVector()[i5].getStart()) != importantPositions2.end() ){
                                                    newStart=importantPositions2[referenceTranscript->getExonVector()[i5].getStart()];
                                                }else{
                                                    newStart=importantPositions[referenceTranscript->getExonVector()[i5].getStart()];
                                                }
                                                if( importantPositions2.find(referenceTranscript->getExonVector()[i5].getEnd()) != importantPositions2.end() ){
                                                    newEnd=importantPositions2[referenceTranscript->getExonVector()[i5].getEnd()];
                                                }else{
                                                    newEnd=importantPositions[referenceTranscript->getExonVector()[i5].getEnd()];
                                                }
                                                GenomeBasicFeature exon(newStart, newEnd );
                                                exon.setLastColumnInformation(
                                                        referenceTranscript->getExonVector()[i5].getLastColumnInformation());
                                                exon.setType(referenceTranscript->getExonVector()[i5].getType());
                                                exon.setCodonFrame(referenceTranscript->getExonVector()[i5].getCodonFrame());
                                                newTranscript2.addExon(exon);
                                            }
                                            //update exon end
                                            newTranscript2.setLastColumnInformation(referenceTranscript->getLastColumnInformation());
                                            if( referenceTranscript->getFivePrimerUtr().size()>0 ){
                                                for(int it6=0; it6 < referenceTranscript->getFivePrimerUtr().size(); ++it6){
                                                    newStart=importantPositions[ referenceTranscript->getFivePrimerUtr()[it6].getStart() ];
                                                    newEnd = importantPositions[referenceTranscript->getFivePrimerUtr()[it6].getEnd()];
                                                    if ( referenceTranscript->getStrand() == POSITIVE ) {
                                                        if( importantPositions2.find( referenceTranscript->getFivePrimerUtr()[it6].getStart() ) != importantPositions2.end() ){
                                                            newStart=importantPositions2[ referenceTranscript->getFivePrimerUtr()[it6].getStart() ];
                                                        }
                                                        if( referenceTranscript->getFivePrimerUtr()[it6].getEnd() ==
                                                            referenceTranscript->getCdsVector()[0].getStart()-1 ){
                                                            newEnd=newTranscript2.getCdsVector()[0].getStart()-1;
                                                        }
                                                    }else{
                                                        if( importantPositions2.find( referenceTranscript->getFivePrimerUtr()[it6].getEnd() ) != importantPositions2.end() ){
                                                            newEnd=importantPositions2[ referenceTranscript->getFivePrimerUtr()[it6].getEnd() ];
                                                        }
                                                        if( referenceTranscript->getFivePrimerUtr()[it6].getStart() ==
                                                            referenceTranscript->getCdsVector()[referenceTranscript->getCdsVector().size()-1].getEnd()+1 ){
                                                            newStart=newTranscript2.getCdsVector()[referenceTranscript->getCdsVector().size()-1].getEnd()+1;
                                                        }
                                                    }
                                                    GenomeBasicFeature fivePrimerUtr(newStart, newEnd);

                                                    fivePrimerUtr.setType(referenceTranscript->getFivePrimerUtr()[it6].getType());
                                                    fivePrimerUtr.setLastColumnInformation(referenceTranscript->getFivePrimerUtr()[it6].getLastColumnInformation());
                                                    fivePrimerUtr.setCodonFrame(referenceTranscript->getFivePrimerUtr()[it6].getCodonFrame());
                                                    newTranscript2.addFivePrimerUtr(fivePrimerUtr);
                                                }
                                            }

                                            if( referenceTranscript->getThreePrimerUtr().size() > 0 ){
                                                for(int it6=0; it6 < referenceTranscript->getThreePrimerUtr().size(); ++it6) {
                                                    newStart = importantPositions[referenceTranscript->getThreePrimerUtr()[it6].getStart()];
                                                    newEnd = importantPositions[referenceTranscript->getThreePrimerUtr()[it6].getEnd()];
                                                    if (referenceTranscript->getStrand() == POSITIVE) {
                                                        if (referenceTranscript->getThreePrimerUtr()[it6].getStart() ==
                                                            referenceTranscript->getCdsVector()[
                                                                    referenceTranscript->getCdsVector().size() -
                                                                    1].getEnd() + 1) {
                                                            newStart = newTranscript2.getCdsVector()[
                                                                            newTranscript2.getCdsVector().size() -
                                                                            1].getEnd() + 1;
                                                        }
                                                        if( importantPositions2.find( referenceTranscript->getThreePrimerUtr()[it6].getEnd() ) != importantPositions2.end() ){
                                                            newEnd=importantPositions2[ referenceTranscript->getThreePrimerUtr()[it6].getEnd() ];
                                                        }

                                                    } else {
                                                        if (referenceTranscript->getThreePrimerUtr()[it6].getEnd() ==
                                                            referenceTranscript->getCdsVector()[0].getStart() - 1) {
                                                            newEnd = newTranscript2.getCdsVector()[0].getStart() - 1;
                                                        }
                                                        if( importantPositions2.find( referenceTranscript->getThreePrimerUtr()[it6].getStart() ) != importantPositions2.end() ){
                                                            newStart=importantPositions2[ referenceTranscript->getThreePrimerUtr()[it6].getStart() ];
                                                        }
                                                    }

                                                    GenomeBasicFeature threePrimerUtr(newStart, newEnd);

                                                    threePrimerUtr.setLastColumnInformation(
                                                            referenceTranscript->getThreePrimerUtr()[it6].getLastColumnInformation());
                                                    threePrimerUtr.setType(referenceTranscript->getThreePrimerUtr()[it6].getType());
                                                    threePrimerUtr.setCodonFrame(referenceTranscript->getThreePrimerUtr()[it6].getCodonFrame());
                                                    newTranscript2.addThreePrimerUtr(threePrimerUtr);
                                                }
                                            }
                                            newTranscript2.setSource("REALIGNMENT");
                                            newGffRecord.transcripts.push_back(newTranscript2);
                                        }
                                    }
                                }
                                newGffRecords.push_back(newGffRecord);
                            }
                        } else {
                            for( Gene * gene : overLappedGenes ) {
                                NewGffRecord newGffRecord;
                                newGffRecord.geneName = gene->getName();
                                for (i=0; i<gene->getTranscriptVector().size();++i ) {
                                    std::string referenceTranscriptId = gene->getTranscriptVector()[i];
                                    Transcript * referenceTranscript = & transcriptHashMap[referenceTranscriptId];
                                    STRAND thisStrand;
                                    if (referenceTranscript->getStrand() == POSITIVE) {
                                        thisStrand = NEGATIVE;
                                    } else {
                                        thisStrand = POSITIVE;
                                    }
                                    newGffRecord.strand=thisStrand;
                                    Transcript newTranscript(referenceTranscript->getName(), alignmentMatch.getQueryChr(), thisStrand);
                                    newTranscriptAddCds((*referenceTranscript), newTranscript, importantPositions, querySequences, nucleotideCodeSubstitutionMatrix, minIntron);
//                                    if( transcriptHashMap[referenceTranscriptId].getCdsVector().size()>0 ){
//                                        newTranscript.setIfOrfShift(true); //todo for debuging
//                                    }
                                    if (!newTranscript.getIfOrfShift() || (    (referenceTranscript->getPEnd() - referenceTranscript->getPStart()) >
                                                                               maxLengthForStructureAlignment    )) {
//                                        std::cout << "line 1038" << std::endl;
                                        updateNewTranscript( newTranscript, referenceTranscript, importantPositions );
//                                        std::cout << "line 1040" << std::endl;
                                        newGffRecord.transcripts.push_back(newTranscript);
                                    } else {
//                                        std::cout << "realign " << referenceTranscript->getName() << std::endl;
                                        Transcript newTranscript2(referenceTranscript->getName(),
                                                                  alignmentMatch.getQueryChr(),
                                                                  thisStrand);

                                        std::map<int, int> importantPositions2;
//                                        std::cout << "line 1092" << std::endl;
                                        geneGeneStruAlnAndGeneRateAnnotation(newTranscript, referenceTranscript, newTranscript2,
                                                                             databaseSequences, querySequences, alignmentMatch,
                                                                             parameters, nucleotideCodeSubstitutionMatrix, importantPositions2);
//                                        std::cout << "line 1096" << std::endl;
                                        TranscriptUpdateCdsInformation(newTranscript2, querySequences);
//                                        std::cout << "line 1098" << std::endl;
                                        checkOrfState(newTranscript2, querySequences, nucleotideCodeSubstitutionMatrix, minIntron);
                                        if (newTranscript2.getIfOrfShift()) {
                                            updateNewTranscript( newTranscript, referenceTranscript, importantPositions );
                                            newGffRecord.transcripts.push_back(newTranscript);
                                        } else {
//                                            std::cout << "line 1103" << std::endl;
                                            newTranscript2.setSource("REALIGNMENT");
                                            newTranscript2.setScore(referenceTranscript->getScore());
                                            //update exon begin
                                            for (size_t i5 = 0;
                                                 i5 < referenceTranscript->getExonVector().size(); ++i5) {
                                                if( importantPositions2.find(referenceTranscript->getExonVector()[i5].getStart()) != importantPositions2.end() ){
                                                    newStart=importantPositions2[referenceTranscript->getExonVector()[i5].getStart()];
                                                }else{
                                                    newStart=importantPositions[referenceTranscript->getExonVector()[i5].getStart()];
                                                }
                                                if( importantPositions2.find(referenceTranscript->getExonVector()[i5].getEnd()) != importantPositions2.end() ){
                                                    newEnd=importantPositions2[referenceTranscript->getExonVector()[i5].getEnd()];
                                                }else{
                                                    newEnd=importantPositions[referenceTranscript->getExonVector()[i5].getEnd()];
                                                }
                                                GenomeBasicFeature exon(newStart, newEnd );
                                                exon.setLastColumnInformation(
                                                        referenceTranscript->getExonVector()[i5].getLastColumnInformation());
                                                exon.setType(referenceTranscript->getExonVector()[i5].getType());
                                                exon.setCodonFrame(referenceTranscript->getExonVector()[i5].getCodonFrame());
                                                newTranscript2.addExon(exon);
                                            }
                                            //update exon end
                                            newTranscript2.setScore(referenceTranscript->getScore());
                                            newTranscript2.setLastColumnInformation(referenceTranscript->getLastColumnInformation());

                                            if(referenceTranscript->getFivePrimerUtr().size()>0){ //todo make sure this is correct
                                                for(int it6=0; it6 < referenceTranscript->getFivePrimerUtr().size(); ++it6) {
                                                    newStart = importantPositions[referenceTranscript->getFivePrimerUtr()[it6].getEnd()];
                                                    newEnd = importantPositions[referenceTranscript->getFivePrimerUtr()[it6].getStart()];
                                                    if (referenceTranscript->getStrand() == POSITIVE) {
                                                        if (referenceTranscript->getFivePrimerUtr()[it6].getEnd() ==
                                                            referenceTranscript->getCdsVector()[0].getStart() - 1) {
                                                            newStart=
                                                                    newTranscript2.getCdsVector()[newTranscript2.getCdsVector().size() - 1].getEnd() + 1;
                                                        }
                                                        if( importantPositions2.find(referenceTranscript->getFivePrimerUtr()[it6].getStart()) != importantPositions2.end() ){
                                                            newEnd=importantPositions2[referenceTranscript->getFivePrimerUtr()[it6].getStart()];
                                                        }
                                                    } else {
                                                        if (referenceTranscript->getFivePrimerUtr()[it6].getStart() ==
                                                            referenceTranscript->getCdsVector()[
                                                                    referenceTranscript->getCdsVector().size() - 1].getEnd() + 1) {
                                                            newEnd= newTranscript2.getCdsVector()[0].getStart() - 1;
                                                        }
                                                        if( importantPositions2.find(referenceTranscript->getFivePrimerUtr()[it6].getEnd()) != importantPositions2.end() ){
                                                            newStart=importantPositions2[referenceTranscript->getFivePrimerUtr()[it6].getEnd()];
                                                        }
                                                    }
                                                    GenomeBasicFeature fivePrimerUtr(newStart, newEnd);
                                                    fivePrimerUtr.setLastColumnInformation(
                                                            referenceTranscript->getFivePrimerUtr()[it6].getLastColumnInformation());
                                                    fivePrimerUtr.setType(
                                                            referenceTranscript->getFivePrimerUtr()[it6].getType());
                                                    fivePrimerUtr.setCodonFrame(referenceTranscript->getFivePrimerUtr()[it6].getCodonFrame());
                                                    newTranscript2.addFivePrimerUtr(fivePrimerUtr);
                                                }
                                            }
                                            if(referenceTranscript->getThreePrimerUtr().size()>0) {
                                                for(int it6=0; it6 < referenceTranscript->getThreePrimerUtr().size(); ++it6) {
                                                    newStart = importantPositions[referenceTranscript->getThreePrimerUtr()[it6].getEnd()];
                                                    newEnd = importantPositions[referenceTranscript->getThreePrimerUtr()[it6].getStart()];

                                                    // update utr start
                                                    if (referenceTranscript->getStrand() == POSITIVE) {
                                                        if (referenceTranscript->getThreePrimerUtr()[it6].getStart() ==
                                                            referenceTranscript->getCdsVector()[
                                                                    referenceTranscript->getCdsVector().size() -
                                                                    1].getEnd() + 1) {
                                                            newEnd = newTranscript2.getCdsVector()[0].getStart() - 1;
                                                        }
                                                        if( importantPositions2.find(referenceTranscript->getThreePrimerUtr()[it6].getEnd() ) != importantPositions2.end() ){
                                                            newStart=importantPositions2[referenceTranscript->getThreePrimerUtr()[it6].getEnd()];
                                                        }

                                                    } else {
                                                        if (referenceTranscript->getThreePrimerUtr()[it6].getEnd() ==
                                                            referenceTranscript->getCdsVector()[0].getStart() - 1) {
                                                            newStart = newTranscript2.getCdsVector()[
                                                                            newTranscript2.getCdsVector().size() -
                                                                            1].getEnd() + 1;
                                                        }
                                                        if( importantPositions2.find(referenceTranscript->getThreePrimerUtr()[it6].getStart() ) != importantPositions2.end() ){
                                                            newEnd=importantPositions2[referenceTranscript->getThreePrimerUtr()[it6].getStart()];
                                                        }
                                                    }// update utr begin
                                                    GenomeBasicFeature threePrimerUtr(newStart, newEnd);
                                                    threePrimerUtr.setCodonFrame(referenceTranscript->getThreePrimerUtr()[it6].getCodonFrame());
                                                    threePrimerUtr.setLastColumnInformation(
                                                            referenceTranscript->getThreePrimerUtr()[it6].getLastColumnInformation());
                                                    threePrimerUtr.setType(
                                                            referenceTranscript->getThreePrimerUtr()[it6].getType());
                                                }
                                            }
                                            newGffRecord.transcripts.push_back(newTranscript2);
                                        }
                                    }
                                }//generate annotation according to sequence alignment end
                                newGffRecords.push_back(newGffRecord);
                            }
                        }
                        outputGffRecords(newGffRecords, geneNameMap, geneHashMap, transcriptHashMap, ofile);
                    }
                }else{
                    std::cout << "could not find" << alignmentMatch.getQueryChr() << " in the query genome sequence" << std::endl;
                }
            }
        }else{
            std::cout << "could not find" << it1->first << " in the database genome sequence" << std::endl;
        }
    }
    ofile.close();
}
