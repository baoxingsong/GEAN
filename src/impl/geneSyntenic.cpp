//
// Created by Baoxing Song on 2019-03-13.
//

#include "geneSyntenic.h"

/**
 * this function try to keep those genes in the syntenic region using a longest path algorithm
 * which is a kind of global alignment method
 *
 * longest path algorithm from here:
 * https://www.geeksforgeeks.org/find-longest-path-directed-acyclic-graph/
 *
 */
// this function put the forward entries in the increasing order and put the reversion entries in a decrease order
void myOrthologPairsSort( std::vector<OrthologPair> & pairedSimilarFragments, const double & score, const double & penalty, const double & scoreThreshold, const bool & keepTandemDuplication){
    std::sort(pairedSimilarFragments.begin(), pairedSimilarFragments.end(), [](OrthologPair a, OrthologPair b) {
        return a < b;
    });
    // the following part is for reversion
    int startIndex=0;
    int endIndex=0;
    double maxScore=0;
    double currentScore=0;
    for( int idx=0;  idx<pairedSimilarFragments.size(); ++idx){
        if( NEGATIVE == pairedSimilarFragments[idx].getStrand() ){  // reverse strand  we need to do something special for inversion
            // look for all following pairs that are reverse strand,
            // ie get all reverse strand entries in this group of reverse strands
            for( int jdx=idx; jdx< pairedSimilarFragments.size(); ++jdx ){
                if( NEGATIVE == pairedSimilarFragments[jdx].getStrand() ) {
                    if( idx == jdx ){ // the first one.
                        currentScore+=score*pairedSimilarFragments[jdx].getScore();
                    }else{ // for the reverse alignments, check the strand of the previous alignment
                        // If both current and previous are reverse strand, and
                        // if previous assembly start is greater than current assembly start, increase score.
                        //  else, apply penalty (they are out of order)
                        if ( pairedSimilarFragments[jdx-1].getQueryMiddlePos() > pairedSimilarFragments[jdx].getQueryMiddlePos() ){ // GOOD INVERSION
                            currentScore+=score*pairedSimilarFragments[jdx].getScore();
                        }else if(keepTandemDuplication && pairedSimilarFragments[jdx-1].getRefMiddlePos() == pairedSimilarFragments[jdx].getRefMiddlePos() ){ // tandem duplication
                            currentScore+=score*pairedSimilarFragments[jdx].getScore();
                        }else{
                            currentScore+=penalty*pairedSimilarFragments[jdx].getScore(); // GIVE PENALTY
                        }
                    }
                }else{
                    currentScore+=penalty*pairedSimilarFragments[jdx].getScore(); // penalty because are now forward strand
                }
                if (maxScore < currentScore) {
                    maxScore = currentScore;
                    endIndex=jdx; // keeps track of where to stop the reverse strand grouping
                }
                // If score is negative, stop the loop.  This will happen as we find more
                // forward vs reverse strands.  If there was just 1 reverse strand followed
                // by a forward strand, currentScore goes from 3 to -1 (3 plus -4)
                if( currentScore<0 ){
                    break;
                }
            }
            // if maxScore is larger than scoreThreshold, it means there were multiple reverse
            // strand entries in this group, or a single reverse alignment of significant length.
            // If we have several reverse alignments, we think it is real.  If just 1, may be false
            // alignment.  If the maxScore is greater than the scoreThreshold, we treat it as real.
            // Flip all elements in the range so assembly coordinates are in increasing order.
            if( maxScore>scoreThreshold ){
                maxScore = 0.0 ;
                currentScore = 0.0;
                // loop to find start index of elements we want to flip
                for( int jdx=endIndex; jdx>=idx; --jdx ){
                    if( NEGATIVE == pairedSimilarFragments[jdx].getStrand() ) {
                        if( jdx>idx ){
                            // Verify the reverse alignments are in order to each other.  If not,
                            // apply penalty.  This doesn't prevent overlaps, which will be dealt with later.
                            if ( pairedSimilarFragments[jdx-1].getQueryMiddlePos() > pairedSimilarFragments[jdx].getQueryMiddlePos() ){
                                currentScore+=score*pairedSimilarFragments[jdx].getScore();
                            }else if( keepTandemDuplication && pairedSimilarFragments[jdx-1].getRefMiddlePos() == pairedSimilarFragments[jdx].getRefMiddlePos() ){ // tandem duplication
                                currentScore+=score*pairedSimilarFragments[jdx].getScore();
                            }else{
                                currentScore+=penalty*pairedSimilarFragments[jdx].getScore(); // GIVE PENALTY
                            }
                        }else{
                            currentScore+=score*pairedSimilarFragments[jdx].getScore();
                        }
                    }else{
                        currentScore+=penalty*pairedSimilarFragments[jdx].getScore(); // GIVE PENALTY
                    }
                    if (maxScore < currentScore) {
                        maxScore = currentScore;
                        startIndex=jdx;
                    }
                    if( currentScore<0 ){
                        break;
                    }
                }

                int length = (endIndex-startIndex+1)/2;
                if(length > 0) {
                    for (int j = 0; j < length; ++j) {
                        OrthologPair temp = pairedSimilarFragments[startIndex + j];
                        pairedSimilarFragments[startIndex + j]=pairedSimilarFragments[endIndex - j];
                        pairedSimilarFragments[endIndex - j] = temp;
                    }
                    // Flip the elements in the list so the asm reverse alignments
                    // all have increasing asm values. (ie the last element of the reverse
                    // grouping is now the first, and the first is the last.
                    bool thereAreReverseAlignments = true;
                    while (thereAreReverseAlignments) {
                        thereAreReverseAlignments = false;
                        for (int j = 1; j < length; ++j) {
                            // If a single reference position has multiple assembly alignments mapping to it,
                            // swap the order until the assembly positions are all increasing.
                            if (pairedSimilarFragments[startIndex + j - 1].getRefMiddlePos() == pairedSimilarFragments[startIndex + j].getRefMiddlePos() &&
                                pairedSimilarFragments[startIndex + j - 1].getQueryMiddlePos() > pairedSimilarFragments[startIndex + j].getQueryMiddlePos() ) {
                                thereAreReverseAlignments = true;
                                OrthologPair temp = pairedSimilarFragments[startIndex + j];
                                pairedSimilarFragments[startIndex + j]=pairedSimilarFragments[startIndex + j - 1];
                                pairedSimilarFragments[startIndex + j - 1]=temp;
                            }
                        }
                    }
                }
                // Flip the elements in the list so the asm reverse ailgnments
                // all have increasing asm values. (ie the last element of the reverse
                // grouping is now the first, and the first is the last.

                idx=endIndex;
            }
            maxScore=0.0;
            currentScore=0.0;
        }
    }
}



void longestPath (std::vector<OrthologPair> & pairedSimilarFragments, std::vector<OrthologPair> & sortedOrthologPairs, const bool & keepTandemDuplication){
    double maxSore = 0;
    int bestEnd = 0;
    double scoreArray [pairedSimilarFragments.size()]; // arrays of scores
    int prev [pairedSimilarFragments.size()];  // index of previous node in longest path
    scoreArray[0] = pairedSimilarFragments[0].getScore();
    prev[0] = -1;
    for (int idx = 1; idx < pairedSimilarFragments.size(); ++idx) {
        scoreArray[idx] = pairedSimilarFragments[idx].getScore();
        prev[idx] = -1;
        for (int jdx = idx - 1; jdx >= 0; --jdx) {// checking all previous nodes
            // Because we swapped asm/query start position so that inversions were all increasing,
            // we should always be on the diagonal.  If not, then we filter it.
            // This gets rid of the noise, while preserving the inversions on
            // the diagonal
            // Are only looking at positions previous to our current "idx" position
            if ( (scoreArray[jdx] + pairedSimilarFragments[idx].getScore()) > scoreArray[idx] &&
                pairedSimilarFragments[jdx].getQueryMiddlePos() < pairedSimilarFragments[idx].getQueryMiddlePos()){
                scoreArray[idx] = scoreArray[jdx] + pairedSimilarFragments[idx].getScore();
                prev[idx] = jdx;
            }else if ( (scoreArray[jdx] + pairedSimilarFragments[idx].getScore()) > scoreArray[idx] &&
                      pairedSimilarFragments[jdx].getQueryMiddlePos() == pairedSimilarFragments[idx].getQueryMiddlePos()
                      && keepTandemDuplication) {
                scoreArray[idx] = scoreArray[jdx] + pairedSimilarFragments[idx].getScore();
                prev[idx] = jdx;
            }
        }
        if (scoreArray[idx] > maxSore){
            bestEnd = idx;
            maxSore = scoreArray[idx];
        }
    }
    int idx=bestEnd; // bestEnd is where to stop the longest path
    sortedOrthologPairs.push_back(pairedSimilarFragments[idx]);
    int jdx = prev[idx]; // prev[] is index on the longest path
    while( jdx>=0 ){
        sortedOrthologPairs.push_back(pairedSimilarFragments[jdx]);
        jdx=prev[jdx];
    }
    // Reversing the order
    std::reverse(std::begin(sortedOrthologPairs), std::end(sortedOrthologPairs));
}


void outputGffRecords(Gene & gene, std::ofstream & ofile, std::set<std::string> & geneNames, std::set<std::string> & transcriptNames  ){
    std::string newGeneName=gene.getName();
    if( geneNames.find(newGeneName) != geneNames.end() ){
        int i=1;
        newGeneName = newGeneName + "_" + std::to_string(i);
        while( geneNames.find(newGeneName) != geneNames.end() ){
            ++i;
            newGeneName = gene.getName() + "_" + std::to_string(i);
        }
        geneNames.insert(newGeneName);
    }else{
        geneNames.insert(newGeneName);
    }
    if( gene.getTranscripts().size() > 0 ){
        // begin of gene output
        // begin update transcript positions
        for(  std::vector<Transcript>::iterator it=gene.getTranscripts().begin();
              it!=gene.getTranscripts().end(); ++it){
            int thisTranscriptStart=std::numeric_limits<int>::max();
            int thisTranscriptEnd=0;
            if( it->getExonVector().size() > 0 ) {
                for (std::vector<GenomeBasicFeature>::iterator it4 = it->getExonVector().begin();
                     it4 != it->getExonVector().end(); ++it4) {
                    if(thisTranscriptStart>it4->getStart()){
                        thisTranscriptStart=it4->getStart();
                    }
                    if(thisTranscriptEnd<it4->getEnd()){
                        thisTranscriptEnd=it4->getEnd();
                    }
                }
            }
            if( it->getCdsVector().size() > 0 ) {
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
        // end of updating transcript positions
        int thisStart = gene.getTranscripts()[0].getStart();
        int thisEnd = gene.getTranscripts()[0].getEnd();
        for(  std::vector<Transcript>::iterator it=gene.getTranscripts().begin();
              it!=gene.getTranscripts().end(); ++it){
            if( thisStart> it->getStart() ){
                thisStart = it->getStart();
            }
            if( thisEnd< it->getEnd() ){
                thisEnd = it->getEnd();
            }
        }
        std::string st = "+";
        if( NEGATIVE == gene.getTranscripts()[0].getStrand() ){
            st="-";
        }
        std::string transcriptResource = gene.getTranscripts()[0].getSource();
        if( transcriptResource.length()<1 ){
            transcriptResource="LIFTOVER";
        }
        ofile << gene.getTranscripts()[0].getChromeSomeName() << "\t"+transcriptResource+"\tgene\t" << thisStart << "\t" <<
              thisEnd << "\t.\t"<< st <<"\t.\tID=" << newGeneName<<";" << std::endl;
        // end of gene output

        for( std::vector<Transcript>::iterator it=gene.getTranscripts().begin();
             it!=gene.getTranscripts().end(); ++it ){
            transcriptResource = it->getSource();
            if( transcriptResource.length()<1 ){
                transcriptResource="LIFTOVER";
            }
            std::string newTranscriptName=it->getName();
            if( transcriptNames.find(newTranscriptName) != transcriptNames.end() ){
                int i=1;
                newTranscriptName = newTranscriptName + "_" + std::to_string(i);
                while( transcriptNames.find(newTranscriptName) != transcriptNames.end() ){
                    ++i;
                    newTranscriptName = it->getName() + "_" + std::to_string(i);
                }
                transcriptNames.insert(newTranscriptName);
            }else{
                transcriptNames.insert(newTranscriptName);
            }
            ofile << it->getChromeSomeName() << "\t"+transcriptResource+"\t"<< it->getType() <<"\t" << it->getStart() << "\t" <<
                  it->getEnd() << "\t" << it->getScore() << "\t"<< st <<"\t.\tID="<< newTranscriptName<<";Parent=" << newGeneName<<";" << std::endl;
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
                    ofile << it->getChromeSomeName() << "\t" + transcriptResource + "\t"+(*it4).getType()+"\t" << (*it4).getStart() << "\t"
                          << (*it4).getEnd() << "\t" << it->getScore() << "\t" << st << "\t" << it4->getCodonFrame() << "\tParent="
                          << newTranscriptName << ";" << std::endl;
                }
            }
            if( it->getCdsVector().size() > 0 ){
                ofile << "#metainformation: " << it->getMetaInformation() << std::endl;
                ofile << "#genome sequence: " << it->getGeneomeSequence() << std::endl;
                ofile << "#CDS sequence: " << it->getCdsSequence() << std::endl;
            }
            std::string cdsSequence = it->getCdsSequence();
        }
    }
}

// in this equation the score and penalty is for the purpose of reversion
void generateLongestOutput( const std::string & referenceGffFile, const std::string & queryNewGffFile,
                                const std::string & queryGenomeFile, const std::string & outputGffFile, const int & minIntron,
                                double score, double penalty, double scoreThreshold, const bool & keepTandemDuplication,
                                std::map<std::string, std::string>& parameters, const double & syntenicScore,
                                const double & orfScore, const double & dropLengthThredshold, const bool & onlySyntenic ){
//    std::cout << "longest line 308" << std::endl;
    std::map<std::string, Fasta> queryGenome;
    readFastaFile(queryGenomeFile, queryGenome);
//    std::cout << "longest line 311" << std::endl;
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
//    std::cout << "longest line 313" << std::endl;
    std::map<std::string, std::vector<Gene> > queryGenes;
//    std::cout << "longest line 315" << std::endl;
    readGffFileWithEveryThing (queryNewGffFile, queryGenes);
//    std::cout << "longest line 317" << std::endl;
    std::map<std::string, std::vector<std::string> > referenceGeneNameMap;
    std::map<std::string, Gene > referenceGeneHashMap;
    std::map<std::string, Transcript> referenceTranscriptHashMap;
//    std::cout << "longest line 321" << std::endl;
    readGffFileWithEveryThing ( referenceGffFile, referenceGeneNameMap, referenceGeneHashMap, referenceTranscriptHashMap);
//    std::cout << "longest line 323" << std::endl;
    std::ofstream ofile;
    ofile.open(outputGffFile);

//    std::cout << "longest line 328" << std::endl;

//filter identical records from query begin
    for( std::map<std::string, std::vector<Gene> >::iterator it=queryGenes.begin(); it!=queryGenes.end(); ++it){
        int i, j, n;
        n=it->second.size();
        for ( i = j = 0;  i < n;  ++i ) {
            if ( it->second[i].getName()==it->second[j].getName() && it->second[i].getStart()==it->second[j].getStart()
                 && it->second[i].getEnd()==it->second[j].getEnd() ) {

            }else{
                ++j;
                if  (i != j){
                    it->second[j] = it->second[i];
                }
            }
        }
        ++j;
        it->second.resize(j);
        queryGenes[it->first]=it->second;
    }
    //filter identical records from query end




    std::map<std::string, std::vector<Gene>> keepGenes;
    std::map<std::string, double> geneScores;

    for(std::map<std::string, std::vector<Gene>>::iterator it1=queryGenes.begin(); it1!=queryGenes.end(); it1++){
        std::vector<OrthologPair> orthologPairs;
        if( referenceGeneNameMap.find(it1->first) != referenceGeneNameMap.end() && queryGenome.find(it1->first)!=queryGenome.end() ){
            for( int i=0; i< it1->second.size(); ++i ){
                if( referenceGeneHashMap.find(it1->second[i].getName())!=referenceGeneHashMap.end()
                    && referenceGeneHashMap[it1->second[i].getName()].getChromeSomeName().compare(it1->second[i].getChromeSomeName())==0  ){
                    // this is a chromosome to chromosome somehow global alignment, so only the same chromosome is aligned
                    int targetIndex=i;
                    uint32_t refStartPos=referenceGeneHashMap[it1->second[i].getName()].getStart();
                    uint32_t refEndPos=referenceGeneHashMap[it1->second[i].getName()].getEnd();
                    uint32_t queryStartPos=it1->second[i].getStart();
                    uint32_t queryEndPos=it1->second[i].getEnd();
                    double thisScore=0;

                    for(  int index=0; index<it1->second[i].getTranscripts().size(); ++index ){
                        Transcript transcript=it1->second[i].getTranscripts()[index];
                        TranscriptUpdateCdsInformation(transcript, queryGenome);
                        checkOrfState(transcript, queryGenome, nucleotideCodeSubstitutionMatrix, minIntron);
                        if( transcript.getIfOrfShift() ){
                            // this score is for gene, if there is one ORF conserved transcript, then give this gene a positive score
                        }else{
                            thisScore=orfScore;
                        }
                        it1->second[i].getTranscripts()[index]=transcript;
                        //break;
                    }

                    double lengthRatio;
                    if( it1->second[i].getEnd() == it1->second[i].getStart() ){
                        lengthRatio=0.0;
                    }else{
                        lengthRatio=(double)((double)it1->second[i].getEnd()-(double)it1->second[i].getStart())/((double)referenceGeneHashMap[it1->second[i].getName()].getEnd()-(double)referenceGeneHashMap[it1->second[i].getName()].getStart());
                    }
                    if( lengthRatio > 1 ){
                        lengthRatio=1/lengthRatio;
                    }
                    thisScore += lengthRatio;
                    STRAND strand;
                    if( it1->second[i].getStrand() == referenceGeneHashMap[it1->second[i].getName()].getStrand() ){
                        strand=POSITIVE;
                    }else{
                        strand=NEGATIVE;
                    }
                    if( lengthRatio > dropLengthThredshold ){
                        OrthologPair orthologPair( targetIndex, refStartPos, refEndPos, queryStartPos, queryEndPos, thisScore, strand );
                        orthologPairs.push_back(orthologPair);
                    }
                }
            }
            myOrthologPairsSort( orthologPairs, score, penalty,  scoreThreshold, keepTandemDuplication);
            std::vector<OrthologPair> sortedOrthologPairs;
            longestPath (orthologPairs, sortedOrthologPairs, keepTandemDuplication);
            for( OrthologPair orthologPair : sortedOrthologPairs ){
                //outputGffRecords(it1->second[orthologPair.getQueryIndex()], ofile, geneNames, transcriptNames  );
                Gene g = it1->second[orthologPair.getQueryIndex()];
                double thisScore = orthologPair.getScore() + syntenicScore;
                if( geneScores.find(g.getName()) != geneScores.end() ){
                    if( geneScores[g.getName()] < thisScore ){
                        geneScores[g.getName()] = thisScore;
                        std::vector<Gene> gs;
                        gs.push_back(g);
                        keepGenes[g.getName()] = gs;
                    }else if( keepTandemDuplication && geneScores[g.getName()] == thisScore ){
                        geneScores[g.getName()] = thisScore;
                        keepGenes[g.getName()].push_back(g);
                    }
                }else{
                    geneScores[g.getName()] = thisScore;
                    std::vector<Gene> gs;
                    gs.push_back(g);
                    keepGenes[g.getName()] = gs;
                }
            }
        }
    }
    if(onlySyntenic){

    }else{
        //std::cout << " line 385" << std::endl;
        for(std::map<std::string,std::vector<Gene>>::iterator it1=queryGenes.begin(); it1!=queryGenes.end(); it1++){
            for( Gene g : it1->second ){
                //if( geneNames.find(g.getName())!=geneNames.end() ){

                //}else{
                    double thisScore=0;
                    for(  int index=0; index<g.getTranscripts().size(); ++index ){
                        Transcript transcript=g.getTranscripts()[index];
                        TranscriptUpdateCdsInformation(transcript, queryGenome);
                        checkOrfState(transcript, queryGenome, nucleotideCodeSubstitutionMatrix, minIntron);
                        if( transcript.getIfOrfShift() ){
                            // this score is for gene, if there is one ORF conserved transcript, then give this gene a positive score
                        }else{
                            thisScore=orfScore;
                        }
                        g.getTranscripts()[index]=transcript;
                    }

                    double lengthRatio;
                    if( g.getEnd() == g.getStart() ){
                        lengthRatio=0.0;
                    }else{
                        lengthRatio=(double)((double)g.getEnd()-(double)g.getStart())/((double)referenceGeneHashMap[g.getName()].getEnd()-(double)referenceGeneHashMap[g.getName()].getStart());
                    }
                    if( lengthRatio > 1 ){
                        lengthRatio=1/lengthRatio;
                    }
                    thisScore += lengthRatio;
                    if( geneScores.find(g.getName()) != geneScores.end() ){
                        if( geneScores[g.getName()] < thisScore ){
                            geneScores[g.getName()] = thisScore;
                            std::vector<Gene> gs;
                            gs.push_back(g);
                            keepGenes[g.getName()] = gs;
                        }else if( keepTandemDuplication && geneScores[g.getName()] == thisScore ){
                            geneScores[g.getName()] = thisScore;
                            keepGenes[g.getName()].push_back(g);
                        }
                    }else{
                        geneScores[g.getName()] = thisScore;
                        std::vector<Gene> gs;
                        gs.push_back(g);
                        keepGenes[g.getName()] = gs;
                    }
                //}
            }
        }
    }
    std::set<std::string> geneNames;
    std::set<std::string> transcriptNames;
    std::vector<Gene> outputGenes;
    for( std::map<std::string, std::vector<Gene>>::iterator it = keepGenes.begin(); it!=keepGenes.end(); ++ it ){
        for( Gene g : it->second ){
            outputGenes.push_back(g);
        }
    }
    std::sort(outputGenes.begin(), outputGenes.end(), [](Gene a, Gene b) {
        if( a.getChromeSomeName() != b.getChromeSomeName() ){ return a.getChromeSomeName() < b.getChromeSomeName(); } else{return a.getStart() < b.getStart();}
    });
    for( Gene g : outputGenes ){
        outputGffRecords(g, ofile, geneNames, transcriptNames  );
    }
    ofile.close();
}




/**
 * my version of Quota-alignment using the longest alignment method
 * This longest alignment is designed as a local alignment model
 * And a lot of ideas were borrowed from the DAGchainder method
 * */
void myOrthologPairsSortQuota( std::vector<OrthologPair2> & pairedSimilarFragments){
    std::sort(pairedSimilarFragments.begin(), pairedSimilarFragments.end(), [](OrthologPair2 a, OrthologPair2 b) {
        return a < b;
    });
}
struct Path{
    double score;
    int index;
};
// since we will change pairedSimilarFragments, so do not use reference C++ data type here
void longestPathQuota (std::vector<OrthologPair2> pairedSimilarFragments, std::vector<std::vector<OrthologPair2>> & sortedOrthologPairChains,
        const bool & keepTandemDuplication, double & INDEL_SCORE, double & GAP_OPEN_PENALTY,
        double & MIN_ALIGNMENT_SCORE, const int & MAX_DIST_BETWEEN_MATCHES, int & refMaximumTimes, int & queryMaximumTimes, std::map<std::string, int> & refTimes,
        std::map<std::string, int> & queryTimes ){
    std::vector <Path> high;
    std::vector <int> ans;
    bool done;
    int n, i, j;
    do{
        done = true;
        n=pairedSimilarFragments.size();
        double scoreArray [n]; // arrays of scores
        int prev [n];  // index of previous node in longest path

        for (int idx = 0; idx < n; ++idx) {
            scoreArray[idx] = pairedSimilarFragments[idx].getScore();
            prev[idx] = -1;
        }
        //std::cout << "line 548" << std::endl;
        for (int idx = 1; idx < n; ++idx) {
            double thisIndexScore = scoreArray[idx];
            for (int jdx = idx - 1; jdx >= 0; --jdx) {// checking all previous nodes
                // Because we swapped asm/query start position so that inversions were all increasing,
                // we should always be on the diagonal.  If not, then we filter it.
                // This gets rid of the noise, while preserving the inversions on
                // the diagonal
                // Are only looking at positions previous to our current "idx" position
                if(  pairedSimilarFragments[idx].getQueryChr()==pairedSimilarFragments[jdx].getQueryChr()
                    && pairedSimilarFragments[idx].getRefChr()==pairedSimilarFragments[jdx].getRefChr() ){
                    if( pairedSimilarFragments[idx].getStrand() == pairedSimilarFragments[jdx].getStrand() ){
                        // the node one the chain should be in the same STRAND, if not this is an INDEL
                        if ( pairedSimilarFragments[idx].getStrand() == POSITIVE  ) { //same strand
                            int ref_del = pairedSimilarFragments[idx].getRefId() - pairedSimilarFragments[jdx].getRefId()-1;
                            //ref_del = ref_del * 2;
                            int query_del = pairedSimilarFragments[idx].getQueryId() - pairedSimilarFragments[jdx].getQueryId()-1;
                            int distance = (int)(((ref_del+query_del)+abs(ref_del-query_del))/(2.5));

                            if ( abs(ref_del) > MAX_DIST_BETWEEN_MATCHES && abs(query_del) > MAX_DIST_BETWEEN_MATCHES) {
                                // if this position is too large then the last node of j could not be i and the chain restart from j
                                break;
                            }
                            if( abs(abs(ref_del) - abs(query_del)) > MAX_DIST_BETWEEN_MATCHES ){
                                break;
                            }
                            if ( abs(ref_del) > MAX_DIST_BETWEEN_MATCHES || abs(query_del) > MAX_DIST_BETWEEN_MATCHES) { // is this for reversion?? I guess
                                continue;
                            }

                            double thisScore = thisIndexScore;
                            if( pairedSimilarFragments[idx].getRefId() != pairedSimilarFragments[jdx].getRefId() ){
                                thisScore += scoreArray[jdx];
                            }
                            if( distance > 0 ){
                                thisScore = thisScore + GAP_OPEN_PENALTY + INDEL_SCORE*distance;
                            }
                            if ( thisScore > scoreArray[idx] &&
                                 pairedSimilarFragments[jdx].getQueryId() < pairedSimilarFragments[idx].getQueryId() ){
                                scoreArray[idx] = thisScore;
                                prev[idx] = jdx;
                            } else if ( pairedSimilarFragments[jdx].getRefId() == pairedSimilarFragments[idx].getRefId()
                                        && keepTandemDuplication
                                        && abs(pairedSimilarFragments[jdx].getQueryId() - pairedSimilarFragments[idx].getQueryId())==1 ) {
                                scoreArray[idx] = scoreArray[jdx];
                                prev[idx] = jdx;
                            }
                        } else if ( pairedSimilarFragments[idx].getStrand() == NEGATIVE ) { // reversion
                            int ref_del = pairedSimilarFragments[idx].getRefId() - pairedSimilarFragments[jdx].getRefId()-1;
                            //ref_del = ref_del * 2;
                            int query_del = pairedSimilarFragments[idx].getQueryId() - pairedSimilarFragments[jdx].getQueryId()-1;
                            int distance = (int)(((ref_del*2+query_del)+abs(ref_del*2-query_del))/(2.5));
                            if ( abs(ref_del) > MAX_DIST_BETWEEN_MATCHES && abs(query_del) > MAX_DIST_BETWEEN_MATCHES) {
                                // if this position is too large then the last node of j could not be i and the chain restart from j
                                break;
                            }
                            if( abs(abs(ref_del) - abs(query_del)) > MAX_DIST_BETWEEN_MATCHES ){
                                break;
                            }
                            if ( abs(ref_del) > MAX_DIST_BETWEEN_MATCHES || abs(query_del) > MAX_DIST_BETWEEN_MATCHES) { // is this for reversion?? I guess
                                continue;
                            }
                            double thisScore = thisIndexScore;
                            if( pairedSimilarFragments[idx].getRefId() != pairedSimilarFragments[jdx].getRefId() ){
                                thisScore += scoreArray[jdx];
                            }
                            if( distance > 0 ){
                                thisScore = thisScore + GAP_OPEN_PENALTY + INDEL_SCORE*distance;
                            }
                            if ( thisScore > scoreArray[idx] &&
                                 pairedSimilarFragments[jdx].getQueryId() > pairedSimilarFragments[idx].getQueryId()){
                                scoreArray[idx] = thisScore;
                                prev[idx] = jdx;
                            }else if ( pairedSimilarFragments[jdx].getRefId() == pairedSimilarFragments[idx].getRefId()
                                       && keepTandemDuplication
                                       && abs(pairedSimilarFragments[jdx].getQueryId() - pairedSimilarFragments[idx].getQueryId())==1 ) {
                                scoreArray[idx] = scoreArray[jdx];
                                prev[idx] = jdx;
                            }
                        }
                    }
                }
            }
        }
        //std::cout << "line 608" << std::endl;
        high.clear();
        for (int idx = 0; idx < n; ++idx) {
            if ( scoreArray[idx] > MIN_ALIGNMENT_SCORE){
                Path p;
                p.index=idx;
                p.score=scoreArray[idx];
                high.push_back(p);
            }
        }
        //std::cout << "line 618" << std::endl;
        std::sort(high.begin(), high.end(), [](Path a, Path b) {
            return a.score > b.score;
        });
        //std::cout << "line 622 high.size(): " << high.size() << " highest value " << high[0].score << std::endl;
        if( high.size() > 0 ){
            done = false;
            i=0;
            if  (prev[high[i].index] != -2) {
                // trace back begin
                double subMaxScore = 0.0;
                int startIndex = high[i].index;
                double currentScore=0.0;
                for  (j = high[i].index;  prev[j] >= 0;  j = prev[j]) {
                    double scoreDifference = scoreArray[j]-scoreArray[prev[j]];
                    currentScore+=scoreDifference;
                    if( subMaxScore < currentScore ){
                        subMaxScore=currentScore;
                        startIndex=j;
                    }
                }
                double scoreDifference = scoreArray[j];  // this one should be the score of the node itself
                currentScore+=scoreDifference;
                if( subMaxScore < currentScore ){
                    subMaxScore=currentScore;
                    startIndex=j;
                }

                assert(abs(subMaxScore-high[i].score)<0.01);
                // trace back end
                ans.clear();
                for  (j = high[i].index;  j>=startIndex;  j = prev[j])  {
                    ans.push_back(j);
                }

                std::vector<OrthologPair2> chain;
                sortedOrthologPairChains.push_back(chain);
                reverse(ans.begin(), ans.end());
                int s=ans.size();
                std::string lastRef="";  // this is for keeping duplication purpose
                std::string lastQuery="";
                for( j=0; j<s; j++ ){
                    prev[ans[j]]=-2;
                    OrthologPair2 orthologPair2 = pairedSimilarFragments[ans[j]];
                    sortedOrthologPairChains[sortedOrthologPairChains.size()-1].push_back(orthologPair2);

                    //for quota begin
                    if( pairedSimilarFragments[ans[j]].getReferenceGeneName() != lastRef ){  // do not count duplications
                        if( refTimes.find(pairedSimilarFragments[ans[j]].getReferenceGeneName()) != refTimes.end() ){
                            refTimes[pairedSimilarFragments[ans[j]].getReferenceGeneName()]=refTimes[pairedSimilarFragments[ans[j]].getReferenceGeneName()]+1;
                        }else{
                            refTimes[pairedSimilarFragments[ans[j]].getReferenceGeneName()]=1;
                        }
                    }
                    if( pairedSimilarFragments[ans[j]].getQueryGeneName() != lastQuery ){
                        if( queryTimes.find( pairedSimilarFragments[ans[j]].getQueryGeneName()) != queryTimes.end() ){
                            queryTimes[ pairedSimilarFragments[ans[j]].getQueryGeneName()] = queryTimes[ pairedSimilarFragments[ans[j]].getQueryGeneName()] + 1;
                        } else {
                            queryTimes[ pairedSimilarFragments[ans[j]].getQueryGeneName()] = 1;
                        }
                    }
                    lastRef = pairedSimilarFragments[ans[j]].getReferenceGeneName();
                    lastQuery = pairedSimilarFragments[ans[j]].getQueryGeneName();
                    //for quota end
                }
            }
        }
        if( !done ){
            for( i=j=0; i<n; i++ ){
                if( prev[i] != -2 &&
                    ( refTimes.find(pairedSimilarFragments[i].getReferenceGeneName()) == refTimes.end() || refTimes[pairedSimilarFragments[i].getReferenceGeneName()]<refMaximumTimes ) &&
                        (queryTimes.find(pairedSimilarFragments[i].getQueryGeneName()) ==queryTimes.end() || queryTimes[pairedSimilarFragments[i].getQueryGeneName()]<queryMaximumTimes )){
                    if( i != j ){ // those elements should be maintained for next loop
                        pairedSimilarFragments[j]=pairedSimilarFragments[i];
                    }
                    ++j;
                }
            }
            pairedSimilarFragments.resize(j);
        }
    }while(! done);
}


void generateLongestQuotaOutput( const std::string & referenceGffFile, const std::string & queryNewGffFile,
                            const std::string & queryGenomeFile, const std::string & outputGffFile, const int & minIntron,
                            const bool & keepTandemDuplication,
                            std::map<std::string, std::string>& parameters, const double & syntenicScore,
                            const double & orfScore, const double & dropLengthThredshold,
                            double & INDEL_SCORE, double & GAP_OPEN_PENALTY, double & MIN_ALIGNMENT_SCORE, int & MAX_DIST_BETWEEN_MATCHES,
                            int & refMaximumTimes, int & queryMaximumTimes, const bool & onlySyntenic, const bool & sortOutPutGffBycoordinate ){
    //std::cout << " line 676" << std::endl;
    std::map<std::string, Fasta> queryGenome;
    readFastaFile(queryGenomeFile, queryGenome);
    //std::cout << " line 679" << std::endl;
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    std::map<std::string, std::vector<Gene> > queryGenes;
    readGffFileWithEveryThing (queryNewGffFile, queryGenes);
    //std::cout << " line 683" << std::endl;
    std::map<std::string, std::vector<std::string> > referenceGeneNameMap;
    std::map<std::string, Gene > referenceGeneHashMap;
    std::map<std::string, Transcript> referenceTranscriptHashMap;
    readGffFileWithEveryThing ( referenceGffFile, referenceGeneNameMap, referenceGeneHashMap, referenceTranscriptHashMap);
    //std::cout << " line 688" << std::endl;
    std::map<std::string, std::vector<Gene> > referenceGenes;
    readGffFileWithEveryThing (referenceGffFile, referenceGenes);
    //std::cout << " line 691" << std::endl;


    //filter identical records from query begin
    for( std::map<std::string, std::vector<Gene> >::iterator it=queryGenes.begin(); it!=queryGenes.end(); ++it){
        int i, j, n;
        n=it->second.size();
        for ( i = j = 0;  i < n;  ++i ) {
            if ( it->second[i].getName()==it->second[j].getName() && it->second[i].getStart()==it->second[j].getStart()
                 && it->second[i].getEnd()==it->second[j].getEnd() ) {

            }else{
                ++j;
                if  (i != j){
                    it->second[j] = it->second[i];
                }
            }
        }
        ++j;
        it->second.resize(j);
        queryGenes[it->first]=it->second;
    }
    //filter identical records from query end


    std::map<std::string, std::map<std::string, int>> referenceGenesIndex;
    for( std::map<std::string, std::vector<Gene> >::iterator it=referenceGenes.begin(); it!=referenceGenes.end(); ++it ){
        std::map<std::string, int> gv;
        referenceGenesIndex[it->first]=gv;
        for( int id=0; id < it->second.size(); ++id ){
            referenceGenesIndex[it->first][it->second[id].getName()]=id;
            //std::cout << "chr: " << it->first << " position " << it->second[id].getStart() << " id: " << id << "  " << std::endl;
        }
    }

    std::ofstream ofile;
    ofile.open(outputGffFile);
    std::map<std::string, std::vector<Gene>> keepGenes;
    std::map<std::string, double> geneScores;
    //std::cout << " line 697" << std::endl;
    std::vector<OrthologPair2> orthologPairs;
    for(std::map<std::string, std::vector<Gene>>::iterator it1=queryGenes.begin(); it1!=queryGenes.end(); ++it1){
        if( queryGenome.find(it1->first)!=queryGenome.end() ){
            for( int i=0; i< it1->second.size(); ++i ){
                if( referenceGeneHashMap.find(it1->second[i].getName())!=referenceGeneHashMap.end() ){
                    int refStartPos=referenceGeneHashMap[it1->second[i].getName()].getStart();
                    int refEndPos   =referenceGeneHashMap[it1->second[i].getName()].getEnd();
                    int queryStartPos=it1->second[i].getStart();
                    int queryEndPos  =it1->second[i].getEnd();
                    double thisScore=0;
                    for(  int index=0; index<it1->second[i].getTranscripts().size(); ++index ){
                        Transcript transcript=it1->second[i].getTranscripts()[index];
                        TranscriptUpdateCdsInformation(transcript, queryGenome);
                        checkOrfState(transcript, queryGenome, nucleotideCodeSubstitutionMatrix, minIntron);
                        if( transcript.getIfOrfShift() ){
                            // this score is for gene, if there is one ORF conserved transcript, then give this gene a positive score
                        }else{
                            thisScore=orfScore;
                        }
                        it1->second[i].getTranscripts()[index]=transcript;
                        //break;
                    }
                    double lengthRatio;
                    if( it1->second[i].getEnd() == it1->second[i].getStart() ){
                        lengthRatio=0.0;
                    }else{
                        lengthRatio=(double)((double)it1->second[i].getEnd()-(double)it1->second[i].getStart())/((double)referenceGeneHashMap[it1->second[i].getName()].getEnd()-(double)referenceGeneHashMap[it1->second[i].getName()].getStart());
                    }
                    if( lengthRatio > 1 ){
                        lengthRatio=1/lengthRatio;
                    }
                    thisScore += lengthRatio;
                    STRAND strand;
                    if( it1->second[i].getStrand() == referenceGeneHashMap[it1->second[i].getName()].getStrand() ){
                        strand=POSITIVE;
                    }else{
                        strand=NEGATIVE;
                    }
                    if( lengthRatio > dropLengthThredshold && referenceGenesIndex.find(referenceGeneHashMap[it1->second[i].getName()].getChromeSomeName())!=referenceGenesIndex.end() &&
                        referenceGenesIndex[referenceGeneHashMap[it1->second[i].getName()].getChromeSomeName()].find(it1->second[i].getName())!=referenceGenesIndex[referenceGeneHashMap[it1->second[i].getName()].getChromeSomeName()].end()){
                        int indexj = referenceGenesIndex[referenceGeneHashMap[it1->second[i].getName()].getChromeSomeName()][it1->second[i].getName()]; //because referenceGenes has been sorted, so it's index makes sense
                        OrthologPair2 orthologPair( referenceGeneHashMap[it1->second[i].getName()].getChromeSomeName(), it1->first,
                            refStartPos, refEndPos, queryStartPos, queryEndPos, thisScore, strand, indexj, i, it1->second[i].getName(), it1->second[i].getName() );
                        orthologPairs.push_back(orthologPair);
                    }
                }
            }
        }
    }
    std::map<std::string, int> refTimes;
    std::map<std::string, int> queryTimes;
    //std::cout << "line 749" << std::endl;
    myOrthologPairsSortQuota( orthologPairs);
    //std::cout << "line 751 orthologPairs.size(): " << orthologPairs.size() << std::endl;
    std::vector<std::vector<OrthologPair2>> sortedOrthologPairChains;
    longestPathQuota (orthologPairs, sortedOrthologPairChains, keepTandemDuplication, INDEL_SCORE, GAP_OPEN_PENALTY,
            MIN_ALIGNMENT_SCORE, MAX_DIST_BETWEEN_MATCHES, refMaximumTimes, queryMaximumTimes, refTimes, queryTimes);
    //std::cout << "line 753 sortedOrthologPairChains.size(): " << sortedOrthologPairChains.size() << std::endl;

    std::set<std::string> geneNames;
    std::set<std::string> transcriptNames;
    // output those syntenic region records
    int chainIndex = 0;
    for( std::vector<OrthologPair2> sortedOrthologPairs : sortedOrthologPairChains ){
        if( onlySyntenic && !sortOutPutGffBycoordinate){
            ofile << "#start chain " << chainIndex << std::endl;
        }
        for( OrthologPair2 orthologPair : sortedOrthologPairs ) {
            Gene g = queryGenes[orthologPair.getQueryChr()][orthologPair.getQueryId()];
            double thisScore = orthologPair.getScore() + syntenicScore;
            if (geneScores.find(g.getName()) != geneScores.end()) {
                if( geneScores[g.getName()] < thisScore ){
                    geneScores[g.getName()] = thisScore;
                }
                std::vector<Gene> gs;
                gs.push_back(g);
                keepGenes[g.getName()] = gs;
            } else {
                geneScores[g.getName()] = thisScore;
                std::vector<Gene> gs;
                gs.push_back(g);
                keepGenes[g.getName()] = gs;
            }
            if( onlySyntenic && !sortOutPutGffBycoordinate ){
                outputGffRecords(g, ofile, geneNames, transcriptNames  );
            }
        }
        if( onlySyntenic && !sortOutPutGffBycoordinate){
            ofile << "#end chain " << chainIndex << std::endl;
        }
        ++chainIndex;
    }
    if( onlySyntenic && !sortOutPutGffBycoordinate ){

    }else{
        if ( onlySyntenic ){

        }else {
            //for those genes not in syntenic region, only output once
            for (std::map<std::string, std::vector<Gene>>::iterator it1 = queryGenes.begin();
                 it1 != queryGenes.end(); ++it1) {
                for (Gene g : it1->second) {
                    double thisScore = 0;
                    for (int index = 0; index < g.getTranscripts().size(); ++index) {
                        Transcript transcript = g.getTranscripts()[index];
                        //TranscriptUpdateCdsInformation(transcript, queryGenome);
                        //checkOrfState(transcript, queryGenome, nucleotideCodeSubstitutionMatrix, minIntron);
                        if (transcript.getIfOrfShift()) {
                            // this score is for gene, if there is one ORF conserved transcript, then give this gene a positive score
                        } else {
                            thisScore = orfScore;
                        }
                        //g.getTranscripts()[index]=transcript;
                    }

                    double lengthRatio;
                    if (g.getEnd() == g.getStart()) {
                        lengthRatio = 0.0;
                    } else {
                        lengthRatio = (double) ((double) g.getEnd() - (double) g.getStart()) /
                                      ((double) referenceGeneHashMap[g.getName()].getEnd() -
                                       (double) referenceGeneHashMap[g.getName()].getStart());
                    }
                    if (lengthRatio > 1) {
                        lengthRatio = 1 / lengthRatio;
                    }
                    thisScore += lengthRatio;
                    if (geneScores.find(g.getName()) != geneScores.end()) {
                        if (geneScores[g.getName()] < thisScore) {
                            geneScores[g.getName()] = thisScore;
                            std::vector<Gene> gs;
                            gs.push_back(g);
                            keepGenes[g.getName()] = gs;
                        }/*else if( keepTandemDuplication && geneScores[g.getName()] == thisScore ){
                            geneScores[g.getName()] = thisScore;
                            keepGenes[g.getName()].push_back(g);
                        }*/
                    } else {
                        geneScores[g.getName()] = thisScore;
                        std::vector<Gene> gs;
                        gs.push_back(g);
                        keepGenes[g.getName()] = gs;
                    }
                }
            }
        }
        std::vector<Gene> outputGenes;
        for( std::map<std::string, std::vector<Gene>>::iterator it = keepGenes.begin(); it!=keepGenes.end(); ++ it ){
            for( Gene g : it->second ){
                outputGenes.push_back(g);
            }
        }
        std::sort(outputGenes.begin(), outputGenes.end(), [](Gene a, Gene b) {
            if( a.getChromeSomeName() != b.getChromeSomeName() ){ return a.getChromeSomeName() < b.getChromeSomeName(); } else{return a.getStart() < b.getStart();}
        });
        for( Gene g : outputGenes ){
            outputGffRecords(g, ofile, geneNames, transcriptNames  );
        }
    }
    ofile.close();
}







/**
 * the DAGchainer algorithm
 * this is another syntenic method
 * reference https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bth397
 * the source code was adopted from the DAGchainer software
 * */

struct Score_t {
//    std::string geneIdx, geneIdy;
    std::string chrx, chry;
    int idx, idy; // database and query gene index in the vector
    double score;
    bool sameStrand;
    bool operator < (const Score_t & node)  const {
        return  ( (idx < node.idx) || (idx == node.idx && idy < node.idy) );
    }
};

struct Path_t {
    double score;
    int rc;  // sum of row and column of last entry
    int sub;
};

bool Descending_Score (const Path_t & a, const Path_t & b) {
    return ( (a.score > b.score) || (a.score == b.score && a.rc > b.rc) );
}

// here score could not be reference, since this function will change it's value
static void chainer (std::vector<Score_t> score, std::vector<std::vector<Score_t>> & chains, const int & MAX_DIST_BETWEEN_MATCHES /*max gap in the term of number of genes*/,
                           const int & BP_GAP_SIZE, const double & INDEL_SCORE, const double & GAP_OPEN_PENALTY, const double & MIN_ALIGNMENT_SCORE, const bool & keepTandemDuplication ) {
    //  Find and output highest scoring chains in score treating
    //  it as a DAG
    std::vector <double> path_score;
    std::vector <int> from, ans;
    std::vector <Path_t> high;
    bool  done;
    int  i, j, m, n, s;
    do {
        done = true;
        n = score.size();
        //std::cout << "line 975 n: " << n << std::endl;
        path_score.resize(n);  // the path score of this node
        from.resize (n);
        for (i = 0;  i < n;  ++i)  {
            path_score[i] = score[i].score;
            from[i] = -1;
//            std::cout << " line 981 " << score[i].score << std::endl;
        }
        //std::cout << "line 982 " << std::endl;
        for (j = 1; j < n; ++j) {
            //std::cout << "line 984 j " << j << std::endl;
            double thisIndexScore = path_score[j];
            for (i=j-1; i >= 0; --i) {
//                for( int z=i; z<j; ++z ){  //todo here I want to check the distance between the current vector index and the gff index distance, to give penalty to those records that have already included in the output list
//
//                }
                //std::cout << "line 990 " << std::endl;
                if ( score[i].chrx==score[j].chrx && score[i].chry==score[j].chry && score[i].sameStrand==score[j].sameStrand && (  // same chromosome. same strand and the position makes sense
                        (score[i].sameStrand==true && score[i].idy<=score[j].idy) || (score[i].sameStrand==false && score[i].idy>=score[j].idy)
                        )){
                    int del_x = score[j].idx - score[i].idx - 1;  // this is the position information
                    int del_y = score[j].idy - score [i].idy - 1;  // this is the position information
                    if ( abs(del_x) > MAX_DIST_BETWEEN_MATCHES && abs(del_y) > MAX_DIST_BETWEEN_MATCHES) {
                        // if this position is too large then the last node of j could not be i and the chain restart from j
                        break;
                    }
                    if( abs(abs(del_x) - abs(del_y)) > MAX_DIST_BETWEEN_MATCHES ){
                        break;
                    }
                    if ( abs(del_x) > MAX_DIST_BETWEEN_MATCHES || abs(del_y) > MAX_DIST_BETWEEN_MATCHES) { // is this for reversion?? I guess
                        continue;
                    }
                    int num_gaps = (int) ( ((del_x*2 + del_y)+ abs(del_x*2-del_y)) / (2 * BP_GAP_SIZE)  + 0.5); // this function is friendly for reversion
                    double x = thisIndexScore;
                    if (  score[i].idx !=  score[j].idx ) { // do not give match score for tandem duplication
                        x += path_score[i];
                    }
                    if (num_gaps > 0) {
                        // Affine gap penalty:
                        // penalty = open + (num_gaps * extension_penalty)
                        x = x + GAP_OPEN_PENALTY + INDEL_SCORE*num_gaps;
                    }
                    //std::cout << "line 1016 " << std::endl;
                    if (x > path_score[j] && score[i].idx!=score[j].idx ) {
                        path_score [j] = x;
                        from [j] = i;
                    } else if ( score[i].idx==score[j].idx && keepTandemDuplication && abs(score[i].idy-score[j].idy)==1 ) {
                        path_score [j] = path_score [i];
                        from [j] = i;
                    }
                }
                //std::cout << "line 1025 " << std::endl;
            }
        }
        //std::cout << "line 1026" << std::endl;
        high.clear();
        //std::cout << "line 1030" << std::endl;
        for  (i = 0;  i < n;  ++i) {
//            std::cout << "i " << i << " line 1032 MIN_ALIGNMENT_SCORE " << MIN_ALIGNMENT_SCORE << " path_score[i] " << path_score[i] << " n " << n<< std::endl;
            if  ( path_score[i] >= MIN_ALIGNMENT_SCORE )  {
  //              std::cout << "line 1033" << std::endl;
                Path_t  p;
                p.score = path_score[i];
                p.sub = i;
                p.rc = score[i].idx + score[i].idy;
                high.push_back(p);
            }
        }
        //std::cout << "line 1035 high.size(): " << high.size() << " highest value " << high[0].score << std::endl;
        sort (high.begin(), high.end(), Descending_Score); // here we sort the score, so we output the highest score path firstly,
        // and once the path has been out putted, all the elements one the score will be removed from the vector

        m = high.size();
        if( m > 0){
            for  (i = 0;  i < m;  ++i) {
                if  (from[high[i].sub] != -2) {
                    // trace back begin
                    double subMaxScore = 0.0;
                    int startIndex = high[i].sub;
                    double currentScore=0.0;
                    for  (j = high[i].sub;  from[j] >= 0;  j = from[j]) {
                        double scoreDifference = path_score[j]-path_score[from[j]];
                        currentScore+=scoreDifference;
                        if( subMaxScore < currentScore ){
                            subMaxScore=currentScore;
                            startIndex=j;
                        }
                    }
                    double scoreDifference = score[j].score;  // this one should be the score of the node itself
                    currentScore+=scoreDifference;
                    if( subMaxScore < currentScore ){
                        subMaxScore=currentScore;
                        startIndex=j;
                    }

                    if  (from[j] == -2)  { // if there is any element in two chain stop the current analysis
                        done = false;
                        break;
                    } else {
                        //if( subMaxScore!=high[i].score ){
                        //std::cout << "subMaxScore: " << subMaxScore << " high[i].score " << high[i].score << std::endl;
                        //}
                        assert(abs(subMaxScore-high[i].score)<0.01);
                        // trace back end
                        ans.clear();
                        for  (j = high[i].sub;  j>=startIndex;  j = from[j])  {
                            ans.push_back(j);
                        }

                        std::vector<Score_t> chain;
                        chains.push_back(chain);
                        reverse(ans.begin(), ans.end ());
                        s = ans.size();
                        for  (j = 0;  j < s;  ++j) {
                            from[ans[j]] = -2; // we should remove those elements
                            chains[chains.size()-1].push_back(score[ans[j]]);
                        }
                    }
                }
            }
            if  (! done)  { // if it is not done,
                for ( i = j = 0;  i < n;  ++i ) {
                    if  (from [i] != -2) {
                        if  (i != j){
                            score[j] = score[i];
                        }
                        j++;
                    }
                }
                score.resize(j); // score is the real input, here we delete those things already on the chain
            }
        }
    } while  (! done);
}


void generateDagChainerOutput( const std::string & referenceGffFile, const std::string & queryNewGffFile,
                            const std::string & queryGenomeFile, const std::string & outputGffFile, const int & minIntron,
                            const bool & keepTandemDuplication,
                            std::map<std::string, std::string>& parameters, const double & syntenicScore, const double & orfScore,
                            const double & dropLengthThredshold, int MAX_DIST_BETWEEN_MATCHES /*max gap in the term of number of genes*/,
                            int BP_GAP_SIZE, double INDEL_SCORE, double GAP_OPEN_PENALTY, double MIN_ALIGNMENT_SCORE, const bool & onlySyntenic, const bool & sortOutPutGffBycoordinate){
//    std::cout << "longest line 308" << std::endl;
    std::map<std::string, Fasta> queryGenome;
    readFastaFile(queryGenomeFile, queryGenome);
//    std::cout << "longest line 311" << std::endl;
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
//    std::cout << "longest line 313" << std::endl;
    std::map<std::string, std::vector<Gene> > queryGenes;
//    std::cout << "longest line 315" << std::endl;
    readGffFileWithEveryThing (queryNewGffFile, queryGenes);
//    std::cout << "longest line 317" << std::endl;
    std::map<std::string, std::vector<std::string> > referenceGeneNameMap;
    std::map<std::string, Gene > referenceGeneHashMap;
    std::map<std::string, Transcript> referenceTranscriptHashMap;
//    std::cout << "longest line 321" << std::endl;
    readGffFileWithEveryThing ( referenceGffFile, referenceGeneNameMap, referenceGeneHashMap, referenceTranscriptHashMap);
//    std::cout << "longest line 323" << std::endl;

    std::map<std::string, std::vector<Gene> > referenceGenes;
    readGffFileWithEveryThing (referenceGffFile, referenceGenes);

    std::ofstream ofile;
    ofile.open(outputGffFile);
//    std::cout << "longest line 328" << std::endl;

    // DAGchainer begin
//    std::map<int, Score_t> map_score;
    std::vector<Score_t>  score; //list holds all inputted matches.
    int  i, j, n;
//    int score_id=0;
    // DAGchainer end


    // construct this data structure, it is needed to get the order of a certain gene on the corresponding chromosome
    std::map<std::string, std::map<std::string, int>> referenceGenesIndex;
    for( std::map<std::string, std::vector<Gene> >::iterator it=referenceGenes.begin(); it!=referenceGenes.end(); ++it ){
        std::map<std::string, int> gv;
        referenceGenesIndex[it->first]=gv;
        for( int id=0; id < it->second.size(); ++id ){
            referenceGenesIndex[it->first][it->second[id].getName()]=id;
            //std::cout << "chr: " << it->first << " geneName " << it->second[id].getName() << " id " << id << std::endl;
        }
    }

    //filter identical records from query begin
    for( std::map<std::string, std::vector<Gene> >::iterator it=queryGenes.begin(); it!=queryGenes.end(); ++it){
        n=it->second.size();
        for ( i = j = 0;  i < n;  ++i ) {
            if ( it->second[i].getName()==it->second[j].getName() && it->second[i].getStart()==it->second[j].getStart()
                 && it->second[i].getEnd()==it->second[j].getEnd() ) {

            }else{
                ++j;
                if  (i != j){
                    it->second[j] = it->second[i];
                }
            }
        }
        ++j;
        it->second.resize(j);
        queryGenes[it->first]=it->second;
    }
    //filter identical records from query end


    for(std::map<std::string,std::vector<Gene>>::iterator it1=queryGenes.begin(); it1!=queryGenes.end(); it1++){
        if( queryGenome.find(it1->first)!=queryGenome.end() ){
            for( i=0; i< it1->second.size(); ++i ){
                if( referenceGeneHashMap.find(it1->second[i].getName())!=referenceGeneHashMap.end() ){
                    std::cout << "line 1207 " << it1->second[i].getName() << std::endl;
                    double thisScore=0;
                    for(  int index=0; index<it1->second[i].getTranscripts().size(); ++index ){
                        Transcript transcript=it1->second[i].getTranscripts()[index];
                        TranscriptUpdateCdsInformation(transcript, queryGenome);
                        checkOrfState(transcript, queryGenome, nucleotideCodeSubstitutionMatrix, minIntron);
                        if( transcript.getIfOrfShift() ){
                            // this score is for gene, if there is one ORF conserved transcript, then give this gene a positive score
                        }else{
                            thisScore=orfScore;
                        }
                        it1->second[i].getTranscripts()[index]=transcript; //break;
                    }
                    double lengthRatio;
                    if( it1->second[i].getEnd() == it1->second[i].getStart() ){
                        lengthRatio=0.0;
                    }else{
                        lengthRatio=(double)((double)it1->second[i].getEnd()-(double)it1->second[i].getStart())/((double)referenceGeneHashMap[it1->second[i].getName()].getEnd()-(double)referenceGeneHashMap[it1->second[i].getName()].getStart());
                    }
                    if( lengthRatio > 1 ){
                        lengthRatio=1/lengthRatio;
                    }
                    thisScore += lengthRatio;
                    if( lengthRatio > dropLengthThredshold && referenceGenesIndex.find(referenceGeneHashMap[it1->second[i].getName()].getChromeSomeName())!=referenceGenesIndex.end() &&
                            referenceGenesIndex[referenceGeneHashMap[it1->second[i].getName()].getChromeSomeName()].find(it1->second[i].getName())!=referenceGenesIndex[referenceGeneHashMap[it1->second[i].getName()].getChromeSomeName()].end() ){
                        // OrthologPair orthologPair( targetIndex, refStartPos, refEndPos, queryStartPos, queryEndPos, thisScore, strand );
                        Score_t  s;
//                        s.geneIdx=it1->second[i].getName();
//                        s.geneIdy=it1->second[i].getName();
                        s.chrx=referenceGeneHashMap[it1->second[i].getName()].getChromeSomeName(); // reference
                        s.chry=it1->first;
                        s.idx = referenceGenesIndex[referenceGeneHashMap[it1->second[i].getName()].getChromeSomeName()][it1->second[i].getName()];
                        s.idy = i;
                        s.score = thisScore;
                        if( it1->second[i].getStrand() == referenceGeneHashMap[it1->second[i].getName()].getStrand() ){
                            s.sameStrand=true;
                        }else{
                            s.sameStrand=false;
                        }
                        score.push_back(s); //copy match to score list. // this is from the input file
                    }
                }
            }
        }
    }
    n = score.size();

    sort (score.begin(), score.end()); // increase order with x

    // for duplication removing purpose
    for  (i = 1, j = 0;  i < n;  i++) {
        if  (score[j].idx == score[i].idx && score[j].idy == score[i].idy  && score[j].chrx == score[i].chrx && score[j].chry == score[i].chry)  {
            // consecutive score entries have the same identities.
            fprintf (stderr,
                     "Duplicate score entries:  (%d,%d,%.1f) and (%d,%d,%.1f)\n",
                     score[j].idx, score[j].idy, score[j].score,
                     score[i].idx, score[i].idy, score[i].score);
            fprintf (stderr, "Discarding latter\n");
        } else {
            // not the same
            j++;
            if  (i != j) {
                score[j] = score[i];
            }
        }
    }
    j++;
    score.resize (j);
    std::vector<std::vector<Score_t>> chains;
    std::cout << "line 1242 score.size(): " << score.size() << std::endl;
    chainer (score, chains, MAX_DIST_BETWEEN_MATCHES, BP_GAP_SIZE, INDEL_SCORE, GAP_OPEN_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication);
    std::cout << "line 1244 chains.size : " << chains.size() << std::endl;
    std::map<std::string, std::vector<Gene>> keepGenes;
    std::map<std::string, double> geneScores;
    std::set<std::string> geneNames;
    std::set<std::string> transcriptNames;
    int chainIndex=0;
    for( std::vector<Score_t> chain : chains ){
        if( onlySyntenic && !sortOutPutGffBycoordinate ){
            ofile << "#start chain " << chainIndex << std::endl;
        }
        for( Score_t s_t : chain ){
            //outputGffRecords(it1->second[orthologPair.getQueryIndex()], ofile, geneNames, transcriptNames  );
            Gene g = (queryGenes[s_t.chry])[s_t.idy];
            //std::cout << "line 1123 g.name: " << g.getName() << std::endl;
            double thisScore = s_t.score + syntenicScore;
            if( geneScores.find(g.getName()) != geneScores.end() ){
                if( geneScores[g.getName()] < thisScore ) {
                    geneScores[g.getName()] = thisScore;
                }/*
                    std::vector<Gene> gs;
                    gs.push_back(g);
                    keepGenes[g.getName()] = gs;
                }else if( keepTandemDuplication && geneScores[g.getName()] == thisScore ){*/
                    //geneScores[g.getName()] = thisScore;
                    keepGenes[g.getName()].push_back(g);
                /*}*/
            }else{
                geneScores[g.getName()] = thisScore;
                std::vector<Gene> gs;
                gs.push_back(g);
                keepGenes[g.getName()] = gs;
            }
            if( onlySyntenic && !sortOutPutGffBycoordinate){
                outputGffRecords(g, ofile, geneNames, transcriptNames  );
            }
        }
        if( onlySyntenic && !sortOutPutGffBycoordinate){
            ofile << "#end chain " << chainIndex << std::endl;
        }
        ++chainIndex;
    }

    if( onlySyntenic && !sortOutPutGffBycoordinate){

    }else {
        if ( onlySyntenic ){

        }else {
            for (std::map<std::string, std::vector<Gene>>::iterator it1 = queryGenes.begin();
                 it1 != queryGenes.end(); ++it1) {
                for (i = 0; i < it1->second.size(); ++i) {
                    Gene g = it1->second[i];
                    double thisScore = 0;
                    for (int index = 0; index < g.getTranscripts().size(); ++index) {
                        Transcript transcript = g.getTranscripts()[index];
                        TranscriptUpdateCdsInformation(transcript, queryGenome);
                        checkOrfState(transcript, queryGenome, nucleotideCodeSubstitutionMatrix, minIntron);
                        if (transcript.getIfOrfShift()) {
                            // this score is for gene, if there is one ORF conserved transcript, then give this gene a positive score
                        } else {
                            thisScore = orfScore;
                        }
                        g.getTranscripts()[index] = transcript;
                    }
                    double lengthRatio;
                    if (g.getEnd() == g.getStart()) {
                        lengthRatio = 0.0;
                    } else {
                        lengthRatio = (double) ((double) g.getEnd() - (double) g.getStart()) /
                                      ((double) referenceGeneHashMap[g.getName()].getEnd() -
                                       (double) referenceGeneHashMap[g.getName()].getStart());
                    }
                    if (lengthRatio > 1) {
                        lengthRatio = 1 / lengthRatio;
                    }
                    thisScore += lengthRatio;
                    if (geneScores.find(g.getName()) != geneScores.end()) {
                        if (geneScores[g.getName()] < thisScore) {
                            geneScores[g.getName()] = thisScore;
                            std::vector<Gene> gs;
                            gs.push_back(g);
                            keepGenes[g.getName()] = gs;
                        } else if (keepTandemDuplication && geneScores[g.getName()] == thisScore) {
                            geneScores[g.getName()] = thisScore;
                            keepGenes[g.getName()].push_back(g);
                        }
                    } else {
                        geneScores[g.getName()] = thisScore;
                        std::vector<Gene> gs;
                        gs.push_back(g);
                        keepGenes[g.getName()] = gs;
                    }
                }
            }
        }

        std::vector<Gene> outputGenes;
        for( std::map<std::string, std::vector<Gene>>::iterator it = keepGenes.begin(); it!=keepGenes.end(); ++ it ){
            for( Gene g : it->second ){
                outputGenes.push_back(g);
            }
        }

        std::sort(outputGenes.begin(), outputGenes.end(), [](Gene a, Gene b) {
            if( a.getChromeSomeName() != b.getChromeSomeName() ){ return a.getChromeSomeName() < b.getChromeSomeName(); } else{return a.getStart() < b.getStart();}
        });

        for( Gene g : outputGenes ){
            outputGffRecords(g, ofile, geneNames, transcriptNames  );
        }
    }
    ofile.close();
}
