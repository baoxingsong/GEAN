//
// Created by Baoxing song on 20.10.18.
//

#include "organizeGffRecords.h"

bool exactlySame( Gene & g1, Gene & g2){
    if( g1.getStart() == g2.getStart() && g1.getEnd() == g2.getEnd() ){
        for( std::vector<Transcript>::iterator t1=g1.getTranscripts().begin();
             t1!=g1.getTranscripts().end(); ++t1 ) {
            bool included = false;
            for( std::vector<Transcript>::iterator t2=g2.getTranscripts().begin();
                 t2!=g2.getTranscripts().end(); ++t2 ) {
                if( t1->getStart()==t2->getStart() && t1->getEnd()==t2->getEnd() &&
                    t1->getCdsSequence().compare(t2->getCdsSequence())==0){
                    included = true;
                }
            }
            if( ! included  ){
                return false;
            }
        }
        return true;
    }
    return false;
}

void mergeToFirstOne( Gene & g1, Gene & g2){
    if( g1.getStrand() ==  g1.getStrand() && (
            ( g1.getStart() <= g2.getStart() && g2.getStart() <=g1.getEnd() ) ||
            ( g1.getStart() <= g2.getEnd() && g2.getEnd() <=g1.getEnd() ) ||
            ( g2.getStart() <= g1.getStart() && g1.getStart() <=g2.getEnd() ) ||
            ( g2.getStart() <= g1.getEnd() && g1.getEnd() <=g2.getEnd() )
    ) ){
        for( size_t i=0; i<g1.getTranscripts().size(); ++i){
            for( size_t j=0; j<g2.getTranscripts().size(); ++j) {
                if( g1.getTranscripts()[i].getName().compare(g2.getTranscripts()[j].getName())==0 && g1.getTranscripts()[i].getIfOrfShift() && (! g2.getTranscripts()[j].getIfOrfShift() ) ){
                    g1.getTranscripts()[i]=g2.getTranscripts()[j];
                }
            }
        }
        int start=g1.getStart();
        int end=g1.getEnd();
        for( size_t i=0; i<g1.getTranscripts().size(); ++i){
            if( start > g1.getTranscripts()[i].getStart() ){
                start = g1.getTranscripts()[i].getStart();
            }
            if( end < g1.getTranscripts()[i].getEnd() ){
                end =g1.getTranscripts()[i].getEnd();
            }
        }
        g1.setStart(start);
        g1.setEnd(end);
    }
}

void updateGeneInformation(std::map<std::string, std::vector<Gene> > & geneMap, NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix,
                           const size_t & minIntron,  std::map<std::string, Fasta> & querySequences){
    for( std::map<std::string, std::vector<Gene> >::iterator it0=geneMap.begin(); it0!=geneMap.end(); ++it0 ){
        for( std::vector<Gene>::iterator it1=it0->second.begin(); it1!=it0->second.end(); ++it1 ){
            if( it1->getTranscripts().size() > 0 ){
                for(  std::vector<Transcript>::iterator it=it1->getTranscripts().begin();
                      it!=it1->getTranscripts().end(); ++it){
                    int thisTranscriptStart=std::numeric_limits<int>::max();
                    int thisTranscriptEnd=0;
                    if( it->getExonVector().size() > 0 ) {
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
                int thisStart = it1->getTranscripts()[0].getStart();
                int thisEnd = it1->getTranscripts()[0].getEnd();
                for(  std::vector<Transcript>::iterator it=it1->getTranscripts().begin();
                      it!=it1->getTranscripts().end(); ++it){
                    if( thisStart> it->getStart() ){
                        thisStart = it->getStart();
                    }
                    if( thisEnd> it->getEnd() ){
                        thisEnd = it->getEnd();
                    }
                }
                it1->setStart(thisStart);
                it1->setEnd(thisEnd);
                if( it1->getTranscripts()[0].getSource().length()<1 ){
                    it1->getTranscripts()[0].setSource("LIFTOVER");
                }

                for( std::vector<Transcript>::iterator it=it1->getTranscripts().begin();
                     it!=it1->getTranscripts().end(); ++it ){

                    if( it->getSource().length()<1 ){
                        it->setSource("LIFTOVER");
                    }
                    if( it->getCdsVector().size() > 0 ){
                        TranscriptUpdateCdsInformation(*it, querySequences);
                        checkOrfState(*it, querySequences, nucleotideCodeSubstitutionMatrix, minIntron);
                    }
                }
            }
        }
    }
}
void removeDuplication(std::map<std::string, std::vector<Gene> > & geneMap, const size_t & minGene, std::map<std::string, std::set<int32_t>> & toRemove){
    for( std::map<std::string, std::vector<Gene> >::iterator it0=geneMap.begin(); it0!=geneMap.end(); ++it0 ){
        toRemove[it0->first]=std::set<int32_t>();
        int32_t j=0;
        if( it0->second.size()>0 ){
            if( it0->second[0].getTranscripts().size() == 0 ){
                toRemove[it0->first].insert(0);
            }else if ( ( it0->second[0].getEnd()-it0->second[0].getStart() ) < minGene ){
                toRemove[it0->first].insert(0);
                //std::cout << "remove zero length gene " << it0->second[0].getName() << std::endl;
            }
        }
        for( int32_t i=1; i<it0->second.size(); ++i ){
            j=i-1;
            if( it0->second[i].getTranscripts().size() == 0 ){
                toRemove[it0->first].insert(i);
            }else if ( ( it0->second[i].getEnd()-it0->second[i].getStart() ) < minGene ){
                toRemove[it0->first].insert(i);
                //std::cout << "remove zero length gene " << it0->second[i].getName() << std::endl;
            }else {
                while( toRemove[it0->first].find(j)!=toRemove[it0->first].end() && j>0 ){
                    --j;
                }
                if( j>=0 && it0->second[i].getName().compare(it0->second[j].getName()) == 0 ){
                    bool allOrfConserve = true;
                    bool allOrfLost = true;
                    for( std::vector<Transcript>::iterator it=it0->second[i].getTranscripts().begin();
                         it!=it0->second[i].getTranscripts().end(); ++it ){
                        if( it->getCdsVector().size() > 0 ){
                            if( it->getIfOrfShift() ){
                                allOrfConserve = false;
                            }else{
                                allOrfLost=false;
                            }
                        }
                    }
                    for( std::vector<Transcript>::iterator it=it0->second[j].getTranscripts().begin();
                         it!=it0->second[j].getTranscripts().end(); ++it ){
                        if( it->getCdsVector().size() > 0 ){
                            if( it->getIfOrfShift() ){
                                allOrfConserve = false;
                            }else{
                                allOrfLost=false;
                            }
                        }
                    }
                    if( allOrfConserve || allOrfLost ){
                        toRemove[it0->first].insert(i);
                    }else if( exactlySame( it0->second[j], it0->second[i]) ) {
                        toRemove[it0->first].insert(i);
                    }else{
                        mergeToFirstOne(it0->second[j], it0->second[i]);
                        toRemove[it0->first].insert(i);
                    }
                }
            }
        }
    }
}
