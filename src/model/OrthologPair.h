//
// Created by Baoxing Song on 2019-03-13.
//

#ifndef GEAN_OrthologPair_H
#define GEAN_OrthologPair_H

#include <string>
#include "STRAND.h"



class OrthologPair {
private:
    int queryIndex;
    uint32_t refStartPos; // start reference coordinate
    uint32_t refEndPos; // end reference coordinate
    uint32_t refMiddlePos; // middle reference coordinate
    uint32_t queryStartPos; // start position of assembly/query
    uint32_t queryEndPos; // end position of assembly/query
    uint32_t queryMiddlePos; // middle position of assembly/query
    double score; // alignment score using as graph edge
    STRAND strand; // positive means same strand and positive means different strand

public:
    OrthologPair(const int & queryIndex, const uint32_t & refStartPos, const uint32_t & refEndPos,
                 const uint32_t & queryStartPos, const uint32_t & queryEndPos, const double & score,
                 const STRAND & strand) ;

    uint32_t getRefStartPos() const;
    void setRefStartPos(uint32_t refStartPos);
    uint32_t getRefEndPos() const;
    void setRefEndPos(uint32_t refEndPos);
    uint32_t getRefMiddlePos() const;
    void setRefMiddlePos(uint32_t refMiddlePos);
    uint32_t getQueryStartPos() const;
    void setQueryStartPos(uint32_t queryStartPos);
    uint32_t getQueryEndPos() const;
    void setQueryEndPos(uint32_t queryEndPos);\
    uint32_t getQueryMiddlePos() const;\
    void setQueryMiddlePos(uint32_t queryMiddlePos);
    double getScore() const;
    void setScore(double score);
    STRAND getStrand() const;
    void setStrand(STRAND &strand);

    int getQueryIndex() const;

    void setQueryIndex(int queryIndex);

    bool operator<( const OrthologPair& orthologPair ) const{
            if( refMiddlePos < orthologPair.refMiddlePos){
                    return true;
            }else if(refMiddlePos == orthologPair.refMiddlePos && queryMiddlePos<orthologPair.queryMiddlePos ){
                    return true;
            }
            return false;
    }
    bool operator>(const OrthologPair& orthologPair )const {
        if( refMiddlePos > orthologPair.refMiddlePos){
            return true;
        }else if(refMiddlePos == orthologPair.refMiddlePos && queryMiddlePos>orthologPair.queryMiddlePos ){
            return true;
        }
        return false;
    }
    bool operator==(const OrthologPair& orthologPair ) const{
        if( refMiddlePos == orthologPair.refMiddlePos && queryMiddlePos==orthologPair.queryMiddlePos ){
                return true;
        }
        return false;
    }
    bool operator!=(const OrthologPair& orthologPair ) const{
        if( refMiddlePos == orthologPair.refMiddlePos && queryMiddlePos==orthologPair.queryMiddlePos ){
                return false;
        }
        return true;
    }
};

#endif //GEAN_OrthologPair_H
