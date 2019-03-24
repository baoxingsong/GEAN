//
// Created by Baoxing Song on 2019-03-13.
//

#ifndef GEAN_OrthologPair_H
#define GEAN_OrthologPair_H

#include <string>
#include "STRAND.h"


class OrthologPair2 {
private:
    std::string refChr;
    std::string queryChr;
    int queryIndex;
    uint32_t refStartPos; // start reference coordinate
    uint32_t refEndPos; // end reference coordinate
    uint32_t refMiddlePos; // middle reference coordinate
    uint32_t queryStartPos; // start position of assembly/query
    uint32_t queryEndPos; // end position of assembly/query
    uint32_t queryMiddlePos; // middle position of assembly/query
    double score; // alignment score using as graph edge
    STRAND strand; // positive means same strand and positive means different strand
    int refId;
    int queryId;
    std::string referenceGeneName;
    std::string queryGeneName;
public:
    OrthologPair2();
    OrthologPair2( const OrthologPair2 & orthologPair2 );
    OrthologPair2(const std::string &refChr, const std::string &queryChr, const int & queryIndex, const uint32_t & refStartPos, const uint32_t & refEndPos,
                 const uint32_t & queryStartPos, const uint32_t & queryEndPos, const double & score,
                 const STRAND & strand, const int & refId, const int & queryId, const std::string & referenceGeneName,
                 const std::string & queryGeneName);
    OrthologPair2(const int & queryIndex, const uint32_t & refStartPos, const uint32_t & refEndPos,
                 const uint32_t & queryStartPos, const uint32_t & queryEndPos, const double & score,
                 const STRAND & strand) ;
    const std::string &getRefChr() const;
    void setRefChr(const std::string &refChr);
    const std::string &getQueryChr() const;
    void setQueryChr(const std::string &queryChr);

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
    int getRefId() const;
    void setRefId(int refId);
    int getQueryId() const;
    void setQueryId(int queryId);

    const std::string &getReferenceGeneName() const;

    void setReferenceGeneName(const std::string &referenceGeneName);

    const std::string &getQueryGeneName() const;

    void setQueryGeneName(const std::string &queryGeneName);

    bool operator<( const OrthologPair2& OrthologPair2 ) const{
            if( refMiddlePos < OrthologPair2.refMiddlePos){
                    return true;
            }else if(refMiddlePos == OrthologPair2.refMiddlePos && queryMiddlePos<OrthologPair2.queryMiddlePos ){
                    return true;
            }
            return false;
    }
    bool operator>(const OrthologPair2& OrthologPair2 )const {
        if( refMiddlePos > OrthologPair2.refMiddlePos){
            return true;
        }else if(refMiddlePos == OrthologPair2.refMiddlePos && queryMiddlePos>OrthologPair2.queryMiddlePos ){
            return true;
        }
        return false;
    }
    bool operator==(const OrthologPair2& OrthologPair2 ) const{
        if( refMiddlePos == OrthologPair2.refMiddlePos && queryMiddlePos==OrthologPair2.queryMiddlePos ){
                return true;
        }
        return false;
    }
    bool operator!=(const OrthologPair2& OrthologPair2 ) const{
        if( refMiddlePos == OrthologPair2.refMiddlePos && queryMiddlePos==OrthologPair2.queryMiddlePos ){
                return false;
        }
        return true;
    }
};



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
