//
// Created by song on 8/4/18.
//

/**
 * the strand of database is always positive
 * */

#include "NucmerRead.h"
void nucmerRead(const std::string & filePath, std::map<std::string, std::vector<AlignmentMatch>> & alignmentMatchsMap){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening file " << filePath << std::endl;
        exit (1);
    }
    std::string line;
    char delim = '\t';
    std::vector<std::string> elems;
    std::string databaseChr;
    size_t databaseStart;
    size_t databaseEnd;
    std::string queryChr;
    size_t queryStart;
    size_t queryEnd;
    STRAND positive=POSITIVE;
    STRAND negative=NEGATIVE;
    size_t windowSize = 0;
    while (std::getline(infile, line)){
        elems.clear();
        split(line, delim, elems);
        databaseStart=stoi(elems[1]);
        databaseEnd=stoi(elems[2]);
        databaseChr=elems[0];
        queryStart=stoi(elems[4]);
        queryEnd=stoi(elems[5]);
        queryChr=elems[3];
        if( elems.size()>6 ){
            windowSize=stoi(elems[6]);
        }else{
            windowSize=0;
        }
        if( alignmentMatchsMap.find(databaseChr) == alignmentMatchsMap.end() ){
            alignmentMatchsMap[databaseChr] = std::vector<AlignmentMatch>();
        }
        if( databaseStart<databaseEnd ){
            if( queryStart<queryEnd ){
                AlignmentMatch alignmentMatch(queryChr, queryStart, queryEnd, positive, databaseChr, databaseStart, databaseEnd, windowSize);
                alignmentMatchsMap[databaseChr].push_back(alignmentMatch);
            }else{
                AlignmentMatch alignmentMatch(queryChr, queryEnd, queryStart, negative, databaseChr, databaseStart, databaseEnd, windowSize);
                alignmentMatchsMap[databaseChr].push_back(alignmentMatch);
            }
        }else{
            if( queryStart<queryEnd ){
                AlignmentMatch alignmentMatch(queryChr, queryStart, queryEnd, negative, databaseChr, databaseEnd, databaseStart, windowSize);
                alignmentMatchsMap[databaseChr].push_back(alignmentMatch);
            }else{
                AlignmentMatch alignmentMatch(queryChr, queryEnd, queryStart, positive, databaseChr, databaseEnd, databaseStart, windowSize);
                alignmentMatchsMap[databaseChr].push_back(alignmentMatch);
            }
        }
    }
}
