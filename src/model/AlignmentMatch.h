//
// Created by song on 8/4/18.
//
//
// the strand of database is always POSITIVE
//
//
#ifndef ANNOTATIONLIFTOVER_ALIGNMENTMATCH_H
#define ANNOTATIONLIFTOVER_ALIGNMENTMATCH_H

#include "./Range.h"
class AlignmentMatch {
private:
    Range query;
    Range database;
    size_t windowSize=0;
public:
    size_t getWindowSize() const;

    void setWindowSize(size_t windowSize);

public:
    AlignmentMatch(const std::string & queryChr, const size_t & queryStart, const size_t & queryEnd, const STRAND & queryStrand,
            const std::string & databaseChr, const size_t & databaseStart, const size_t & databaseEnd/*, const STRAND & databaseStrand*/ );
    AlignmentMatch(const std::string & queryChr, const size_t & queryStart, const size_t & queryEnd, const STRAND & queryStrand,
                   const std::string & databaseChr, const size_t & databaseStart, const size_t & databaseEnd, const size_t & _windowSize);
    const Range &getQuery() const;
    void setQuery(const Range &query);
    const Range &getDatabase() const;
    void setDatabase(const Range &database);

    const std::string & getQueryChr();
    void setQueryChr(const std::string & queryChr);
    const size_t & getQueryStart() const;
    void setQueryStart(const size_t & queryStart);
    const size_t & getQueryEnd() const;
    void setQueryEnd( const size_t & queryEnd );
    const STRAND & getQueryStrand() const;
    void setQueryStrand(const STRAND & queryStrand);
    const std::string & getDatabaseChr() const;
    void setDatabaseChr(const std::string & databaseChr);
    const size_t & getDatabaseStart() const;
    void setDatabaseStart(const size_t & databaseStart);
    const size_t & getDatabaseEnd() const;
    void setDatabaseEnd(const size_t & databaseEnd);
    /*
    const STRAND & getDatabaseStrand() const;
    void setDatabaseStrand(const STRAND & strand);*/
};


#endif //ANNOTATIONLIFTOVER_ALIGNMENTMATCH_H
