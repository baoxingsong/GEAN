/*
 * =====================================================================================
 *
 *       Filename:  InputParser.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 01:13:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include <iostream>
#include "InputParser.h"

InputParser::InputParser (int &argc, char **argv){
    for (int i=1; i < argc; ++i)
        this->tokens.push_back(std::string(argv[i]));
}

const std::string InputParser::getCmdOption( std::string &option) {
    std::vector<std::string>::iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
        return *itr;
    }
    return "";
}

std::string InputParser::getCmdOption( const char* o) {
    std::string option = o;
    return getCmdOption(option);
}

bool InputParser::cmdOptionExists(std::string &option) {
    return std::find(this->tokens.begin(), this->tokens.end(), option)
           != this->tokens.end();
}

bool InputParser::cmdOptionExists( const char* o){
    std::string option = o;
    return cmdOptionExists(option);
}

void usage( ){
    std::string progName = "gean";
    std::cout << "Program " << progName << std::endl <<
    "Usage: "<<progName<<" <command> [options]"<< std::endl <<
    "Commands:"<< std::endl <<
        " -- variant calling:" << std::endl <<
        "    pseudogeno  create pseudo genome sequence" << std::endl <<
        "    lift        transform coordinate to another accession" << std::endl <<
        "    revlift     transform coordinate of another accession to reference" << std::endl <<
        "    liftgff     transform all the GFF/GTF coordinates" << std::endl <<
        "    revliftgff  transform all the GFF/GTF coordinates back to reference" << std::endl <<
        "    reanva      update variants records for functional annotation" << std::endl <<
        "    gff2seq     get the protein/CDS/gene sequence of GFF/GTF file" << std::endl<<
        "    annowgr     annotate re-sequenced genome" << std::endl <<
        "    randomVar   assign a random position for each variant" << std::endl << std::endl <<
//        "annotationLiftOver             transform the reference GFF/GTF coordinate to the coordinate of another accession," << std::endl <<
//        "                               complementing with trying to keep the ORF/splice sites and complete as possible by " << std::endl <<
//        "                               genome sequence alignment" << std::endl <<
//        "annotationLiftOverAndOrth      transform the reference GFF/GTF coordinate to the coordinate of another accession," << std::endl <<
//        "                               complementing with trying to keep the ORF/splice sites and complete as possible by " << std::endl <<
//        "                               genome sequence alignment and then complementing with aligning CDS sequence of " << std::endl <<
//        "                               reference to the genome sequence of target line and then complementing with " << std::endl <<
//        "                               aligning protein sequence of reference to the genome sequence of target line" << std::endl <<
        " -- whole genome wide MSA:" << std::endl <<
        "    premsa      cut the whole genome sequence into fragments" << std::endl <<
        "    msatosdi    generate sdi files from MSA results" <<  std::endl << std::endl <<
        " -- de novo assembly genome:" << std::endl <<
        "    transgff    trans reference gff/gtf to de novo assembly genome" <<  std::endl <<
        "    spltogff    trans reference gff/gtf to de novo assembly genome using splice alignment sam file" <<  std::endl <<
        "    purifygff   purify the result from transgff" <<  std::endl <<
        "    sinsyn      (beta version) purify by keeping syntenic genes priorly and single copy genes (for inner species) " <<  std::endl <<
        //"    sinsyn2     (beta version) purify by keeping syntenic genes priorly and single copy genes (for inter species) " <<  std::endl <<
        //"    quotasyn    (beta version) quota syntenic blocks  " <<  std::endl <<
        "    orf         keep only ORF conserved genes " <<  std::endl <<
        "    varcall     variant calling for de novo genome sequence " <<  std::endl;
}
