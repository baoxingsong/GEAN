//
// Created by baoxing on 10/10/17.
//

#include "myCoordinateLiftOver.h"



int myGetChangedFromReference( std::string& sdiFile, std::string& chromosome, int& position, std::string& vcfFix, const std::string& referenceGenomeFastaFile ){
    std::map<std::string, Fasta> referenceSequences;
    readFastaFile(referenceGenomeFastaFile, referenceSequences);

    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix, referenceSequences);
    return getChangedFromBasement( chromosome, position, sdiMaps );
}

int myGetReferenceFromChanged( std::string& sdiFile, std::string& chromosome, int& position, std::string& vcfFix, const std::string& referenceGenomeFastaFile ){
    std::map<std::string, Fasta> referenceSequences;
    readFastaFile(referenceGenomeFastaFile, referenceSequences);

    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix, referenceSequences);
    return getBasementFromChanged( chromosome, position, sdiMaps );
}

void myGffCoordinateLiftOver( std::string& sdiFile, std::string& gffFile, std::string& outputFile, std::string& vcfFix, const std::string& referenceGenomeFastaFile ){
    std::map<std::string, Fasta> referenceSequences;
    readFastaFile(referenceGenomeFastaFile, referenceSequences);

    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix, referenceSequences);

    std::ofstream ofile;
    ofile.open(outputFile);
    std::ifstream infile(gffFile);

    std::regex reg("^([\\s\\S]+?)\t([\\s\\S]+?)\t([\\s\\S]+?)\t(\\S+?)\t(\\S+?)\t([\\s\\S]+)$");
    std::string line;
    while (std::getline(infile, line)){
        if( line.size() > 8 ){
            std::smatch match;
            std::regex_search(line, match, reg);
            if( match.size()>5 ){
                int start = std::stoi(match[4]);
                int end = std::stoi(match[5]);
                int startL = getChangedFromBasement( match[1], start, sdiMaps );
                int endL = getChangedFromBasement( match[1], end, sdiMaps );
                ofile << match[1] << "\t" << match[2] << "\t" << match[3] <<
                      "\t" << startL << "\t" << endL <<"\t" << match[6] << std::endl;
                continue;
            }
        }
        ofile << line << std::endl;
    }
    ofile.close();
}

void myRevGffCoordinateLiftOver( std::string& sdiFile, std::string& gffFile, std::string& outputFile, std::string& vcfFix, const std::string& referenceGenomeFastaFile ){
    std::map<std::string, Fasta> referenceSequences;
    readFastaFile(referenceGenomeFastaFile, referenceSequences);

    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix, referenceSequences);

    std::ofstream ofile;
    ofile.open(outputFile);
    std::ifstream infile(gffFile);

    std::regex reg("^([\\s\\S]*?)\t([\\s\\S]*?)\t([\\s\\S]*?)\t(\\S*?)\t(\\S*?)\t([\\s\\S]+)$");
    std::string line;
    while (std::getline(infile, line)){
        if( line.size() > 8 ){
            std::smatch match;
            std::regex_search(line, match, reg);
            if( match.size()>5 ){
                int start = std::stoi(match[4]);
                int end = std::stoi(match[5]);
                int startL = getBasementFromChanged( match[1], start, sdiMaps );
                int endL = getBasementFromChanged( match[1], end, sdiMaps );
                ofile << match[1] << "\t" << match[2] << "\t" << match[3] <<
                      "\t" << startL << "\t" << endL <<"\t" << match[6] << std::endl;
                continue;
            }
        }
        ofile << line << std::endl;
    }
    ofile.close();
}