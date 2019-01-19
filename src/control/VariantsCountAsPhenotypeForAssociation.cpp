//
// Created by baoxing on 10/10/17.
//

#include "VariantsCountAsPhenotypeForAssociation.h"


void countNumberOfTwoneighborSNP( std::string& sdiFile, std::string & outputPrefix, int & rangeLength, std::string& vcfFix, std::map<std::string, std::string>& parameters, const std::string& referenceGenomeFastaFile){
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    std::set<std::string>& legalNasString = nucleotideCodeSubstitutionMatrix.getLegalNasString();

    std::map<std::string, Fasta> referenceSequences;
    readFastaFile(referenceGenomeFastaFile, referenceSequences);

    std::ofstream ofilealldouble;
    ofilealldouble.open(outputPrefix+".alldouble.sdi");

    std::ofstream ofileaacc;
    ofileaacc.open(outputPrefix+".aacc.sdi");

    std::ofstream ofileaacg;
    ofileaacg.open(outputPrefix+".aacg.sdi");

    std::ofstream ofileaact;
    ofileaact.open(outputPrefix+".aact.sdi");

    std::ofstream ofileaagc;
    ofileaagc.open(outputPrefix+".aagc.sdi");

    std::ofstream ofileaagg;
    ofileaagg.open(outputPrefix+".aagg.sdi");

    std::ofstream ofileaagt;
    ofileaagt.open(outputPrefix+".aagt.sdi");

    std::ofstream ofileaatc;
    ofileaatc.open(outputPrefix+".aatc.sdi");

    std::ofstream ofileaatg;
    ofileaatg.open(outputPrefix+".aatg.sdi");

    std::ofstream ofileaatt;
    ofileaatt.open(outputPrefix+".aatt.sdi");

    std::ofstream ofileacca;
    ofileacca.open(outputPrefix+".acca.sdi");

    std::ofstream ofileaccg;
    ofileaccg.open(outputPrefix+".accg.sdi");

    std::ofstream ofileacct;
    ofileacct.open(outputPrefix+".acct.sdi");

    std::ofstream ofileacga;
    ofileacga.open(outputPrefix+".acga.sdi");

    std::ofstream ofileacgg;
    ofileacgg.open(outputPrefix+".acgg.sdi");

    std::ofstream ofileacgt;
    ofileacgt.open(outputPrefix+".acgt.sdi");

    std::ofstream ofileacta;
    ofileacta.open(outputPrefix+".acta.sdi");

    std::ofstream ofileactg;
    ofileactg.open(outputPrefix+".actg.sdi");

    std::ofstream ofileagca;
    ofileagca.open(outputPrefix+".agca.sdi");

    std::ofstream ofileagcc;
    ofileagcc.open(outputPrefix+".agcc.sdi");

    std::ofstream ofileagct;
    ofileagct.open(outputPrefix+".agct.sdi");

    std::ofstream ofileagga;
    ofileagga.open(outputPrefix+".agga.sdi");

    std::ofstream ofileaggc;
    ofileaggc.open(outputPrefix+".aggc.sdi");

    std::ofstream ofileagta;
    ofileagta.open(outputPrefix+".agta.sdi");

    std::ofstream ofileagtc;
    ofileagtc.open(outputPrefix+".agtc.sdi");

    std::ofstream ofileatca;
    ofileatca.open(outputPrefix+".atca.sdi");

    std::ofstream ofileatcc;
    ofileatcc.open(outputPrefix+".atcc.sdi");

    std::ofstream ofileatcg;
    ofileatcg.open(outputPrefix+".atcg.sdi");

    std::ofstream ofileatga;
    ofileatga.open(outputPrefix+".atga.sdi");

    std::ofstream ofileatgc;
    ofileatgc.open(outputPrefix+".atgc.sdi");

    std::ofstream ofileatta;
    ofileatta.open(outputPrefix+".atta.sdi");

    std::ofstream ofilecagc;
    ofilecagc.open(outputPrefix+".cagc.sdi");

    std::ofstream ofilecagg;
    ofilecagg.open(outputPrefix+".cagg.sdi");

    std::ofstream ofilecatc;
    ofilecatc.open(outputPrefix+".catc.sdi");

    std::ofstream ofilecatg;
    ofilecatg.open(outputPrefix+".catg.sdi");

    std::ofstream ofileccga;
    ofileccga.open(outputPrefix+".ccga.sdi");

    std::ofstream ofileccgg;
    ofileccgg.open(outputPrefix+".ccgg.sdi");

    std::ofstream ofileccta;
    ofileccta.open(outputPrefix+".ccta.sdi");

    std::ofstream ofilecgga;
    ofilecgga.open(outputPrefix+".cgga.sdi");

    std::ofstream ofilecggc;
    ofilecggc.open(outputPrefix+".cggc.sdi");

    std::ofstream ofilecgta;
    ofilecgta.open(outputPrefix+".cgta.sdi");

    std::ofstream ofilegatc;
    ofilegatc.open(outputPrefix+".gatc.sdi");

    std::ofstream ofilegcta;
    ofilegcta.open(outputPrefix+".gcta.sdi");

    std::ofstream ofileunclassified;
    ofileunclassified.open(outputPrefix+".unclassified.sdi");


    std::ofstream ofileRaatt;
    ofileRaatt.open(outputPrefix+".Raatt.sdi");
    std::ofstream ofileRacgt;
    ofileRacgt.open(outputPrefix+".Racgt.sdi");
    std::ofstream ofileRatat;
    ofileRatat.open(outputPrefix+".Ratat.sdi");
    std::ofstream ofileRagct;
    ofileRagct.open(outputPrefix+".Ragct.sdi");
    std::ofstream ofileRcatg;
    ofileRcatg.open(outputPrefix+".Rcatg.sdi");
    std::ofstream ofileRccgg;
    ofileRccgg.open(outputPrefix+".Rccgg.sdi");
    std::ofstream ofileRcgcg;
    ofileRcgcg.open(outputPrefix+".Rcgcg.sdi");
    std::ofstream ofileRgatc;
    ofileRgatc.open(outputPrefix+".Rgatc.sdi");
    std::ofstream ofileRgcgc;
    ofileRgcgc.open(outputPrefix+".Rgcgc.sdi");
    std::ofstream ofileRtata;
    ofileRtata.open(outputPrefix+".Rtata.sdi");

    std::map<std::string, std::vector<Variant> > sdiMaps;
    //std::cout << "495" <<std::endl;
    readSdiFile(sdiFile, sdiMaps, vcfFix, referenceSequences);
    //std::cout << "497" <<std::endl;
    for( std::map<std::string, std::vector<Variant> >::iterator it=sdiMaps.begin();
            it!=sdiMaps.end(); ++it){
        size_t thisChromosomeRecordsNumber = sdiMaps[it->first].size();
        size_t i=1;
        std::vector<size_t> wantedIds;
        //  std::cout << "501" <<std::endl;
        if( sdiMaps[it->first][i].getPosition() == sdiMaps[it->first][i-1].getPosition()+1
            && sdiMaps[it->first][i].getChanginglength()==0 &&
            sdiMaps[it->first][i-1].getChanginglength() == 0) {
            if( legalNasString.find(sdiMaps[it->first][i].getReference())!= legalNasString.end() &&
                legalNasString.find(sdiMaps[it->first][i-1].getReference())!= legalNasString.end()
                && legalNasString.find(sdiMaps[it->first][i].getAlternative())!= legalNasString.end() &&
                legalNasString.find(sdiMaps[it->first][i-1].getAlternative())!= legalNasString.end() ){
                bool addIt = true;
                for( int j = 0; j <= rangeLength; ++j){
                    if (sdiMaps[it->first][i + 1].overlap(sdiMaps[it->first][i].getPosition() + j)) { //
                        addIt = false;
                    }
                }
                if( addIt ){
                    wantedIds.push_back(i);
                }
            }
        }
        //std::cout << "519" <<std::endl;
        for( i=2; i < thisChromosomeRecordsNumber-1; ++i ){
            if( sdiMaps[it->first][i].getPosition() == sdiMaps[it->first][i-1].getPosition()+1
                && sdiMaps[it->first][i].getChanginglength()==0 &&
                sdiMaps[it->first][i-1].getChanginglength() == 0){
                if( legalNasString.find(sdiMaps[it->first][i].getReference())!= legalNasString.end() &&
                    legalNasString.find(sdiMaps[it->first][i-1].getReference())!= legalNasString.end()
                    && legalNasString.find(sdiMaps[it->first][i].getAlternative())!= legalNasString.end() &&
                    legalNasString.find(sdiMaps[it->first][i-1].getAlternative())!= legalNasString.end()) {
                    bool addIt = true;
                    for( int j = 0; j <= rangeLength; ++j){
                        if (sdiMaps[it->first][i + 1].overlap(sdiMaps[it->first][i].getPosition() + j)) { //
                            addIt = false;
                        }
                        if (sdiMaps[it->first][i - 2].overlap(sdiMaps[it->first][i - 1].getPosition() - j)) {
                            addIt = false;
                        }
                    }
                    if( addIt ){
                        wantedIds.push_back(i);
                    }
                }
            }
        }
        //std::cout << "543" <<std::endl;
        i = thisChromosomeRecordsNumber-1;
        if( sdiMaps[it->first][i].getPosition() == sdiMaps[it->first][i-1].getPosition()+1
            && sdiMaps[it->first][i].getChanginglength()==0 &&
            sdiMaps[it->first][i-1].getChanginglength() == 0) {
            if( legalNasString.find(sdiMaps[it->first][i].getReference())!= legalNasString.end() &&
                legalNasString.find(sdiMaps[it->first][i-1].getReference())!= legalNasString.end()
                && legalNasString.find(sdiMaps[it->first][i].getAlternative())!= legalNasString.end() &&
                legalNasString.find(sdiMaps[it->first][i-1].getAlternative())!= legalNasString.end()){
                bool addIt = true;
                for( int j = 0; j <= rangeLength; ++j){
                    if (sdiMaps[it->first][i - 2].overlap(sdiMaps[it->first][i - 1].getPosition() - j)) {
                        addIt = false;
                    }
                }
                if( addIt ){
                    wantedIds.push_back(i);
                }
            }
        }
        for( std::vector<size_t>::iterator itj=wantedIds.begin();
             itj!=wantedIds.end(); ++itj){
            size_t j = (*itj);
            printSdiln(ofilealldouble, sdiMaps[it->first][j - 1]);
            printSdiln(ofilealldouble, sdiMaps[it->first][j]);
            //  care about both the reference and alternative
            if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileaacc, sdiMaps[it->first][j - 1]);
                printSdiln(ofileaacc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileaacg, sdiMaps[it->first][j - 1]);
                printSdiln(ofileaacg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileaact, sdiMaps[it->first][j - 1]);
                printSdiln(ofileaact, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileaagc, sdiMaps[it->first][j - 1]);
                printSdiln(ofileaagc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileaagg, sdiMaps[it->first][j - 1]);
                printSdiln(ofileaagg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileaagt, sdiMaps[it->first][j - 1]);
                printSdiln(ofileaagt, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileaatc, sdiMaps[it->first][j - 1]);
                printSdiln(ofileaatc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileaatg, sdiMaps[it->first][j - 1]);
                printSdiln(ofileaatg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") )
                    ){
                printSdiln(ofileaatt, sdiMaps[it->first][j - 1]);
                printSdiln(ofileaatt, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileacca, sdiMaps[it->first][j - 1]);
                printSdiln(ofileacca, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileaccg, sdiMaps[it->first][j - 1]);
                printSdiln(ofileaccg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileacct, sdiMaps[it->first][j - 1]);
                printSdiln(ofileacct, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileacga, sdiMaps[it->first][j - 1]);
                printSdiln(ofileacga, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileacgg, sdiMaps[it->first][j - 1]);
                printSdiln(ofileacgg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") )
                    ){
                printSdiln(ofileacgt, sdiMaps[it->first][j - 1]);
                printSdiln(ofileacgt, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileacta, sdiMaps[it->first][j - 1]);
                printSdiln(ofileacta, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileactg, sdiMaps[it->first][j - 1]);
                printSdiln(ofileactg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileagca, sdiMaps[it->first][j - 1]);
                printSdiln(ofileagca, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileagcc, sdiMaps[it->first][j - 1]);
                printSdiln(ofileagcc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdiln(ofileagct, sdiMaps[it->first][j - 1]);
                printSdiln(ofileagct, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileagga, sdiMaps[it->first][j - 1]);
                printSdiln(ofileagga, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileaggc, sdiMaps[it->first][j - 1]);
                printSdiln(ofileaggc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileagta, sdiMaps[it->first][j - 1]);
                printSdiln(ofileagta, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileagtc, sdiMaps[it->first][j - 1]);
                printSdiln(ofileagtc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileatca, sdiMaps[it->first][j - 1]);
                printSdiln(ofileatca, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileatcc, sdiMaps[it->first][j - 1]);
                printSdiln(ofileatcc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileatcg, sdiMaps[it->first][j - 1]);
                printSdiln(ofileatcg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileatga, sdiMaps[it->first][j - 1]);
                printSdiln(ofileatga, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileatgc, sdiMaps[it->first][j - 1]);
                printSdiln(ofileatgc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdiln(ofileatta, sdiMaps[it->first][j - 1]);
                printSdiln(ofileatta, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdiln(ofilecagc, sdiMaps[it->first][j - 1]);
                printSdiln(ofilecagc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdiln(ofilecagg, sdiMaps[it->first][j - 1]);
                printSdiln(ofilecagg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdiln(ofilecatc, sdiMaps[it->first][j - 1]);
                printSdiln(ofilecatc, sdiMaps[it->first][j]);
            } else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                       ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") )
                    ){
                printSdiln(ofilecatg, sdiMaps[it->first][j - 1]);
                printSdiln(ofilecatg, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdiln(ofileccga, sdiMaps[it->first][j - 1]);
                printSdiln(ofileccga, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") )
                    ){
                printSdiln(ofileccgg, sdiMaps[it->first][j - 1]);
                printSdiln(ofileccgg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdiln(ofileccta, sdiMaps[it->first][j - 1]);
                printSdiln(ofileccta, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdiln(ofilecgga, sdiMaps[it->first][j - 1]);
                printSdiln(ofilecgga, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdiln(ofilecggc, sdiMaps[it->first][j - 1]);
                printSdiln(ofilecggc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdiln(ofilecgta, sdiMaps[it->first][j - 1]);
                printSdiln(ofilecgta, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") )
                    ){
                printSdiln(ofilegatc, sdiMaps[it->first][j - 1]);
                printSdiln(ofilegatc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") )
                    ){
                printSdiln(ofilegcta, sdiMaps[it->first][j - 1]);
                printSdiln(ofilegcta, sdiMaps[it->first][j]);
            }else{
                printSdiln(ofileunclassified, sdiMaps[it->first][j - 1]);
                printSdiln(ofileunclassified, sdiMaps[it->first][j]);
            }

            //  only care about what is the reference
            if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A")  ) ||

                ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") ) ){
                printSdiln(ofileRaatt, sdiMaps[it->first][j - 1]);
                printSdiln(ofileRaatt, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C")  ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") ) ){
                printSdiln(ofileRacgt, sdiMaps[it->first][j - 1]);
                printSdiln(ofileRacgt, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T")  ) ){
                printSdiln(ofileRatat, sdiMaps[it->first][j - 1]);
                printSdiln(ofileRatat, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G")  ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") ) ){
                printSdiln(ofileRagct, sdiMaps[it->first][j - 1]);
                printSdiln(ofileRagct, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A")  ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") ) ){
                printSdiln(ofileRcatg, sdiMaps[it->first][j - 1]);
                printSdiln(ofileRcatg, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C")  ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") ) ){
                printSdiln(ofileRccgg, sdiMaps[it->first][j - 1]);
                printSdiln(ofileRccgg, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G")  )  ){
                printSdiln(ofileRcgcg, sdiMaps[it->first][j - 1]);
                printSdiln(ofileRcgcg, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A")  ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") ) ){
                printSdiln(ofileRgatc, sdiMaps[it->first][j - 1]);
                printSdiln(ofileRgatc, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C")  ) ){
                printSdiln(ofileRgcgc, sdiMaps[it->first][j - 1]);
                printSdiln(ofileRgcgc, sdiMaps[it->first][j]);
            }else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A")  ) ){
                printSdiln(ofileRtata, sdiMaps[it->first][j - 1]);
                printSdiln(ofileRtata, sdiMaps[it->first][j]);
            }
        }
    }
    ofilealldouble.close();


    ofileaacc.close();
    ofileaacg.close();
    ofileaact.close();
    ofileaagc.close();
    ofileaagg.close();
    ofileaagt.close();
    ofileaatc.close();
    ofileaatg.close();
    ofileaatt.close();
    ofileacca.close();
    ofileaccg.close();
    ofileacct.close();
    ofileacga.close();
    ofileacgg.close();
    ofileacgt.close();
    ofileacta.close();
    ofileactg.close();
    ofileagca.close();
    ofileagcc.close();
    ofileagct.close();
    ofileagga.close();
    ofileaggc.close();
    ofileagta.close();
    ofileagtc.close();
    ofileatca.close();
    ofileatcc.close();
    ofileatcg.close();
    ofileatga.close();
    ofileatgc.close();
    ofileatta.close();
    ofilecagc.close();
    ofilecagg.close();
    ofilecatc.close();
    ofilecatg.close();
    ofileccga.close();
    ofileccgg.close();
    ofileccta.close();
    ofilecgga.close();
    ofilecggc.close();
    ofilecgta.close();
    ofilegatc.close();
    ofilegcta.close();
    ofileunclassified.close();


    ofileRaatt.close();
    ofileRacgt.close();
    ofileRatat.close();
    ofileRagct.close();
    ofileRcatg.close();
    ofileRccgg.close();
    ofileRcgcg.close();
    ofileRgatc.close();
    ofileRgcgc.close();
    ofileRtata.close();

}


void countNumberSNPAndIndel( std::string& sdiFile, std::string & outputPrefix, int & rangeLength, std::string& vcfFix, std::map<std::string, std::string>& parameters, const std::string& referenceGenomeFastaFile) {

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    std::set<std::string> &legalNasString = nucleotideCodeSubstitutionMatrix.getLegalNasString();

    std::map<std::string, Fasta> referenceSequences;
    readFastaFile(referenceGenomeFastaFile, referenceSequences);

    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix, referenceSequences);

    std::ofstream ofileallsnps;
    ofileallsnps.open(outputPrefix+".allsnps.sdi");
    std::ofstream ofileallindels;
    ofileallindels.open(outputPrefix+".allindels.sdi");

    std::ofstream ofilelegalsnps;
    ofilelegalsnps.open(outputPrefix+".legalsnps.sdi");

    std::ofstream ofilelegalsnpsNoNearByIndel;
    ofilelegalsnpsNoNearByIndel.open(outputPrefix+".legalsnpsNoNearByIndel.sdi");

    std::ofstream ofilelegalsnpsIsolate;
    ofilelegalsnpsIsolate.open(outputPrefix+".legalsnpsIsolate.sdi");

    std::ofstream ofileillegalsnpsNoNearByIndel;
    ofileillegalsnpsNoNearByIndel.open(outputPrefix+".illegalsnpsNoNearByIndel.sdi");

    std::ofstream ofileillegalsnpsIsolate;
    ofileillegalsnpsIsolate.open(outputPrefix+".illegalsnpsIsolate.sdi");

    std::ofstream ofileillegalsnps;
    ofileillegalsnps.open(outputPrefix+".illegalsnps.sdi");

    for( std::map<std::string, std::vector<Variant> >::iterator it=sdiMaps.begin();
            it!=sdiMaps.end(); ++it) {
        size_t thisChromosomeRecordsNumber = sdiMaps[it->first].size();

        for( size_t i=0; i < thisChromosomeRecordsNumber; ++i ){
            if(  sdiMaps[it->first][i].getChanginglength()==0 && sdiMaps[it->first][i].getReference().size()==1 ){
                if( legalNasString.find(sdiMaps[it->first][i].getReference())!= legalNasString.end()
                    && legalNasString.find(sdiMaps[it->first][i].getAlternative())!= legalNasString.end() ) {
                    printSdiln(ofilelegalsnps, sdiMaps[it->first][i]);

                    //ofilelegalsnpsNoNearByIndel begin
                    bool addIt = true;
                    for( int j = 0; j <= rangeLength; ++j){
                        if ( i>0 && sdiMaps[it->first][i-1].getChanginglength()!=0 && sdiMaps[it->first][i-1].overlap(sdiMaps[it->first][i].getPosition() - j)) { //
                            addIt = false;
                        }
                        if (i<(thisChromosomeRecordsNumber-1) && sdiMaps[it->first][i+1].getChanginglength()!=0 && sdiMaps[it->first][i +1].overlap(sdiMaps[it->first][i].getPosition() + j)) {
                            addIt = false;
                        }
                    }
                    if( addIt ){
                        printSdiln(ofilelegalsnpsNoNearByIndel, sdiMaps[it->first][i]);
                    }//ofilelegalsnpsNoNearByIndel end


                    //ofilelegalsnpsIsolate begin
                    addIt = true;
                    for( int j = 0; j <= rangeLength; ++j){
                        if ( i>0  && sdiMaps[it->first][i-1].overlap(sdiMaps[it->first][i].getPosition() - j)) { //
                            addIt = false;
                        }
                        if (i<(thisChromosomeRecordsNumber-1) && sdiMaps[it->first][i +1].overlap(sdiMaps[it->first][i].getPosition() + j)) {
                            addIt = false;
                        }
                    }
                    if( addIt ){
                        printSdiln(ofilelegalsnpsIsolate, sdiMaps[it->first][i]);
                    }
                    //ofilelegalsnpsIsolate end

                }else{
                    printSdiln(ofileillegalsnps, sdiMaps[it->first][i]);

                    //ofileillegalsnpsNoNearByIndel begin
                    bool addIt = true;
                    for( int j = 0; j <= rangeLength; ++j){
                        if ( i>0 && sdiMaps[it->first][i-1].getChanginglength()!=0 && sdiMaps[it->first][i-1].overlap(sdiMaps[it->first][i].getPosition() - j)) { //
                            addIt = false;
                        }
                        if (i<(thisChromosomeRecordsNumber-1) && sdiMaps[it->first][i+1].getChanginglength()!=0 && sdiMaps[it->first][i +1].overlap(sdiMaps[it->first][i].getPosition() + j)) {
                            addIt = false;
                        }
                    }
                    if( addIt ){
                        printSdiln(ofileillegalsnpsNoNearByIndel, sdiMaps[it->first][i]);
                    }//ofileillegalsnpsNoNearByIndel end


                    //ofileillegalsnpsIsolate begin
                    addIt = true;
                    for( int j = 0; j <= rangeLength; ++j){
                        if ( i>0  && sdiMaps[it->first][i-1].overlap(sdiMaps[it->first][i].getPosition() - j)) { //
                            addIt = false;
                        }
                        if (i<(thisChromosomeRecordsNumber-1) && sdiMaps[it->first][i +1].overlap(sdiMaps[it->first][i].getPosition() + j)) {
                            addIt = false;
                        }
                    }
                    if( addIt ){
                        printSdiln(ofileillegalsnpsIsolate, sdiMaps[it->first][i]);
                    }
                    //ofileillegalsnpsIsolate end

                }
                printSdiln(ofileallsnps, sdiMaps[it->first][i]);


            }else{
                printSdiln(ofileallindels, sdiMaps[it->first][i]);
            }
        }
    }
    ofileallsnps.close();
    ofileallindels.close();
    ofilelegalsnps.close();
    ofileillegalsnps.close();
    ofilelegalsnpsIsolate.close();
    ofilelegalsnpsNoNearByIndel.close();
    ofileillegalsnpsNoNearByIndel.close();
    ofileillegalsnpsIsolate.close();
}


int ranDomItDeletion( Variant & obvervedVariant, std::map<std::string, std::vector<Variant> >& randomSdiMaps, std::string& chr, int& chrSize,
        std::map<std::string, Fasta> & refGenome ){
    bool notInserted = true; // the record is not inserted
    int timesTried = 0;
    while( notInserted ){

        int position = (rand() % chrSize) + 1;
        std::string alt = obvervedVariant.getAlternative();
        std::string ref = getSubsequence(refGenome, obvervedVariant.getChromosome(), position, position+obvervedVariant.getReference().length()-1);
        Variant variant(chr, position, ref, alt); // the (*it2) variant records but assigned with a random position
        if( variant.getChanginglength() <0 && variant.getPosition()-variant.getChanginglength() > chrSize+1  ){
            // if this is deletion and deletion is beyond the chromosome, it does not make sense, so avoid such records
        }else{
            bool hasOverLap = false;
            for( std::vector<Variant>::iterator it3=randomSdiMaps[chr].begin(); it3!=randomSdiMaps[chr].end(); ++it3 ){
                // if the new record conflicts with any of the previous records, re-do position random
                int start1 = it3->getPosition();
                int end1 = it3->getPosition() + abs(it3->getChanginglength());
                int start2 = variant.getPosition();
                int end2 = variant.getPosition() + abs(variant.getChanginglength());
                if(  end1 < start2 || end2 < start1  ){
                    continue; // for such case, there is definitly no confliction
                }
                if(
                        ( (it3->getPosition()<=variant.getPosition() && variant.getPosition()<=it3->getPosition()-it3->getChanginglength()) ||
                          ( variant.getPosition()<=it3->getPosition() && it3->getPosition()<=variant.getPosition()-variant.getChanginglength() ) ) ){
                    // deletion VS deletion
                    //two deletion could not neighbor with each other, else they should be marged as one single records
                    hasOverLap = true;
                    break;
                }
            }
            if( ! hasOverLap  ){
                randomSdiMaps[chr].push_back(variant); // insert the variants with random position into random records data structure
                //printSdiln(std::cout, variant);
                notInserted = false;
                return 0; // good
            }
        }
        timesTried++;
        if( timesTried>100000 ){
            //after 100000 fails, give up and redo everything
            return 1;
        }
    }
    return 0; // good
}

int ranDomItSnp( Variant & obvervedVariant, std::map<std::string, std::vector<Variant> >& randomSdiMaps, std::string& chr, int& chrSize,
                 std::map<std::string, Fasta> & refGenome ){
    bool notInserted = true; // the record is not inserted
    int timesTried = 0;
    while( notInserted ){
        int position = (rand() % chrSize) + 1;
        std::string alt = obvervedVariant.getAlternative();
        std::string ref = getSubsequence(refGenome, obvervedVariant.getChromosome(), position, position+obvervedVariant.getReference().length()-1);
        Variant variant(chr, position, ref, alt); // the (*it2) variant records but assigned with a random position
        if(  ( obvervedVariant.getReference().compare("A") == 0
            || obvervedVariant.getReference().compare("C") == 0
            || obvervedVariant.getReference().compare("G") == 0
            || obvervedVariant.getReference().compare("T") == 0) && ref.compare(obvervedVariant.getReference()) != 0   ){

        }else if ( ref.compare(alt) == 0 ) { // this is not a variants

        }else{
            bool hasOverLap = false;
            for( std::vector<Variant>::iterator it3=randomSdiMaps[chr].begin(); it3!=randomSdiMaps[chr].end(); ++it3 ){
                // if the new record conflicts with any of the previous records, re-do position random
                int start1 = it3->getPosition();
                int end1 = it3->getPosition() + abs(it3->getChanginglength());
                int start2 = variant.getPosition();
                int end2 = variant.getPosition() + abs(variant.getChanginglength());
                if(  end1 < start2 || end2 < start1  ){
                    continue; // for such case, there is definitly no confliction
                }
                if( it3->getChanginglength()==0 && variant.getPosition()==it3->getPosition() ){
                    //SNP VS SNP
                    //two SNP could not share the same position
                    hasOverLap = true;
                    break;
                }else if( it3->getChanginglength()<0  && it3->getPosition()<=variant.getPosition() && variant.getPosition() <it3->getPosition()-it3->getChanginglength() ){
                    //deletion VS SNP
                    hasOverLap = true;
                    break;
                }// there is no confliction between insertion and SNP
            }
            if( ! hasOverLap  ){
                randomSdiMaps[chr].push_back(variant); // insert the variants with random position into random records data structure
                notInserted = false;
                return 0; // good
            }
        }
        timesTried++;
        if( timesTried>1000000 ){
            //after 1000000 fails, give up and redo everything
            return 1;
        }
    }
    return 0; // good
}


int ranDomItInsertion( Variant & obvervedVariant, std::map<std::string, std::vector<Variant> >& randomSdiMaps, std::string& chr, int& chrSize ){
    bool notInserted = true; // the record is not inserted
    int timesTried = 0;
    while( notInserted ){
        int position = (rand() % chrSize) + 1;
        std::string ref = obvervedVariant.getReference();
        std::string alt = obvervedVariant.getAlternative();
        Variant variant(chr, position, ref, alt); // the (*it2) variant records but assigned with a random position
        //if( variant.getChanginglength() <0 && variant.getPosition()-variant.getChanginglength() > chrSize+1  ){
            // if this is deletion and deletion is beyond the chromosome, it does not make sense, so avoid such records
        //}else{
            bool hasOverLap = false;
            for( std::vector<Variant>::iterator it3=randomSdiMaps[chr].begin(); it3!=randomSdiMaps[chr].end(); ++it3 ){
                // if the new record conflicts with any of the previous records, re-do position random
                int start1 = it3->getPosition();
                int end1 = it3->getPosition() + abs(it3->getChanginglength());
                int start2 = variant.getPosition();
                int end2 = variant.getPosition() + abs(variant.getChanginglength());
                if(  end1 < start2 || end2 < start1  ){
                    continue; // for such case, there is definitly no confliction
                }

                if( it3->getChanginglength()>0 && variant.getChanginglength()>0 && it3->getPosition()==position ){
                    //insertion VS insertion
                    //two insertion could not at the same position, else they should be marged as one single records
                    hasOverLap = true;
                    break;
                }else if( it3->getChanginglength()<0 && variant.getPosition()>it3->getPosition() &&  variant.getPosition()< it3->getPosition()-it3->getChanginglength() ){
                    //deletion VS insertion 1
                    //could no insert at deleted sequence
                    hasOverLap = true;
                    break;
                }// there is no confliction between insertion and SNP
            }
            if( ! hasOverLap  ){
                randomSdiMaps[chr].push_back(variant); // insert the variants with random position into random records data structure
//                printSdiln(std::cout, variant);
                notInserted = false;
                return 0; // good
            }
        //}
        timesTried++;
        if( timesTried>100000 ){
            //after 100000 fails, give up and redo everything
            return 1;
        }
    }
    return 0; // good
}

int ranDomIt( Variant & obvervedVariant, std::map<std::string, std::vector<Variant> >& randomSdiMaps, std::string& chr, int& chrSize ){
    bool notInserted = true; // the record is not inserted
    int timesTried = 0;
    while( notInserted ){

        int position = (rand() % chrSize) + 1;
        std::string ref = obvervedVariant.getReference();
        std::string alt = obvervedVariant.getAlternative();
        Variant variant(chr, position, ref, alt); // the (*it2) variant records but assigned with a random position
        if( variant.getChanginglength() <0 && variant.getPosition()-variant.getChanginglength() > chrSize+1  ){
            // if this is deletion and deletion is beyond the chromosome, it does not make sense, so avoid such records
        }else{
            bool hasOverLap = false;
            for( std::vector<Variant>::iterator it3=randomSdiMaps[chr].begin(); it3!=randomSdiMaps[chr].end(); ++it3 ){
                // if the new record conflicts with any of the previous records, re-do position random
                int start1 = it3->getPosition();
                int end1 = it3->getPosition() + abs(it3->getChanginglength());
                int start2 = variant.getPosition();
                int end2 = variant.getPosition() + abs(variant.getChanginglength());
                if(  end1 < start2 || end2 < start1  ){
                    continue; // for such case, there is definitly no confliction
                }
                if( it3->getChanginglength()<0 && variant.getChanginglength()<0 &&
                    ( (it3->getPosition()<=variant.getPosition() && variant.getPosition()<=it3->getPosition()-it3->getChanginglength()) ||
                      ( variant.getPosition()<=it3->getPosition() && it3->getPosition()<=variant.getPosition()-variant.getChanginglength() ) ) ){
                    // deletion VS deletion
                    //two deletion could not neighbor with each other, else they should be marged as one single records
                    hasOverLap = true;
                    break;
                }else if( it3->getChanginglength()>0 && variant.getChanginglength()>0 && it3->getPosition()==position ){
                    //insertion VS insertion
                    //two insertion could not at the same position, else they should be marged as one single records
                    hasOverLap = true;
                    break;
                }else if( it3->getChanginglength()<0 && variant.getChanginglength()>0 && variant.getPosition()>it3->getPosition() &&  variant.getPosition()< it3->getPosition()-it3->getChanginglength() ){
                    //deletion VS insertion 1
                    //could no insert at deleted sequence
                    hasOverLap = true;
                    break;
                }else if( variant.getChanginglength()<0 && it3->getChanginglength()>0 && it3->getPosition()>variant.getPosition() &&  it3->getPosition()< variant.getPosition()-variant.getChanginglength() ){
                    //deletion VS insertion 2
                    hasOverLap = true;
                    break;
                }else if( variant.getChanginglength()==0 && it3->getChanginglength()==0 && variant.getPosition()==it3->getPosition() ){
                    //SNP VS SNP
                    //two SNP could not share the same position
                    hasOverLap = true;
                    break;
                }else if( variant.getChanginglength()<0 && it3->getChanginglength()==0 && variant.getPosition()<=it3->getPosition() && it3->getPosition() <variant.getPosition()-variant.getChanginglength() ){
                    //deletion VS SNP 1
                    //could no SNP at deleted sequence
                    hasOverLap = true;
                    break;
                }else if( it3->getChanginglength()<0 && variant.getChanginglength()==0 && it3->getPosition()<=variant.getPosition() && variant.getPosition() <it3->getPosition()-variant.getChanginglength() ){
                    //deletion VS SNP 2
                    hasOverLap = true;
                    break;
                }// there is no confliction between insertion and SNP
            }
            if( ! hasOverLap  ){
                randomSdiMaps[chr].push_back(variant); // insert the variants with random position into random records data structure
                //printSdiln(std::cout, variant);
                notInserted = false;
                return 0; // good
            }
        }
        timesTried++;
        if( timesTried>100000 ){
            //after 100000 fails, give up and redo everything
            return 1;
        }
    }
    return 0; // good
}

void generateRandomSdi( std::string& sdiFile, std::string & outputPrefix, std::string& vcfFix, const std::string& referenceGenomeFastaFile){

    std::map<std::string, Fasta> refGenome;
    readFastaFile(referenceGenomeFastaFile, refGenome);

    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix, refGenome);

    std::map<std::string, int> chrSizeMap;
    std::map<std::string, std::vector<Variant> > randomSdiMaps;
    for( std::map<std::string, Fasta>::iterator it=refGenome.begin(); it!=refGenome.end(); ++it ){
        chrSizeMap[it->first]=it->second.getSequence().size();
        randomSdiMaps[it->first]=std::vector<Variant>();
    }
    std::cout << "begin to random" << std::endl;
    for( std::map<std::string, std::vector<Variant> >::iterator it=sdiMaps.begin(); it!=sdiMaps.end(); ++it) {
        std::string chr = it->first;
        int chrSize = chrSizeMap[chr];
        std::cout << "begin to random " << chr << std::endl;

        reDoThisChrLable:
        int deletionNumber = 0;

        for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){ // the observed variant records
            if( it2->getChanginglength() < -10000 ){
                int results = ranDomItDeletion( (*it2), randomSdiMaps, chr, chrSize, refGenome);
                if( 1==results ){
                    randomSdiMaps[chr].clear();
                    std::cout << "begin to re-random " << chr << std::endl;
                    goto reDoThisChrLable;
                }
                ++deletionNumber;
                if( deletionNumber % 5000 ==0 ) {
                    std::cout << deletionNumber << " deletion finished" << std::endl;
                }
            }
        }

        for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){ // the observed variant records
            if( it2->getChanginglength() < -1000 &&it2->getChanginglength() >= -10000  ){
                int results = ranDomItDeletion( (*it2), randomSdiMaps, chr, chrSize, refGenome);
                if( 1==results ){
                    randomSdiMaps[chr].clear();
                    std::cout << "begin to re-random " << chr << std::endl;
                    goto reDoThisChrLable;
                }
                ++deletionNumber;
                if( deletionNumber % 5000 ==0 ) {
                    std::cout << deletionNumber << " deletion finished" << std::endl;
                }
            }
        }
        for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){ // the observed variant records
            if( it2->getChanginglength() < 0 && it2->getChanginglength() >= -1000 ){
                int results = ranDomItDeletion( (*it2), randomSdiMaps, chr, chrSize, refGenome );
                if( 1==results ){
                    randomSdiMaps[chr].clear();
                    std::cout << "begin to re-random " << chr << std::endl;
                    goto reDoThisChrLable;
                }
                ++deletionNumber;
                if( deletionNumber % 5000 ==0 ) {
                    std::cout << deletionNumber << " deletion finished" << std::endl;
                }
            }
        }

        int snpNumber = 0;
        for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){ // the observed variant records
            if( it2->getChanginglength() == 0 ){
                int results = ranDomItSnp( (*it2), randomSdiMaps, chr, chrSize, refGenome );
                if( 1==results ){
                    randomSdiMaps[chr].clear();
                    std::cout << "begin to re-random " << chr << std::endl;
                    goto reDoThisChrLable;
                }
                ++snpNumber;
                if( snpNumber % 5000 ==0 ) {
                    std::cout << snpNumber << " SNP finished" << std::endl;
                }
            }
        }

        int insertNumber = 0;
        for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){ // the observed variant records
            if( it2->getChanginglength() > 0 ) {
                int results = ranDomItInsertion((*it2), randomSdiMaps, chr, chrSize);
                if (1 == results) {
                    randomSdiMaps[chr].clear();
                    std::cout << "begin to re-random " << chr << std::endl;
                    goto reDoThisChrLable;
                }
                ++insertNumber;
                if (insertNumber % 5000 == 0) {
                    std::cout << insertNumber << " insert finished" << std::endl;
                }
            }
        }
    }

    std::ofstream ofileRandomSdi;
    ofileRandomSdi.open(outputPrefix);
    for(std::map<std::string, std::vector<Variant> >::iterator it=randomSdiMaps.begin(); it!=randomSdiMaps.end(); ++it){
        std::sort(it->second.begin(), it->second.end(), compare_sdi_record());
        for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){
            printSdiln(ofileRandomSdi, (*it2));
        }
    }
    ofileRandomSdi.close();
}

//
//**
// *
// * AA TT
// * AC GT
// * AT AT
// * AG CT
// *
// * CA TG
// * CC GG
// * CG CG
// * CT AG
// *
// * GA TC
// * GC GC
// * GG CC
// * GT AC
// *
// * TA TA
// * TC GA
// * TG CA
// * TT AA
// *
// * THEN WHAT IS LEFT
// *
// * AA/TT
// * AC/GT
// * AT/AT
// * AG/CT
// * CA/TG
// * CC/GG
// * CG/CG
// * GA/TC
// * GC/GC
// * TA/TA
// * /

/*
double snp combinations:
AA - CC  CC - AA  TT - GG  GG - TT
AA - CG  CG - AA  TT - CG  CG - TT
AA - CT  CT - AA  TT - AG  AG - TT
AA - GC  GC - AA  TT - GC  GC - TT
AA - GG  GG - AA  TT - CC  CC - TT
AA - GT  GT - AA  TT - AC  AC - TT
AA - TC  TC - AA  TT - GA  GA - TT
AA - TG  TG - AA  TT - CA  CA - TT
AA - TT  TT - AA

AC - CA  CA - AC  TG - GT  GT - TG
AC - CG  CG - AC  CG - GT  GT - CG
AC - CT  CT - AC  AG - GT  GT - AG
AC - GA  GA - AC  TC - GT  GT - TC
AC - GG  GG - AC  CC - GT  GT - CC
AC - GT  GT - AC
AC - TA  TA - AC  TA - GT  GT - TA
AC - TG  TG - AC  CA - GT  GT - CA

AG - CA  CA - AG  CT - TG  TG - CT
AG - CC  CC - AG  CT - GG  GG - CT
AG - CT  CT - AG
AG - GA  GA - AG  CT - TC  TC - CT
AG - GC  GC - AG  CT - GC  GC - CT
AG - TA  TA - AG  CT - TA  TA - CT
AG - TC  TC - AG  CT - GA  GA - CT

AT - CA  CA - AT  AT - TG  TG - AT
AT - CC  CC - AT  AT - GG  GG - AT
AT - CG  CG - AT
AT - GA  GA - AT  AT - TC  TC - AT
AT - GC  GC - AT
AT - TA  TA - AT

CA - GC  GC - CA  TG - GC  GC - TG
CA - GG  GG - CA  TG - CC  CC - TG
CA - TC  TC - CA  TG - GA  GA - TG
CA - TG  TG - CA

CC - GA  GA - CC  GG - TC  TC - GG
CC - GG  GG - CC
CC - TA  TA - CC  GG - TA  TA - GG

CG - GA  GA - CG  CG - TC  TC - CG
CG - GC  GC - CG
CG - TA  TA - CG

GA - TC  TC - GA

GC - TA  TA - GC


 */