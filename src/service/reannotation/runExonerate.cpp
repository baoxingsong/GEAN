//
// Created by baoxing on 10/10/17.
//

#include "runExonerate.h"

std::mutex gmutexRunExonerateEst;


void runExonerateEst(std::string& transcriptName, std::string& cdsSequence, std::string& targetSequence,
                     NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                     std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                     STRAND strand, std::string& tchromeSomeName, std::string& prefixUuid , std::map<std::string, Fasta>& targetGenome,
                     std::map<std::string, std::string>& parameters, int& minIntron){
    try{
        std::string tempFolder = get_parameters("tempFolder", parameters);
        gmutexRunExonerateEst.lock();
        std::string tempFile = tempFolder + "/" + generateUUID( prefixUuid );
        gmutexRunExonerateEst.unlock();
        std::string cdsFile = tempFile + "cds.fasta";
        std::string targetFile = tempFile + "target.fasta";

        std::ofstream ofile2;
        ofile2.open(cdsFile);
        ofile2 << ">cds" << std::endl << cdsSequence << std::endl;
        ofile2.close();

        std::ofstream ofile;
        ofile.open(targetFile);
        ofile << ">target" << std::endl << targetSequence << std::endl;
        ofile.close();
        std::string maxintron = get_parameters("maxintronlength", parameters);
        std::string minintron = get_parameters("minintronlength", parameters);
        std::string command =
                "exonerate --maxintron " + maxintron + " --model est2genome -i -10 --score 10 --bestn 1 --minintron " + minintron + " " +
                cdsFile + " " + targetFile + " --showtargetgff true >" + tempFile;
        system(&command[0]);
        readExonerateEstResult(transcriptName, cdsSequence, targetSequence,
                               nucleotideCodeSubstitutionMatrix,
                               targetTranscriptsHashMap, startTarget, endTarget,
                               strand, tchromeSomeName, tempFile, targetGenome, parameters, minIntron);
        std::string cleanFileCommand = "rm " + cdsFile;
        system(&cleanFileCommand[0]);
        cleanFileCommand = "rm " + targetFile;
        system(&cleanFileCommand[0]);
        cleanFileCommand = "rm " + tempFile;
        system(&cleanFileCommand[0]);
    }catch (...){
        return;
    }
}

void readExonerateEstResult( std::string& transcriptName, std::string& cdsSequence, std::string& targetSequence,
                             NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                             std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                             STRAND& strand, std::string& tchromeSomeName, std::string& fileLocation, std::map<std::string, Fasta>& targetGenome,
                             std::map<std::string, std::string>& parameters, int& minIntron){

    if( (!if_file_exists(fileLocation)) || GetFileSize(fileLocation)<1 ){
        return;
    }

    Transcript targetTranscript(transcriptName, tchromeSomeName, strand );
    int cdsNumber = 0;
    if( POSITIVE ==strand ){
        std::ifstream infile2(fileLocation);
        std::regex reg2("^(target)\\t(\\S*)\\texon\\t(\\S*)\\t(\\S*)\\t(\\S*)\\t\\+");
        std::string line2;
        while (std::getline(infile2, line2)){
            std::smatch match;
            regex_search(line2, match, reg2);

            if( match.empty() ){
            }else{
                int start = stoi(match[3]);
                int end = stoi(match[4]);
                if(start>end){
                    int temp=start;
                    start = end;
                    end = temp;
                }
                start = start + startTarget -1;
                end = end + startTarget - 1;
                GenomeBasicFeature cds(start, end);
                targetTranscript.addCds(cds);
                ++cdsNumber;
            }
        }
        infile2.close();
    }else{
        std::ifstream infile2(fileLocation);
        std::regex reg2("^(target)\\t(\\S*)\\texon\\t(\\S*)\\t(\\S*)\\t(\\S*)\\t\\-");
        std::string line2;
        while (std::getline(infile2, line2)){
            std::smatch match;
            regex_search(line2, match, reg2);

            if( match.empty() ){
            }else{
                int start = stoi(match[3]);
                int end = stoi(match[4]);
                if(start>end){
                    int temp=start;
                    start = end;
                    end = temp;
                }
                start = start + startTarget -1;
                end = end + startTarget - 1;
                GenomeBasicFeature cds(start, end);
                targetTranscript.addCds(cds);
                ++cdsNumber;
            }
        }
        infile2.close();
    }
    if( cdsNumber > 0 && !targetTranscript.getCdsVector().empty() ){
//        std::cout << transcriptName << " exonerate est updateInfor begin" << std::endl;
        TranscriptUpdateCdsInformation(targetTranscript, targetGenome);
        //targetTranscript.updateInfor(targetGenome);

//        std::cout << transcriptName << " exonerate est checkOrfState begin" << std::endl;
//        checkOrfState( targetTranscript, referenceTranscript,
//                       targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix);
        checkOrfState( targetTranscript, targetGenome, nucleotideCodeSubstitutionMatrix, minIntron);
//        std::cout << transcriptName << " exonerate est checkOrfState finished" << std::endl;

        if( targetTranscript.getIfOrfShift() ){

//            if( targetTranscript.getMetaInformation().find("spliceSitesDestroyed")!=std::string::npos ){
//                std::cout << "879 there is something wrong with the splice sites of cds based alignment " << targetTranscript.getMetaInformation() << std::endl;
//            }else{
//                std::cout << "881 " << targetTranscript.getMetaInformation() << std::endl;
//            }
            std::string protenSequene = nA2AA(cdsSequence, nucleotideCodeSubstitutionMatrix);
            runExonerateProtein(transcriptName, protenSequene, targetSequence,
                                nucleotideCodeSubstitutionMatrix,
                                targetTranscriptsHashMap, startTarget, endTarget,
                                strand, tchromeSomeName, fileLocation, targetGenome, parameters, minIntron);
        }else{
            targetTranscript.setSource("CDSALIGNMENT");
            gmutexRunExonerateEst.lock();
            targetTranscriptsHashMap[transcriptName]=targetTranscript;
            gmutexRunExonerateEst.unlock();
        }
    }else{
        std::string protenSequene = nA2AA(cdsSequence, nucleotideCodeSubstitutionMatrix);
        runExonerateProtein(transcriptName, protenSequene, targetSequence,
                            nucleotideCodeSubstitutionMatrix,
                            targetTranscriptsHashMap, startTarget, endTarget,
                            strand, tchromeSomeName, fileLocation, targetGenome, parameters, minIntron );
    }
}

void runExonerateProtein(std::string& transcriptName, std::string& protenSequene, std::string& targetSequence,
                         NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                         std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                         STRAND& strand, std::string& tchromeSomeName, std::string& fileLocation , std::map<std::string, Fasta>& targetGenome,
                         std::map<std::string, std::string>& parameters, int& minIntron ) {
    std::string proteinFile = fileLocation + "protein.fasta";
    std::string targetFile = fileLocation + "target.fasta";
    std::string tempFile=fileLocation + "protein.aln";

    std::ofstream ofile;
    ofile.open(proteinFile);
    ofile << ">protein" << std::endl << protenSequene << std::endl;
    ofile.close();

    std::string maxintron = get_parameters("maxintronlength", parameters);
    std::string minintron = get_parameters("minintronlength", parameters);

    std::string command = "exonerate --bestn 1 --maxintron " + maxintron + " --intronpenalty -10 --model protein2genome --percent 10 --score 10 --minintron " + minintron + " "+ proteinFile + " " + targetFile + " --showtargetgff true >" +tempFile;
    system(&command[0]);
    readExonerateProteinResult(tempFile, nucleotideCodeSubstitutionMatrix, startTarget, endTarget, strand, targetTranscriptsHashMap, transcriptName, tchromeSomeName, targetGenome, minIntron) ;
    std::string cleanFileCommand = "rm " + proteinFile;
    system(&cleanFileCommand[0]);
    cleanFileCommand = "rm " + tempFile;
    system(&cleanFileCommand[0]);
}

void readExonerateProteinResult( std::string& fileLocation, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                 int& startTarget, int& endTarget, STRAND& strand, std::map<std::string, Transcript>& targetTranscriptsHashMap,
                                 std::string& transcriptName, std::string& tchromeSomeName, std::map<std::string, Fasta>& targetGenome, int& minIntron){
    if( (!if_file_exists(fileLocation)) || GetFileSize(fileLocation)<1 ){
        return;
    }

    Transcript targetTranscript(transcriptName, tchromeSomeName, strand );
    targetTranscript.setSource("PROTEINALIGNMENT");
    int cdsNumber = 0;
    if( POSITIVE ==strand ){
        std::ifstream infile(fileLocation);
        std::regex reg("^(target)\t([\\s\\S]*)\tcds\t(\\S*)\t(\\S*)\t(\\S*)\t\\+");
        std::string line;
        while (std::getline(infile, line)){
            std::smatch match;
            regex_search(line, match, reg);

            if( match.empty() ){
            }else{
                int start = stoi(match[3]);
                int end = stoi(match[4]);
                if(start>end){
                    int temp=start;
                    start = end;
                    end = temp;
                }
                start = start + startTarget -1;
                end = end + startTarget - 1;
                GenomeBasicFeature cds(start, end);
                targetTranscript.addCds(cds);
                ++cdsNumber;
                //std::cout << start << " " << end << std::endl;
            }
        }
        infile.close();
    }else {
        std::ifstream infile(fileLocation);
        std::regex reg("^(target)\t([\\s\\S]*)\tcds\t(\\S*)\t(\\S*)\t(\\S*)\t\\-");
        std::string line;
        while (std::getline(infile, line)){
            std::smatch match;
            regex_search(line, match, reg);

            if( match.empty() ){
            }else{
                int start = stoi(match[3]);
                int end = stoi(match[4]);
                if(start>end){
                    int temp=start;
                    start = end;
                    end = temp;
                }
                start = start + startTarget -1;
                end = end + startTarget - 1;
                GenomeBasicFeature cds(start, end);
                targetTranscript.addCds(cds);
                ++cdsNumber;
            }
        }
        infile.close();
    }
    if( cdsNumber > 0 && !targetTranscript.getCdsVector().empty() ) {
        TranscriptUpdateCdsInformation(targetTranscript, targetGenome);
        checkOrfState(targetTranscript, targetGenome, nucleotideCodeSubstitutionMatrix, minIntron);
        if (!targetTranscript.getIfOrfShift()) {
            gmutexRunExonerateEst.lock();
            targetTranscriptsHashMap[transcriptName] = targetTranscript;
            gmutexRunExonerateEst.unlock();
        }
//        else{
//            if( targetTranscript.getMetaInformation().find("spliceSitesDestroyed")!=std::string::npos ){
//                std::cout << "991 there is something wrong with the splice sites of cds based alignment " << targetTranscript.getMetaInformation() << std::endl;
//            }else{
//                std::cout << "993 " << targetTranscript.getMetaInformation() << std::endl;
//            }
//        }
    }
}
