//
// Created by baoxing on 10/10/17.
//

#include "cutWindow.h"
std::mutex g_num_mutex;


void writeToFileParallelVersion_ref(int const index, int const & binSize, std::vector<Window>& windows, int const & thisChrLength, int const & msaWindowOverlap, std::string const & chr,
                                    std::map<std::string, Fasta> & sequences, std::string const & accessionName, std::string const & folder,
                                    std::atomic_int & number_of_runing_threads){
    //std::cout << "there are " << number_of_runing_threads << " threads running for reference output" << std::endl;
//    try{
        int checked = 0;
        for( int i=index; checked<binSize && i<windows.size(); ++i, ++checked ){
            int start;
            if ((windows[i]).getPreExtand()) {
                start = (windows[i]).getStart() - msaWindowOverlap;
            } else {
                start = (windows[i]).getStart();
            }
            if (start < 1) {
                start = 1;
            }
            int end;
            if ((windows[i]).getPostExtand()) {
                end = (windows[i]).getEnd() + msaWindowOverlap;
            } else {
                end = (windows[i]).getEnd();
            }
            if( end > thisChrLength ){
                end = thisChrLength;
            }

            std::ofstream ofile2;
            ofile2.open( folder+"/"+std::to_string( windows[i].getStart() ) + "_" + std::to_string( windows[i].getEnd()),  std::ofstream::app | std::ofstream::out );
            ofile2 << ">" << accessionName << "_" << std::to_string(start) << "_" + std::to_string(end) << std::endl << getSubsequence(sequences, chr, start, end) << std::endl;
            ofile2.close();
        }
//    }catch (...) {
//        std::cerr << " there is something wrong with the reference sequence output" << std::endl;
//        --number_of_runing_threads;
//        return;
//    }
    --number_of_runing_threads;
}


void outputReferenceSequence_one_chromosome( std::map<std::string, std::vector<Window>>& chromosomeToWindowMap, std::string & tempFolder, int & maxThread,
    int & msaWindowOverlap, std::map<std::string, Fasta> referenceGenome,
    std::string thisChr, std::atomic_int& number_of_runing_threads_chromosome, std::atomic_int& number_of_runing_threads_window){

    std::string accessionName = "ref";
    std::string outputFolder = tempFolder + "/" + thisChr;
    int thisChrLength = referenceGenome[thisChr].getSequence().length();

    int binSize = 500;

//    for (std::vector<Window>::iterator itwin = chromosomeToWindowMap[thisChr].begin();
//         itwin != chromosomeToWindowMap[thisChr].end(); ++itwin) {
    for( int index=0; index < chromosomeToWindowMap[thisChr].size(); index+=binSize ){
        bool isThisThreadUnrun = true;
        while (isThisThreadUnrun) {
            if ( (number_of_runing_threads_window+number_of_runing_threads_chromosome) < maxThread) {
                std::thread t(writeToFileParallelVersion_ref, index,  std::ref(binSize), std::ref(chromosomeToWindowMap[thisChr]), std::ref(thisChrLength), std::ref(msaWindowOverlap), std::ref(thisChr), std::ref(referenceGenome),
                              std::ref(accessionName), std::ref(outputFolder), std::ref(number_of_runing_threads_window));
                ++number_of_runing_threads_window;
                t.detach();
                isThisThreadUnrun = false;
                break;
            }else {
                usleep(10);
            }
        }
    }
    while (number_of_runing_threads_window > 0) {// wait for all the thread
        usleep(100);
    }
    --number_of_runing_threads_chromosome;
}

void outputReferenceSequence(std::map<std::string, std::vector<Window>>& chromosomeToWindowMap, std::string & tempFolder, int & maxThread, int & msaWindowOverlap, std::map<std::string, Fasta> referenceGenome){
    // output reference file begin
    std::atomic_int number_of_runing_threads_chromosome(0);
    std::atomic_int number_of_runing_threads_window(0);

    for( std::map<std::string, std::vector<Window>>::iterator ittm=chromosomeToWindowMap.begin();
         ittm!=chromosomeToWindowMap.end(); ++ittm ) {
        std::string thisChr = ittm->first;

        bool isThisThreadUnrun = true;
        while (isThisThreadUnrun) {
            if ((number_of_runing_threads_chromosome) < maxThread) {
                std::thread t(outputReferenceSequence_one_chromosome, std::ref(chromosomeToWindowMap), std::ref(tempFolder), std::ref(maxThread), std::ref(msaWindowOverlap), std::ref(referenceGenome),
                              thisChr, std::ref(number_of_runing_threads_chromosome), std::ref(number_of_runing_threads_window));
                ++number_of_runing_threads_chromosome;
                t.detach();
                isThisThreadUnrun = false;
                break;
            } else {
                usleep(10);
            }
        }
    }// output reference file end
    while (number_of_runing_threads_chromosome > 0) {// wait for all the thread
        usleep(100);
    }
}


// this function is used to output the sequence for non-reference accessions
void writeToFileParallelVersion_non_ref_intergenetic( int wini, std::string& thisChr, std::map<std::string, Fasta>& sequences,
     std::string& accessionName, std::string & folder, std::map<std::string, std::vector<Window>>& chromosomeToWindowMap, int& msaWindowOverlap,
     std::map<std::string, std::vector<Variant>> & variantsMaps, std::map<std::string, Transcript>& targetTranscriptsHashMap,
     std::atomic_int & number_of_runing_threads, FileWriteCatch & fileWriteCatch, int lastOrfConservedTranscriptEndPlus_1, int nextOrfConservedTranscriptStartMinut_1){
    try{
        int start;
        int end;
        if( wini == 0 ){
            start = 1;
        }else if( chromosomeToWindowMap[thisChr][wini-1].get_transcript_name().length() > 1 &&
                  !targetTranscriptsHashMap[chromosomeToWindowMap[thisChr][wini-1].get_transcript_name()].getIfOrfShift()){
            start = targetTranscriptsHashMap[chromosomeToWindowMap[thisChr][wini-1].get_transcript_name()].getPEnd() + 1;
            //start = intergeneticMinStart;
            //std::cout << "1502" << std::endl;
        }else if ( chromosomeToWindowMap[thisChr][wini-1].get_transcript_name().length() > 1 &&
                   targetTranscriptsHashMap[ chromosomeToWindowMap[thisChr][wini-1].get_transcript_name()].getIfOrfShift() ){
            start = getChangedFromBasement(thisChr,  chromosomeToWindowMap[thisChr][wini].getStart(), variantsMaps) - msaWindowOverlap;
            //std::cout << "1506" << std::endl;
        }else{
            start = getChangedFromBasement(thisChr, chromosomeToWindowMap[thisChr][wini].getStart(), variantsMaps) - msaWindowOverlap;
            //std::cout << "1509" << std::endl;
        }

        if( wini == ( chromosomeToWindowMap[thisChr].size()-1) ){ // this is the last window
            end = sequences[thisChr].getSequence().length();
            //std::cout << "1519" << std::endl;
        }else if( chromosomeToWindowMap[thisChr][wini+1].get_transcript_name().length() > 1 &&
                  !targetTranscriptsHashMap[ chromosomeToWindowMap[thisChr][wini+1].get_transcript_name()].getIfOrfShift()){
            end = targetTranscriptsHashMap[ chromosomeToWindowMap[thisChr][wini+1].get_transcript_name()].getPStart() - 1;

            //std::cout << "1523" << std::endl;
        }else if (  chromosomeToWindowMap[thisChr][wini+1].get_transcript_name().length() > 1 &&
                    targetTranscriptsHashMap[ chromosomeToWindowMap[thisChr][wini+1].get_transcript_name()].getIfOrfShift() ){
            end = getChangedFromBasement(thisChr, chromosomeToWindowMap[thisChr][wini].getEnd(), variantsMaps) + msaWindowOverlap;
            //std::cout << "1527" << std::endl;
        }else{
            end = getChangedFromBasement(thisChr, chromosomeToWindowMap[thisChr][wini].getEnd(), variantsMaps) + msaWindowOverlap;
            //std::cout << "1529" << std::endl;
        }

        if( start < lastOrfConservedTranscriptEndPlus_1){
            start = lastOrfConservedTranscriptEndPlus_1;
        }
        if( end > nextOrfConservedTranscriptStartMinut_1 ){
            end = nextOrfConservedTranscriptStartMinut_1;
        }

        int chrLength = sequences[thisChr].getSequence().length();
        if( end > chrLength ){
            end = chrLength;
        }
        if( start > end ){
            --start;
            end = start ;
            if( end > chrLength ){
                end = chrLength;
                start = end;
            }
        }
        std::string fileName = folder+"/"+std::to_string(chromosomeToWindowMap[thisChr][wini].getStart()) + "_" + std::to_string(chromosomeToWindowMap[thisChr][wini].getEnd());
        std::stringstream ofile2;
        ofile2 << ">" << accessionName << "_" << std::to_string(start) << "_" << std::to_string(end) << std::endl << getSubsequence(sequences, thisChr, start, end) << std::endl;
        std::string content = ofile2.str();
        g_num_mutex.lock();
        fileWriteCatch.addContent(fileName, content);
        g_num_mutex.unlock();

    }catch (...) {
        std::cerr << " there is something wrong with the " << accessionName << "intergenetic sequence output" << std::endl;
        --number_of_runing_threads;
        return;
    }
    --number_of_runing_threads;
}

void writeToFileParallelVersion_non_ref_genetic_orf_shift(int wini, std::map<std::string, std::vector<Window>>& chromosomeToWindowMap, std::string& thisChr, std::map<std::string, Fasta>& sequences,
     std::string& accessionName, std::string & folder, int& msaWindowOverlap, std::map<std::string, std::vector<Variant>> & variantsMaps,
     std::atomic_int & number_of_runing_threads, FileWriteCatch & fileWriteCatch, int lastOrfConservedTranscriptEndPlus_1, int nextOrfConservedTranscriptStartMinut_1){
    try{
        int windowsSize =  chromosomeToWindowMap[thisChr][wini].getEnd() -  chromosomeToWindowMap[thisChr][wini].getStart() + 1;
        int thisOverLapSize = ceil(windowsSize *1.5);
        if( thisOverLapSize > msaWindowOverlap ){
            thisOverLapSize = msaWindowOverlap;
        }
        int start = getChangedFromBasement(thisChr,  chromosomeToWindowMap[thisChr][wini].getStart(), variantsMaps) - thisOverLapSize;
        int end = getChangedFromBasement(thisChr,  chromosomeToWindowMap[thisChr][wini].getEnd(), variantsMaps) + thisOverLapSize;

        if( start < lastOrfConservedTranscriptEndPlus_1){
            start = lastOrfConservedTranscriptEndPlus_1;
        }
        if( end > nextOrfConservedTranscriptStartMinut_1 ){
            end = nextOrfConservedTranscriptStartMinut_1;
        }

        int chrLength = sequences[thisChr].getSequence().length();
        if( end > chrLength ){
            end = chrLength;
        }
        if( start > end ){
            --start;
            end = start ;
            if( end > chrLength ){
                end = chrLength;
                start = end;
            }
        }
//        std::string file_name = ;
  //      std::string content = ">" + accessionName + "_" + std::to_string(start) + "_" + std::to_string(end) + "\n" + getSubsequence(sequences, thisChr, start, end)+"\n";

        std::string fileName = folder+"/"+std::to_string( chromosomeToWindowMap[thisChr][wini].getStart()) + "_" + std::to_string( chromosomeToWindowMap[thisChr][wini].getEnd());
        std::stringstream ofile2;
        ofile2 << ">" << accessionName << "_" << std::to_string(start) << "_" << std::to_string(end) << std::endl << getSubsequence(sequences, thisChr, start, end) << std::endl;
        std::string content = ofile2.str();
        g_num_mutex.lock();
        fileWriteCatch.addContent(fileName, content);
        g_num_mutex.unlock();

    }catch (...) {
        std::cerr << " there is something wrong with the " << accessionName << "intergenetic sequence output" << std::endl;
        --number_of_runing_threads;
        return;
    }
    --number_of_runing_threads;
}

void writeToFileParallelVersion_non_ref_genetic_orf_conserved(int wini, std::map<std::string, std::vector<Window>>& chromosomeToWindowMap, std::string& thisChr, std::map<std::string, Fasta>& sequences,
    std::string& accessionName, std::string & folder, int& msaWindowOverlap, std::map<std::string, std::vector<Variant>> & variantsMaps,
    std::atomic_int & number_of_runing_threads, std::map<std::string, Transcript>& targetTranscriptsHashMap, std::map<std::string, std::vector<Window>>& transcriptToWindowMap,
    FileWriteCatch & fileWriteCatch, int cdsId, int winIndexInTranscript, std::string transcriptName){ // the last three parameters should not be reference
    try{

        int start;
        int end;

        if( winIndexInTranscript == 0){
//                                std::cout << "1556" << std::endl;
            start = (targetTranscriptsHashMap[transcriptName]).getPStart();
        } else {
            if( CDSINTRON == transcriptToWindowMap[transcriptName][winIndexInTranscript].get_region() &&
                CDS == transcriptToWindowMap[transcriptName][winIndexInTranscript-1].get_region() ){
                //                                  std::cout << "1561" << std::endl;
                start = (targetTranscriptsHashMap[transcriptName].getCdsVector())[cdsId].getEnd() + 1;
            }else if( CDS == transcriptToWindowMap[transcriptName][winIndexInTranscript].get_region() &&
                      CDSINTRON == transcriptToWindowMap[transcriptName][winIndexInTranscript-1].get_region()  ){
                //                                std::cout << "1565" << std::endl;
                start = (targetTranscriptsHashMap[transcriptName].getCdsVector())[cdsId].getStart();
            }else{
                //                              std::cout << "1568" << std::endl;
                start = getChangedFromBasement(thisChr, transcriptToWindowMap[transcriptName][winIndexInTranscript].getStart(), variantsMaps) - msaWindowOverlap;
                if( CDS == transcriptToWindowMap[transcriptName][winIndexInTranscript].get_region()  ){
                    int minStart = (targetTranscriptsHashMap[transcriptName].getCdsVector())[cdsId].getStart();
                    if( start < minStart ){
                        //                                    std::cout << "1573" << std::endl;
                        start = minStart;
                    }
                }else if( CDSINTRON == transcriptToWindowMap[transcriptName][winIndexInTranscript].get_region()  ){
                    int minStart = (targetTranscriptsHashMap[transcriptName].getCdsVector())[cdsId].getEnd() + 1;
                    if( start < minStart ){
                        //                                    std::cout << "1581" << std::endl;
                        start = minStart;
                    }
                }
            }
        }

        if( winIndexInTranscript == (transcriptToWindowMap[transcriptName].size()-1) ){ // last window of this
//                                std::cout << "1589" << std::endl;
            end = targetTranscriptsHashMap[transcriptName].getPEnd();
        } else {
            if ( CDS == transcriptToWindowMap[transcriptName][winIndexInTranscript].get_region() &&
                 CDSINTRON == transcriptToWindowMap[transcriptName][winIndexInTranscript+1].get_region() ){
//                                    std::cout << "1594" << std::endl;
                end = (targetTranscriptsHashMap[transcriptName].getCdsVector())[cdsId].getEnd();
            } else if ( CDSINTRON == transcriptToWindowMap[transcriptName][winIndexInTranscript].get_region() &&
                        CDS == transcriptToWindowMap[transcriptName][winIndexInTranscript+1].get_region()  ){
                //                                  std::cout << "1598" << std::endl;
                end = (targetTranscriptsHashMap[transcriptName].getCdsVector())[cdsId+1].getStart() -1;
            } else {
//                                    std::cout << "1601" << std::endl;
                end = getChangedFromBasement(thisChr, transcriptToWindowMap[transcriptName][winIndexInTranscript].getEnd(), variantsMaps) + msaWindowOverlap;
                if( CDS == transcriptToWindowMap[transcriptName][winIndexInTranscript].get_region()  ){
                    int maxEnd = (targetTranscriptsHashMap[transcriptName].getCdsVector())[cdsId].getEnd();
                    if( end > maxEnd ){
//                                            std::cout << "1606" << std::endl;
                        end = maxEnd;
                    }
                }else if( CDSINTRON == transcriptToWindowMap[transcriptName][winIndexInTranscript].get_region()  ){
                    int maxEnd = (targetTranscriptsHashMap[transcriptName].getCdsVector())[cdsId+1].getStart() -1;
                    if( end > maxEnd ){
//                                            std::cout << "1613" << std::endl;
                        end = maxEnd;
                    }
                }
            }
        }

        int chrLength = sequences[thisChr].getSequence().length();
        if( end > chrLength ){
            end = chrLength;
        }
        if( start > end ){
            --start;
            end = start ;
            if( end > chrLength ){
                end = chrLength;
                start = end;
            }
        }
//        std::string file_name = ;
//        std::string content = ">" + accessionName + "_" + std::to_string(start) + "_" + std::to_string(end) + "\n" + getSubsequence(sequences, thisChr, start, end)+"\n";

        std::string fileName = folder+"/"+std::to_string( chromosomeToWindowMap[thisChr][wini].getStart()) + "_" + std::to_string( chromosomeToWindowMap[thisChr][wini].getEnd());
        std::stringstream ofile2;
        ofile2 << ">" << accessionName << "_" << std::to_string(start) << "_" << std::to_string(end) << std::endl << getSubsequence(sequences, thisChr, start, end) << std::endl;
        std::string content = ofile2.str();
        g_num_mutex.lock();
        fileWriteCatch.addContent(fileName, content);
        g_num_mutex.unlock();

    }catch (...) {
        std::cerr << " there is something wrong with the " << accessionName << "inter-genetic sequence output" << std::endl;
        --number_of_runing_threads;
        return;
    }
    --number_of_runing_threads;
}

void writeToFileParallelVersion_certain_chromosome_non_ref(std::map<std::string, std::vector<Window>> & chromosomeToWindowMap, std::string& accessionName,
      std::string &tempFolder, int &maxThread, std::map<std::string, Transcript>& targetTranscriptsHashMap, std::map<std::string, Fasta>& targetGenome,
      std::string thisChr, std::map<std::string, std::vector<Variant>> & variantsMaps, int & msaWindowOverlap,
      std::map<std::string, std::vector<Window>>& transcriptToWindowMap, std::atomic_int & number_of_runing_threads_chromosome,
      std::atomic_int& number_of_runing_threads_window, FileWriteCatch& fileWriteCatch){

    std::string outputFolder = tempFolder + "/" + thisChr;
    int lastOrfConservedTranscriptEndPlus_1 = 1;

    int nextOrfConservedTranscriptStartMinut_1 = targetGenome[thisChr].getSequence().size();
    int nextOrfConservedTranscriptStartMinut_1_index = 0;

    for ( int wini=0; wini < chromosomeToWindowMap[thisChr].size(); ){
        if( nextOrfConservedTranscriptStartMinut_1_index <= wini ){
            for( int tempWini = wini; tempWini< chromosomeToWindowMap[thisChr].size(); ++tempWini ){
                if( (chromosomeToWindowMap[thisChr][tempWini]).get_transcript_name().length() > 1
                    && (!targetTranscriptsHashMap[(chromosomeToWindowMap[thisChr][tempWini]).get_transcript_name()].getIfOrfShift()) ){
                    nextOrfConservedTranscriptStartMinut_1 = targetTranscriptsHashMap[(chromosomeToWindowMap[thisChr][tempWini]).get_transcript_name()].getPStart() -1;
                    nextOrfConservedTranscriptStartMinut_1_index = tempWini;
                    tempWini = chromosomeToWindowMap[thisChr].size(); // stop the loop
                    break;
                }
                nextOrfConservedTranscriptStartMinut_1 = targetGenome[thisChr].getSequence().size();
                // if could not go out of the loop, then there is no ORF conserved transcript left, and take the chromosome length as maximum end value
            }
        }

        if( (chromosomeToWindowMap[thisChr][wini]).get_transcript_name().length() < 2 ){ // windows not in transcript region
            bool isThisThreadUnrun = true;
            while (isThisThreadUnrun) {
                if ( (number_of_runing_threads_window+number_of_runing_threads_chromosome) < maxThread) {
                    std::thread t(writeToFileParallelVersion_non_ref_intergenetic, wini, std::ref(thisChr), std::ref(targetGenome),
                    std::ref(accessionName), std::ref(outputFolder), std::ref(chromosomeToWindowMap), std::ref(msaWindowOverlap),
                    std::ref(variantsMaps), std::ref(targetTranscriptsHashMap), std::ref(number_of_runing_threads_window), std::ref(fileWriteCatch),
                    lastOrfConservedTranscriptEndPlus_1, nextOrfConservedTranscriptStartMinut_1 ); // the last two parameter should not be reference
                    ++number_of_runing_threads_window;
                    t.detach();
                    isThisThreadUnrun = false;
                    break;
                }else {
                    usleep(10);
                }
            }
            ++wini;
        } else {
            std::string transcriptName = (chromosomeToWindowMap[thisChr][wini]).get_transcript_name();
            if( targetTranscriptsHashMap[transcriptName].getIfOrfShift() ){
                for ( int i=0; i< transcriptToWindowMap[transcriptName].size(); ++i ) {
                    //Window refWindow = transcriptToWindowMap[transcriptName][i];
                    bool isThisThreadUnrun = true;
                    while (isThisThreadUnrun) {
                        if ((number_of_runing_threads_window+number_of_runing_threads_chromosome) < maxThread) {
                            std::thread t(writeToFileParallelVersion_non_ref_genetic_orf_shift, wini,  std::ref(chromosomeToWindowMap), std::ref(thisChr), std::ref(targetGenome),
                                std::ref(accessionName), std::ref(outputFolder), std::ref(msaWindowOverlap), std::ref(variantsMaps),
                                std::ref(number_of_runing_threads_window), std::ref(fileWriteCatch), lastOrfConservedTranscriptEndPlus_1, nextOrfConservedTranscriptStartMinut_1);

                            ++number_of_runing_threads_window;
                            t.detach();
                            isThisThreadUnrun = false;
                            break;
                        }else {
                            usleep(10);
                        }
                    }
                    ++wini;
                }
            } else {
                lastOrfConservedTranscriptEndPlus_1 = targetTranscriptsHashMap[transcriptName].getPEnd(); // updata this value for intergenetic region and ORF-shift region windows
                int cdsId = -1;
                for( int i=0; i < transcriptToWindowMap[transcriptName].size(); ++i ){
                    if( 0 == i  ){
                        ++cdsId;
                    } else if ( CDS == transcriptToWindowMap[transcriptName][i].get_region()
                                && CDSINTRON == transcriptToWindowMap[transcriptName][i-1].get_region()){
                        ++cdsId;
                    }
                    bool isThisThreadUnrun = true;
                    while (isThisThreadUnrun) {
                        if ((number_of_runing_threads_window+number_of_runing_threads_chromosome) < maxThread) {

                            std::thread t(writeToFileParallelVersion_non_ref_genetic_orf_conserved, wini, std::ref(chromosomeToWindowMap), std::ref(thisChr),
                                 std::ref(targetGenome), std::ref(accessionName), std::ref(outputFolder), std::ref(msaWindowOverlap), std::ref(variantsMaps),
                                 std::ref(number_of_runing_threads_window), std::ref(targetTranscriptsHashMap),std::ref(transcriptToWindowMap), std::ref(fileWriteCatch),
                                 cdsId, i, transcriptName);

                            ++number_of_runing_threads_window;
                            t.detach();
                            isThisThreadUnrun = false;
                            break;
                        } else {
                            usleep(10);
                        }
                    }
                    ++wini;
                }
            }
        }
    }
    while (number_of_runing_threads_window > 0) {// wait for all the thread
        usleep(100);
    }
    std::cout << accessionName << " " << thisChr << "output done" << std::endl;
    --number_of_runing_threads_chromosome;
}


void outputNonreferenceSequence(std::map<std::string, std::vector<Window>> & chromosomeToWindowMap, std::string& accessionName,
    std::string& tempFolder, int& maxThread, std::map<std::string, Transcript>& targetTranscriptsHashMap, std::map<std::string, Fasta>& targetGenome,
    std::map<std::string, std::vector<Variant>>& variantsMaps, int& msaWindowOverlap, std::map<std::string, std::vector<Window>>& transcriptToWindowMap, FileWriteCatch & fileWriteCatch ){

    std::atomic_int number_of_runing_threads_chromosome(0);
    std::atomic_int number_of_runing_threads_window(0);

    for ( std::map<std::string, std::vector<Window>>::iterator chrIt=chromosomeToWindowMap.begin();
          chrIt!=chromosomeToWindowMap.end(); ++chrIt){
        std::string thisChr = chrIt->first;// this is a local variable, should not be used as reference for multiple threads function

        bool isThisThreadUnrun = true;
        while (isThisThreadUnrun) {
            if ((number_of_runing_threads_window+number_of_runing_threads_chromosome) < maxThread) {
                std::thread t(writeToFileParallelVersion_certain_chromosome_non_ref, std::ref(chromosomeToWindowMap), std::ref(accessionName),
                    std::ref(tempFolder), std::ref(maxThread), std::ref(targetTranscriptsHashMap), std::ref(targetGenome),
                    thisChr, std::ref(variantsMaps), std::ref(msaWindowOverlap),
                    std::ref(transcriptToWindowMap), std::ref(number_of_runing_threads_chromosome), std::ref(number_of_runing_threads_window), std::ref(fileWriteCatch) );

                ++number_of_runing_threads_chromosome;
                t.detach();
                isThisThreadUnrun = false;
                break;
            }else {
                usleep(10);
            }
        }
    }
    while (number_of_runing_threads_chromosome > 0) {// wait for all the thread
        usleep(100);
    }
    while (number_of_runing_threads_window > 0) {// wait for all the thread
        usleep(100);
    }
}


void prepareForMsa( std::string& referenceGenomeFilePath, std::string& referenceGffFilePath, std::map<std::string, std::string>& sdiFiles,
                    int& maxThread, int & lengthThread, std::string & vcfFix, std::map<std::string, std::string>& parameters,
                    bool& append, int& msaWindowSize, int& msaWindowOverlap, int& minIntron, int & outputPoolSize) {
    --msaWindowSize;
    if( msaWindowSize<=0 ){
        std::cerr << "window size for MSA should be large than 1" << std::endl;
    }
    if( msaWindowOverlap<=0 ){
        std::cerr << "window overlap for MSA should be positive value" << std::endl;
    }
    if( msaWindowOverlap >= msaWindowSize ){
        std::cerr << "window overlap for MSA should be be smaller than window size" << std::endl;
    }
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);
    std::map<std::string, Fasta> referenceGenome;
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    readFastaFile(referenceGenomeFilePath, referenceGenome);
    std::cout << "reference genome sequence reading done" << std::endl;
    std::string regex = get_parameters("cdsParentRegex", parameters);
    readGffFile(referenceGffFilePath, referenceTranscriptHashSet, regex);
    std::cout << "reference genome annotation reading done" << std::endl;

    // remove those transcripts that have no reference chromosome avaliable begin
    std::set<std::string> toRemoveChromosomes;
    for (std::map<std::string, std::vector<Transcript> >::iterator it = referenceTranscriptHashSet.begin();
            it != referenceTranscriptHashSet.end(); ++it) {
        if (referenceGenome.find(it->first) == referenceGenome.end()) {
            toRemoveChromosomes.insert(it->first);
        }
    }
    for (std::set<std::string>::iterator it = toRemoveChromosomes.begin();
         it != toRemoveChromosomes.end(); ++it) {
        if (referenceTranscriptHashSet.find(*it) != referenceTranscriptHashSet.end()) {
            referenceTranscriptHashSet.erase(*it);
        }
    } // remove those transcripts that have no reference chromosome avaliable end

    // remove ORF shift allele from reference dataset begin
    int binSize = 600; // for parallel, how many transcripts to check for each thread
    for ( std::map<std::string, std::vector<Transcript> >::iterator it1=referenceTranscriptHashSet.begin();
            it1!=referenceTranscriptHashSet.end(); ++it1){

//        std::cout << it1->first << " check ORF begin" << std::endl;
        //check ORF begin
        std::atomic_int number_of_runing_threads0(0);
        for (int index=0; index < (it1->second.size()); index+=binSize){

//        for (std::vector<Transcript>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
            bool isThisThreadUnrun = true;
            while (isThisThreadUnrun) {
                if (number_of_runing_threads0 < maxThread) {
                    std::thread t(checkOrfPversionBin, index, std::ref(binSize), std::ref(it1->second), std::ref(referenceGenome), std::ref(nucleotideCodeSubstitutionMatrix), std::ref(minIntron), std::ref(number_of_runing_threads0));
                    ++number_of_runing_threads0;
                    t.detach();
                    isThisThreadUnrun = false;
                    break;
                } else {
                    usleep(10);
                }
            }
        }
        while (number_of_runing_threads0 > 0) {// wait for all the thread
            usleep(100);
        }//check ORF end
//        std::cout << it1->first << " check ORF end" << std::endl;

        std::vector<int> transcriptToRemove;
        for (int index=0; index < (it1->second.size()); ++index){
            if ((it1->second)[index].getIfOrfShift()) {
                transcriptToRemove.push_back(index);
                //std::cout << it1->second[index].getName() << " ORF shift " << it1->second[index].getMetaInformation() << std::endl;
            }//else{
               // std::cout << it1->second[index].getName() << " ORF conserved" << std::endl;
            //}
        }// loop from end to begin
        for(int index=transcriptToRemove.size()-1; index>=0; --index){
            referenceTranscriptHashSet[it1->first].erase(referenceTranscriptHashSet[it1->first].begin()+transcriptToRemove[index]);
        }
    } // remove ORF shift allele from reference data set end
    std::cout << "remove ORF shift allele done" << std::endl;

    // sort reference transcripts and remove overlapped transcripts begin
    for ( std::map<std::string, std::vector<Transcript> >::iterator it1=referenceTranscriptHashSet.begin();
            it1!=referenceTranscriptHashSet.end(); ++it1) {
        std::sort(referenceTranscriptHashSet[it1->first].begin(), referenceTranscriptHashSet[it1->first].end(), [](Transcript a, Transcript b) {
            return a.getPStart() < b.getPStart();
        });
    }
    for ( std::map<std::string, std::vector<Transcript> >::iterator it1=referenceTranscriptHashSet.begin();
            it1!=referenceTranscriptHashSet.end(); ++it1) {
        if( it1->second.size()>0 ){
            std::vector<int> transcriptToRemove;
            Transcript lastTranscript = (it1->second)[0];
            for (int index=1; index < (it1->second).size(); ++index){
                if( lastTranscript.ifOverLapIgnorStrand((it1->second)[index]) ){
                    transcriptToRemove.push_back(index);
                    //std::cout << "delete overlap " << referenceTranscriptHashSet[it1->first][index].getName() << std::endl;
                }else{
                    lastTranscript = (it1->second)[index];
                }
            }
            // loop from end to begin
            for(int index=transcriptToRemove.size()-1; index>=0; --index){
                referenceTranscriptHashSet[it1->first].erase(referenceTranscriptHashSet[it1->first].begin()+transcriptToRemove[index]);
            }
        }
    } // sort reference transcripts and remove overlapped transcripts end
    std::cout << "remove overlapped transcripts done" << std::endl;

    std::map<std::string, std::vector<Window>> transcriptToWindowMap;
    std::map<std::string, std::vector<Window>> chromosomeToWindowMap;
    std::string tempFolder = createdMsaPreFolder(parameters); // create tempFolder

    //create folder for each chromosome begin
    for (std::map<std::string, Fasta>::iterator folderit=referenceGenome.begin();
         folderit!=referenceGenome.end(); ++folderit){
        createdFolder(tempFolder+"/"+folderit->first);
    } //create folder for each chromosome end

    // prepare window begin
    for (std::map<std::string, Fasta>::iterator chrIt=referenceGenome.begin(); chrIt!=referenceGenome.end(); ++chrIt){
        int chrSize = chrIt->second.getSequence().length();
        // transcript windows begin
        std::vector<Window> cdsWindows;
        if( referenceTranscriptHashSet.find(chrIt->first) != referenceTranscriptHashSet.end() ){
            for ( std::vector<Transcript>::iterator transcriptIt = referenceTranscriptHashSet[chrIt->first].begin();
                  transcriptIt!=referenceTranscriptHashSet[chrIt->first].end(); ++transcriptIt ){

                //std::cout << transcriptIt->getName() << " begin to prepare window for this transcript" << std::endl;
                /*if( transcriptIt->getIfOrfShift() ){
                    std::cout << transcriptIt->getName() << " ORF shift" << std::endl;
                }else{
                    std::cout << transcriptIt->getName() << " ORF conserved" << std::endl;
                }*/
                transcriptToWindowMap[transcriptIt->getName()] = std::vector<Window>();
                for ( int i=0; i< transcriptIt->getCdsVector().size(); ++i ){
                    if( i>0 ){
                        // intron region begin
                        int thisRegionStart = (transcriptIt->getCdsVector())[i-1].getEnd()+1;
                        int thisRegionEnd = (transcriptIt->getCdsVector())[i].getStart()-1;

                        bool _postExtand = true;
                        int thinWindowStart = thisRegionStart;
                        int thinWindowEnd = thinWindowStart + msaWindowSize;
                        if( thinWindowEnd >= thisRegionEnd){
                            thinWindowEnd = thisRegionEnd;
                            _postExtand = false;
                        }

                        // the first window in this region does not need pre- extend
                        Window window1 = Window(thinWindowStart, thinWindowEnd, transcriptIt->getName(), CDSINTRON, false, _postExtand);
                        //std::cout << thinWindowStart << " " << thinWindowEnd << " CDSINTRON" << std::endl;
                        cdsWindows.push_back(window1);
                        transcriptToWindowMap[transcriptIt->getName()].push_back(window1);
                        while( thinWindowEnd < thisRegionEnd ){
                            thinWindowStart = thinWindowStart + msaWindowSize + 1;
                            thinWindowEnd = thinWindowStart + msaWindowSize;
                            if( thinWindowEnd >= thisRegionEnd){
                                thinWindowEnd = thisRegionEnd;
                                _postExtand = false;
                            } // the middle window in this region need pre- extend
                            Window window2 = Window(thinWindowStart, thinWindowEnd, transcriptIt->getName(),  CDSINTRON, true, _postExtand);
                            //std::cout << thinWindowStart << " " << thinWindowEnd << " CDSINTRON" << std::endl;
                            cdsWindows.push_back(window2);
                            transcriptToWindowMap[transcriptIt->getName()].push_back(window2);
                        }// intron region end
                    }

                    // for exon begin
                    int thisRegionStart = (transcriptIt->getCdsVector())[i].getStart();
                    int thisRegionEnd = (transcriptIt->getCdsVector())[i].getEnd();

                    bool _postExtand = true;
                    int thinWindowStart = thisRegionStart;
                    int thinWindowEnd = thinWindowStart + msaWindowSize;
                    if( thinWindowEnd >= thisRegionEnd){
                        thinWindowEnd = thisRegionEnd;
                        _postExtand = false;
                    }
                    // the first window in this region does not need pre- extend
                    Window window1 = Window(thinWindowStart, thinWindowEnd, transcriptIt->getName(), CDS, false, _postExtand);
                    //std::cout << thinWindowStart << " " << thinWindowEnd << " CDS" << std::endl;
                    cdsWindows.push_back(window1);
                    transcriptToWindowMap[transcriptIt->getName()].push_back(window1);
                    while( thinWindowEnd < thisRegionEnd ){
                        thinWindowStart = thinWindowStart + msaWindowSize + 1;
                        thinWindowEnd = thinWindowStart + msaWindowSize;
                        if( thinWindowEnd >= thisRegionEnd){
                            thinWindowEnd = thisRegionEnd;
                            _postExtand = false;
                        }
                        // the middle window in this region need pre- extend
                        Window window2 = Window(thinWindowStart, thinWindowEnd, transcriptIt->getName(), CDS, true, _postExtand);
                        //std::cout << thinWindowStart << " " << thinWindowEnd << " CDS" << std::endl;
                        cdsWindows.push_back(window2);
                        transcriptToWindowMap[transcriptIt->getName()].push_back(window2);
                    } // for exon end
                }
                std::sort(transcriptToWindowMap[transcriptIt->getName()].begin(), transcriptToWindowMap[transcriptIt->getName()].end());
            }
        } // transcript windows end

        std::sort(cdsWindows.begin(), cdsWindows.end());

        bool _postExtand = false;
        bool _preExtand = false;
        std::vector<Window> otherWindows;
        Window lastWindow = Window(0, 0, _preExtand, _postExtand);
        for ( std::vector<Window>::iterator winit=cdsWindows.begin();
              winit!=cdsWindows.end(); ++winit){
            if( lastWindow.getEnd() < winit->getStart()-1 ){
                _preExtand = false;
                _postExtand = true;

                int thisRegionStart = lastWindow.getEnd() + 1;
                int thisRegionEnd = (*winit).getStart() - 1;
                int thinWindowStart = thisRegionStart;
                int thinWindowEnd = thinWindowStart + msaWindowSize;
                if( thinWindowEnd >= thisRegionEnd){
                    thinWindowEnd = thisRegionEnd;
                    _postExtand = false;
                }
                Window window1 = Window(thinWindowStart, thinWindowEnd, _preExtand, _postExtand);
                otherWindows.push_back(window1);
                _preExtand = true;
                while( thinWindowEnd < thisRegionEnd ){
                    thinWindowStart = thinWindowStart + msaWindowSize + 1;
                    thinWindowEnd = thinWindowStart + msaWindowSize;
                    if( thinWindowEnd >= thisRegionEnd){
                        thinWindowEnd = thisRegionEnd;
                        _postExtand = false;
                    }
                    Window window2 = Window(thinWindowStart, thinWindowEnd, _preExtand, _postExtand);
                    otherWindows.push_back(window2);
                }
            }
            lastWindow = (*winit);
        }

        // deal with the region after last gene begin
        _preExtand = false;
        _postExtand = true;
        int start = cdsWindows[cdsWindows.size()-1].getEnd() + 1;
        int end = start+msaWindowSize;
        if( end >= chrSize ){
            end = chrSize;
            _postExtand = false;
        }

        Window window1 = Window(start, end, _preExtand, _postExtand);
        otherWindows.push_back(window1);
        _preExtand = true;
        while( end < chrSize ){
            start = start + msaWindowSize + 1;
            end = start + msaWindowSize;
            if( end >= chrSize ){
                end = chrSize;
                _postExtand = false;
            }
            Window window2 = Window(start, end, _preExtand, _postExtand);
            otherWindows.push_back(window2);
        }// deal with the region after last gene end


        for ( std::vector<Window>::iterator winit=cdsWindows.begin();
              winit!=cdsWindows.end(); ++winit){
            otherWindows.push_back(*winit);
        }
        cdsWindows.empty();
        std::sort(otherWindows.begin(), otherWindows.end());
        chromosomeToWindowMap[chrIt->first] = otherWindows;
    }// prepare window end

    if( append ){
        // just add new sequence to avaliable multiple sequence fasta files
    }else{
        std::cout << "output reference sequence" << std::endl;
        outputReferenceSequence(chromosomeToWindowMap, tempFolder, maxThread, msaWindowOverlap, referenceGenome);
    }
    // begin to deal with non-reference accession/line
    std::cout << "begin to deal with non-reference accession/line" << std::endl;

    FileWriteCatch fileWriteCatch;
    int numberInCatch = 0;

    for (std::map<std::string, std::string >::iterator it0=sdiFiles.begin(); it0!=sdiFiles.end(); ++it0){
        std::string accessionName = it0->first;

        std::map<std::string, std::vector<Variant> > variantsMaps;
        readSdiFile(it0->second, variantsMaps, vcfFix, referenceGenome);
        std::cout << it0->first << " sdi reading done" << std::endl;
        std::map<std::string, Fasta> targetGenome;
        getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome, parameters);
        std::cout << it0->first << " pseudo genome sequence done" << std::endl;
        std::map<std::string, Transcript> targetTranscriptsHashMap;// transcriptName Transcript
     // this function is very slow
        annotationLiftOverPversion(referenceTranscriptHashSet, targetTranscriptsHashMap, variantsMaps, targetGenome, referenceGenome, parameters, minIntron, maxThread);
        /*
         std::cout << it0->first << " annotationLiftOver done. Gene structure alignment begin. maxThread: " << maxThread << std::endl;
        //std::vector<std::thread> threads;

        std::atomic_int number_of_runing_threads2(0);
        for (std::map<std::string, std::vector<Transcript> >::iterator it = referenceTranscriptHashSet.begin();
                it != referenceTranscriptHashSet.end(); ++it) {
            for (std::vector<Transcript>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                Transcript *it3 = &targetTranscriptsHashMap[it2->getName()];
                if ((*it3).getIfOrfShift()) {
                    bool isThisThreadUnrun = true;
                    while (isThisThreadUnrun) {
                        if (number_of_runing_threads2 < (maxThread))  {
                            std::thread t(transcriptRealignmentPversion, std::ref(*it3), std::ref(*it2),
                                          std::ref(nucleotideCodeSubstitutionMatrix), std::ref(targetGenome), std::ref(number_of_runing_threads2),
                                          std::ref(referenceGenome), it->first, std::ref(targetTranscriptsHashMap), std::ref(lengthThread), std::ref(parameters), std::ref(minIntron));
                            ++number_of_runing_threads2;
                            t.detach();
                            isThisThreadUnrun = false;
                            break;
                        } else {
                            usleep(10000);
                        }
                    }
                }
            }
        }
        while (number_of_runing_threads2 > 0) {// wait for all the thread
            usleep(100000);
        }
        std::cerr << it0->first << " realignment done" << std::endl;
*/
        outputNonreferenceSequence( chromosomeToWindowMap, accessionName, tempFolder, maxThread, targetTranscriptsHashMap, targetGenome,
            variantsMaps, msaWindowOverlap, transcriptToWindowMap, fileWriteCatch );
        std::cerr << it0->first << " window output done" << std::endl;

        numberInCatch++;
        if( numberInCatch%outputPoolSize == 0 ){
            fileWriteCatch.writeOut(maxThread);
        }
    }
    fileWriteCatch.writeOut(maxThread);
}










































































































































void constructSdiFromMsa_v_beta(std::vector<std::string>& chromosomes, std::set<std::string>& accessionNames, std::string& folder, std::string& outputFolder, std::string & referenceGenomeFilePath,
                                std::map<std::string, std::string>& sdiFiles, std::string & vcfFix, std::map<std::string, std::string>& parameters){
    std::map<std::string, Fasta> referenceGenome;
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    readFastaFile(referenceGenomeFilePath, referenceGenome);

    int len = folder.length();
    if (folder[len-1] != '/'){
        folder = folder + "/";
    }
    createdFolder("newSdis");

    for( int chromosomei=0; chromosomei<chromosomes.size(); ++chromosomei ){
        // read file begin
        std::string chrName = chromosomes[chromosomei];

        createdFolder("newSdis/" + chrName);
        std::vector<std::string> files;
        std::string thisFoler = folder + chrName;
        getMafftResultListFromFolder(thisFoler, files);
        std::vector<MsaFileRecord> msaFileRecords;
        std::regex reg("^(\\d+)_(\\d+)\\.mafft");
        for(std::vector<std::string>::iterator itf=files.begin();
            itf!=files.end(); ++ itf ){
            std::smatch match;
            regex_search((*itf), match, reg);
            if( !match.empty() ){
                MsaFileRecord msaFileRecord(stoi(match[1]), stoi(match[2]));
                std::string msaFileLocation = thisFoler + "/"+(*itf);
                msaFileRead(msaFileRecord, msaFileLocation, accessionNames);
                msaFileRecords.push_back(msaFileRecord);
            }
        }
        std::sort(msaFileRecords.begin(), msaFileRecords.end());
        // read files end
        std::cout << "get file list done" << std::endl;
        for( std::set<std::string>::iterator itName=accessionNames.begin(); itName!=accessionNames.end(); ++itName ){

            std::map<std::string, std::vector<Variant> > variantsMaps;
            readSdiFile(sdiFiles[*itName], variantsMaps, vcfFix, referenceGenome);
            std::map<std::string, Fasta> targetGenome;
            getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome, parameters);
            std::string targetGenomeThisChrSequence = targetGenome[chrName].getSequence();
            int thisTargetChromosomeLength = targetGenomeThisChrSequence.length();

            FirstLastList sdiRecords;
//            std::vector<Variant>  newVariants;
//            std::vector<Variant> newSdiRecords;
            std::vector<TwoSeqOfMsaResult> twoSeqOfMsaResults;


            //delete the extend sequence only keep the wanted region begin
            for( int i=0; i<msaFileRecords.size(); ++i ) {
                MsaFileRecord msaFileRecord = msaFileRecords[i];
                int start = msaFileRecord.getStart();
                int end = msaFileRecord.getEnd();
//                --end;
                int targetStart = 0;
                int targetEnd = 0;

                MsaSingleRecord refMsaSingleRecord = msaFileRecord.getMsaSingleRecordRecords()["ref"];
                MsaSingleRecord targetMsaSingleRecord = msaFileRecord.getMsaSingleRecordRecords()[*itName];
                int msaRefStart = refMsaSingleRecord.getStart();
                int msaTargetStart = targetMsaSingleRecord.getStart();

                //delete the extend sequence only keep the wanted region
                int refLetterNumber = 0;
                int targetLetterNumber = 0;
                std::stringstream refSeq;
                std::stringstream targetSeq;


                //for debug begin

                std::string thisResultSeqS = targetMsaSingleRecord.getSequence();
                thisResultSeqS.erase(std::remove(thisResultSeqS.begin(), thisResultSeqS.end(), '-'),
                                     thisResultSeqS.end());
                //if( targetMsaSingleRecord.getStart() < targetMsaSingleRecord.getEnd() ){

                //}else{
                //}
                    std::cout << thisResultSeqS << " " << thisResultSeqS.size() << " " << targetMsaSingleRecord.getStart() << " " << targetMsaSingleRecord.getEnd() << std::endl <<
                              refMsaSingleRecord.getStart() << " " << refMsaSingleRecord.getEnd() << std::endl;
                assert( thisResultSeqS.size() == (targetMsaSingleRecord.getEnd()-targetMsaSingleRecord.getStart()+1) );

                std::string thisRefSeqS = refMsaSingleRecord.getSequence();
                thisRefSeqS.erase(std::remove(thisRefSeqS.begin(), thisRefSeqS.end(), '-'), thisRefSeqS.end());

                int thisSubStart1 = refMsaSingleRecord.getStart();
                int thisSubEnd1 = refMsaSingleRecord.getEnd();

                int thisSubStart2 = targetMsaSingleRecord.getStart();
                int thisSubEnd2 = targetMsaSingleRecord.getEnd();
                std::string correctRefSeqs = getSubsequence(referenceGenome, chrName, thisSubStart1, thisSubEnd1);
                std::string correctResultSeqs = getSubsequence(targetGenome, chrName, thisSubStart2, thisSubEnd2);
                if( targetMsaSingleRecord.getStart() > targetMsaSingleRecord.getEnd()){
                    correctResultSeqs="";
                }
//                std::cout << "1696 ref: " << thisRefSeqS << " correct: " << correctRefSeqs << " " << refMsaSingleRecord.getStart() << " " << refMsaSingleRecord.getEnd() << std::endl;
//                std::cout << "1697 result:" << thisResultSeqS << "correct: " << correctResultSeqs << " " << targetMsaSingleRecord.getStart() << " " << targetMsaSingleRecord.getEnd() << std::endl;
                assert(thisRefSeqS.compare(correctRefSeqs) == 0);
                assert(thisResultSeqS.compare(correctResultSeqs) == 0);

                //for debug end

                for( int ai=0; ai<refMsaSingleRecord.getSequence().length(); ++ai){
                    if( refMsaSingleRecord.getSequence()[ai] != '-' ){
                        ++refLetterNumber;
                    }
                    if( targetMsaSingleRecord.getSequence()[ai] != '-' ){
                        ++targetLetterNumber;
                    }
                    if( start == (msaRefStart+refLetterNumber-1) && refMsaSingleRecord.getSequence()[ai] != '-' ){
                        targetStart = msaTargetStart+targetLetterNumber-1;
                        if( targetMsaSingleRecord.getSequence()[ai] == '-' ){
                            ++targetStart;
                        }
                    }
                    if( (msaRefStart+refLetterNumber-1) == end && refMsaSingleRecord.getSequence()[ai] != '-'){
                        targetEnd = msaTargetStart + targetLetterNumber -1;
                    }
                    if( start <= (msaRefStart+refLetterNumber-1) && (msaRefStart+refLetterNumber-1)<end ){
                        refSeq << refMsaSingleRecord.getSequence()[ai];
                        targetSeq << targetMsaSingleRecord.getSequence()[ai];
                    }else if ( (msaRefStart+refLetterNumber-1)==end &&  refMsaSingleRecord.getSequence()[ai] != '-' ){
                        refSeq << refMsaSingleRecord.getSequence()[ai];
                        targetSeq << targetMsaSingleRecord.getSequence()[ai];
                    }
                }
                if( ( i == msaFileRecords.size()-1 ) && targetEnd ==0){
                    targetEnd = thisTargetChromosomeLength;
                    std::cout << "1724 never never run here" << std::endl;
                    assert(0);
                }
                std::string refSequence = refSeq.str();
                transform(refSequence.begin(), refSequence.end(), refSequence.begin(),::toupper);
                std::string targetSequence = targetSeq.str();
                transform(targetSequence.begin(), targetSequence.end(), targetSequence.begin(),::toupper);
                if( targetStart==0 && targetEnd==0  ){
                    targetEnd = -1;
                }
                TwoSeqOfMsaResult twoSeqOfMsaResult( start, end, refSequence, targetStart, targetEnd, targetSequence );
                twoSeqOfMsaResults.push_back(twoSeqOfMsaResult);
//                std::cout << start << " " << end << std::endl;
//                std::cout << targetStart << " " << targetEnd << std::endl;
                assert(start <= end);
                std::string thisResultSeq = targetSequence;
                thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());
//                std::cout << targetStart << " " << targetEnd << " " << thisResultSeq << std::endl;
                assert( thisResultSeq.size() == (targetEnd-targetStart+1) );


                std::string thisResultSeq2 = twoSeqOfMsaResult.getResultSeq();
                thisResultSeq2.erase(std::remove(thisResultSeq2.begin(), thisResultSeq2.end(), '-'), thisResultSeq2.end());
                assert( thisResultSeq2.size() == (twoSeqOfMsaResult.getResultEnd()-twoSeqOfMsaResult.getResultStart()+1) );
                std::string thisRefSeq = twoSeqOfMsaResult.getRefSeq();
                thisResultSeq = targetSequence;
                thisRefSeq.erase(std::remove(thisRefSeq.begin(), thisRefSeq.end(), '-'), thisRefSeq.end());
                thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());

                int thisSubStart_1 = twoSeqOfMsaResult.getRefStart();
                int thisSubEnd_1 = twoSeqOfMsaResult.getRefEnd();

                int thisSubStart_2 = twoSeqOfMsaResult.getResultStart();
                int thisSubEnd_2 = twoSeqOfMsaResult.getResultEnd();

                std::string correctRefSeq = getSubsequence(referenceGenome, chrName, thisSubStart_1, thisSubEnd_1);
                std::string correctResultSeq = getSubsequence(targetGenome, chrName, thisSubStart_2, thisSubEnd_2);
                if( twoSeqOfMsaResult.getResultStart() > twoSeqOfMsaResult.getResultEnd() ){
                    correctResultSeq="";
                }
//                std::cout << "1754 ref: " << thisRefSeq << " correct: " << correctRefSeq << " " << twoSeqOfMsaResult.getRefStart() << " " << twoSeqOfMsaResult.getRefEnd() << std::endl;
//                std::cout << "1755 result:" << thisResultSeq << "correct: " << correctResultSeq << " " << twoSeqOfMsaResult.getResultStart() << " " << twoSeqOfMsaResult.getResultEnd() << std::endl;
                assert(thisRefSeq.compare(correctRefSeq) == 0);
                assert(thisResultSeq.compare(correctResultSeq) == 0);
            }
            std::sort(twoSeqOfMsaResults.begin(), twoSeqOfMsaResults.end());
            //delete the extend sequence end
            std::cout << (*itName) << " delete the extend sequence  done" << std::endl;

            for ( int iIndex = 0; iIndex < twoSeqOfMsaResults.size(); ++iIndex) {
                std::string thisResultSeq2 = twoSeqOfMsaResults[iIndex].getResultSeq();
                thisResultSeq2.erase(std::remove(thisResultSeq2.begin(), thisResultSeq2.end(), '-'), thisResultSeq2.end());
                assert( thisResultSeq2.size() == (twoSeqOfMsaResults[iIndex].getResultEnd()-twoSeqOfMsaResults[iIndex].getResultStart()+1) );
                std::string thisRefSeq = twoSeqOfMsaResults[iIndex].getRefSeq();
                std::string thisResultSeq = twoSeqOfMsaResults[iIndex].getResultSeq();
                thisRefSeq.erase(std::remove(thisRefSeq.begin(), thisRefSeq.end(), '-'), thisRefSeq.end());
                thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());
                int thisSubStart1 = twoSeqOfMsaResults[iIndex].getRefStart();
                int thisSubEnd1 = twoSeqOfMsaResults[iIndex].getRefEnd();

                int thisSubStart2 = twoSeqOfMsaResults[iIndex].getResultStart();
                int thisSubEnd2 = twoSeqOfMsaResults[iIndex].getResultEnd();
                std::string correctRefSeq = getSubsequence(referenceGenome, chrName, thisSubStart1, thisSubEnd1);
                std::string correctResultSeq = getSubsequence(targetGenome, chrName, thisSubStart2, thisSubEnd2);
                if( twoSeqOfMsaResults[iIndex].getResultStart() > twoSeqOfMsaResults[iIndex].getResultEnd() ){
                    correctResultSeq="";
                }
//                std::cout << "1737ref: " << thisRefSeq << " correct: " << correctRefSeq << " " << twoSeqOfMsaResults[iIndex].getRefStart() << " " << twoSeqOfMsaResults[iIndex].getRefEnd() << std::endl;
//                std::cout << "1738 result:" << thisResultSeq << "correct: " << correctResultSeq << " " << twoSeqOfMsaResults[iIndex].getResultStart() << " " << twoSeqOfMsaResults[iIndex].getResultEnd() << std::endl;
                assert(thisRefSeq.compare(correctRefSeq) == 0);
                assert(thisResultSeq.compare(correctResultSeq) == 0);
            }

            // merge overlapped neighbor window begin
            for( int j=1; j<twoSeqOfMsaResults.size(); j++ ){
                //if( twoSeqOfMsaResults.size() > 1 ){
                if( j<1 ){
                    j=1;
                }
                std::string thisResultSeq = twoSeqOfMsaResults[j-1].getResultSeq();
                thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());
                assert( thisResultSeq.size() == (twoSeqOfMsaResults[j-1].getResultEnd()-twoSeqOfMsaResults[j-1].getResultStart()+1) );

                thisResultSeq = twoSeqOfMsaResults[j].getResultSeq();
                thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());
                assert( thisResultSeq.size() == (twoSeqOfMsaResults[j].getResultEnd()-twoSeqOfMsaResults[j].getResultStart()+1) );

                if( twoSeqOfMsaResults[j].getResultEnd() < twoSeqOfMsaResults[j].getResultStart() ){
                    std::cout <<"1732 to be complete , This should never appear. If you find it, please contact me." << std::endl;
//                        std::cout << "j: " << j <<std::endl;
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::string thisResultSeq = twoSeqOfMsaResults[j].getResultSeq();
//                        std::cout << "j: " << thisResultSeq <<std::endl;

                    if ( (twoSeqOfMsaResults[j-1].getRefEnd() >= twoSeqOfMsaResults[j].getRefStart()) ){
                        int overLapLength = twoSeqOfMsaResults[j-1].getRefEnd() - twoSeqOfMsaResults[j].getRefStart();
                        overLapLength++;
                        std::string newRef = twoSeqOfMsaResults[j].getRefSeq();
                        int corrected = 0;
                        for( int t=0; t<newRef.length(); t++ ){
                            if(corrected<overLapLength && newRef[t]!='-' ){
                                newRef=newRef.substr(0, t)+"-"+newRef.substr(t+1, newRef.length()-t-1);
                                corrected++;
                            }
                        }
                        int largerRefEnd;
                        if( twoSeqOfMsaResults[j].getRefEnd() > twoSeqOfMsaResults[j-1].getRefEnd() ){
                            largerRefEnd = twoSeqOfMsaResults[j].getRefEnd();
                        }else{
                            largerRefEnd = twoSeqOfMsaResults[j-1].getRefEnd();
                        }
//
//                            std::cout << "j: " << j <<std::endl;
//                            std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                            std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;

                        twoSeqOfMsaResults[j-1].setRefEnd(largerRefEnd);
                        twoSeqOfMsaResults[j-1].setRefSeq(twoSeqOfMsaResults[j-1].getRefSeq()+newRef);

                        twoSeqOfMsaResults[j-1].setResultSeq(twoSeqOfMsaResults[j-1].getResultSeq()+twoSeqOfMsaResults[j].getResultSeq());
//                            std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                            std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
                        twoSeqOfMsaResults.erase(twoSeqOfMsaResults.begin()+j);
//                            std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                            std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                            std::cout <<"1762 to be complete , This should never appear. If you find it, please contact me." << std::endl;
                        --j;
                        --j;
                    } else if ( twoSeqOfMsaResults[j-1].getResultEnd() < twoSeqOfMsaResults[j-1].getResultStart()){
                        twoSeqOfMsaResults[j].setResultStart(twoSeqOfMsaResults[j-1].getResultStart());
                        twoSeqOfMsaResults[j].setResultEnd(twoSeqOfMsaResults[j-1].getResultStart()-1);
                    } else {
                        twoSeqOfMsaResults[j].setResultStart(twoSeqOfMsaResults[j-1].getResultEnd());
                        twoSeqOfMsaResults[j].setResultEnd(twoSeqOfMsaResults[j-1].getResultEnd()-1);
                    }
                } else if ( (twoSeqOfMsaResults[j-1].getResultEnd() >= twoSeqOfMsaResults[j].getResultStart()) && (twoSeqOfMsaResults[j-1].getRefEnd() >= twoSeqOfMsaResults[j].getRefStart())
                            && twoSeqOfMsaResults[j-1].getResultStart() <= twoSeqOfMsaResults[j-1].getResultEnd() ){
                    int overLapLength = twoSeqOfMsaResults[j-1].getResultEnd() - twoSeqOfMsaResults[j].getResultStart();
                    ++overLapLength;
                    std::string newSequence = twoSeqOfMsaResults[j].getResultSeq();
                    int corrected = 0;
                    for( int t=0; t<newSequence.length(); t++ ){
                        if(corrected<overLapLength && newSequence[t]!='-' ){
                            newSequence=newSequence.substr(0, t)+"-"+newSequence.substr(t+1, newSequence.length()-t-1);
                            corrected++;
                        }
                    }
                    int largerEnd;
                    if( twoSeqOfMsaResults[j].getResultEnd() > twoSeqOfMsaResults[j-1].getResultEnd() ){
                        largerEnd = twoSeqOfMsaResults[j].getResultEnd();
                    } else {
                        largerEnd = twoSeqOfMsaResults[j-1].getResultEnd();
                    }

                    overLapLength = twoSeqOfMsaResults[j-1].getRefEnd() - twoSeqOfMsaResults[j].getRefStart();
                    overLapLength++;
                    std::string newRef = twoSeqOfMsaResults[j].getRefSeq();
                    corrected = 0;
                    for( int t=0; t<newRef.length(); t++ ){
                        if(corrected<overLapLength && newRef[t]!='-' ){
                            newRef=newRef.substr(0, t)+"-"+newRef.substr(t+1, newRef.length()-t-1);
                            corrected++;
                        }
                    }
                    int largerRefEnd;
                    if( twoSeqOfMsaResults[j].getRefEnd() > twoSeqOfMsaResults[j-1].getRefEnd() ){
                        largerRefEnd = twoSeqOfMsaResults[j].getRefEnd();
                    }else{
                        largerRefEnd = twoSeqOfMsaResults[j-1].getRefEnd();
                    }
//                        std::cout << "j: " << j <<std::endl;
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;

                    twoSeqOfMsaResults[j-1].setRefEnd(largerRefEnd);
                    twoSeqOfMsaResults[j-1].setRefSeq(twoSeqOfMsaResults[j-1].getRefSeq()+newRef);

                    twoSeqOfMsaResults[j-1].setResultEnd(largerEnd);
                    twoSeqOfMsaResults[j-1].setResultSeq(twoSeqOfMsaResults[j-1].getResultSeq()+newSequence);
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
                    twoSeqOfMsaResults.erase(twoSeqOfMsaResults.begin()+j);
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout <<"1777 to be complete , This should never appear. If you find it, please contact me." << std::endl;
                    --j;
                    --j;
                }else if ( (twoSeqOfMsaResults[j-1].getResultStart() >= twoSeqOfMsaResults[j].getResultStart()) && (twoSeqOfMsaResults[j-1].getRefEnd() >= twoSeqOfMsaResults[j].getRefStart())
                           && twoSeqOfMsaResults[j-1].getResultStart() > twoSeqOfMsaResults[j-1].getResultEnd() ){
                    int overLapLength = twoSeqOfMsaResults[j-1].getRefEnd() - twoSeqOfMsaResults[j].getRefStart();
                    overLapLength++;
                    std::string newRef = twoSeqOfMsaResults[j].getRefSeq();
                    int corrected = 0;
                    for( int t=0; t<newRef.length(); t++ ){
                        if(corrected<overLapLength && newRef[t]!='-' ){
                            newRef=newRef.substr(0, t)+"-"+newRef.substr(t+1, newRef.length()-t-1);
                            corrected++;
                        }
                    }
                    int largerRefEnd;
                    if( twoSeqOfMsaResults[j].getRefEnd() > twoSeqOfMsaResults[j-1].getRefEnd() ){
                        largerRefEnd = twoSeqOfMsaResults[j].getRefEnd();
                    }else{
                        largerRefEnd = twoSeqOfMsaResults[j-1].getRefEnd();
                    }

                    overLapLength = twoSeqOfMsaResults[j-1].getResultStart() - twoSeqOfMsaResults[j].getResultStart();
                    ++overLapLength;
                    std::string newSequence = twoSeqOfMsaResults[j].getResultSeq();
                    corrected = 0;
                    for( int t=0; t<newSequence.length(); t++ ){
                        if(corrected<overLapLength && newSequence[t]!='-' ){
                            newSequence=newSequence.substr(0, t)+"-"+newSequence.substr(t+1, newSequence.length()-t-1);
                            corrected++;
                        }
                    }
                    int end = twoSeqOfMsaResults[j].getResultEnd(); // the current one is not empty
                    int start =twoSeqOfMsaResults[j].getResultStart() + overLapLength;


                    twoSeqOfMsaResults[j-1].setRefEnd(largerRefEnd);
                    twoSeqOfMsaResults[j-1].setRefSeq(twoSeqOfMsaResults[j-1].getRefSeq()+newRef);

                    if( start < end ){
                        twoSeqOfMsaResults[j-1].setResultEnd(end);
                        twoSeqOfMsaResults[j-1].setResultStart(start);
                        twoSeqOfMsaResults[j-1].setResultSeq(twoSeqOfMsaResults[j-1].getResultSeq()+newSequence);
                    } else {
                        twoSeqOfMsaResults[j-1].setResultSeq(twoSeqOfMsaResults[j-1].getResultSeq()+twoSeqOfMsaResults[j].getResultSeq());
                    }

//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
                    twoSeqOfMsaResults.erase(twoSeqOfMsaResults.begin()+j);

//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout <<"line 1952" << std::endl;
                    --j;
                    --j;

                } else if( (twoSeqOfMsaResults[j-1].getResultEnd() >= twoSeqOfMsaResults[j].getResultStart()) && twoSeqOfMsaResults[j-1].getResultStart() <= twoSeqOfMsaResults[j-1].getResultEnd() ) {

                    int overLapLength =  twoSeqOfMsaResults[j - 1].getResultEnd() - twoSeqOfMsaResults[j].getResultStart();
                    ++overLapLength;
                    std::string newSequence = twoSeqOfMsaResults[j].getResultSeq();
                    int corrected = 0;
                    for (int t = 0; t < newSequence.length(); ++t) {
                        if (corrected < overLapLength && newSequence[t] != '-') {
                            newSequence = newSequence.substr(0, t) + "-" + newSequence.substr(t + 1, newSequence.length() - t - 1);
                            corrected++;
                        }
                    }
//
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
                    twoSeqOfMsaResults[j].setResultSeq(newSequence);
                    twoSeqOfMsaResults[j].setResultStart(twoSeqOfMsaResults[j].getResultStart()+overLapLength);
                    if( twoSeqOfMsaResults[j].getResultStart() > twoSeqOfMsaResults[j].getResultEnd() ){
                        twoSeqOfMsaResults[j].setResultStart( twoSeqOfMsaResults[j-1].getResultEnd() );
                        twoSeqOfMsaResults[j].setResultEnd(twoSeqOfMsaResults[j].getResultStart() - 1);
                    }
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout <<"line 1985" << std::endl;
                } else if( (twoSeqOfMsaResults[j-1].getResultStart()>=twoSeqOfMsaResults[j].getResultStart()) && twoSeqOfMsaResults[j-1].getResultStart()>twoSeqOfMsaResults[j-1].getResultEnd() ){

                    int overLapLength =  twoSeqOfMsaResults[j - 1].getResultStart() - twoSeqOfMsaResults[j].getResultStart();

                    ++overLapLength;
                    std::string newSequence = twoSeqOfMsaResults[j].getResultSeq();
                    int corrected = 0;
                    for (int t = 0; t < newSequence.length(); ++t) {
                        if (corrected < overLapLength && newSequence[t] != '-') {
                            newSequence = newSequence.substr(0, t) + "-" + newSequence.substr(t + 1, newSequence.length() - t - 1);
                            corrected++;
                        }
                    }
//
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
                    twoSeqOfMsaResults[j].setResultSeq(newSequence);
                    twoSeqOfMsaResults[j].setResultStart(twoSeqOfMsaResults[j].getResultStart()+overLapLength);
                    if( twoSeqOfMsaResults[j].getResultStart() > twoSeqOfMsaResults[j].getResultEnd() ){
                        twoSeqOfMsaResults[j].setResultStart( twoSeqOfMsaResults[j-1].getResultStart() );
                        twoSeqOfMsaResults[j].setResultEnd(twoSeqOfMsaResults[j].getResultStart() - 1);
                    }

//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout <<"merged 2008" << std::endl;

                }else if ( (twoSeqOfMsaResults[j-1].getRefEnd() >= twoSeqOfMsaResults[j].getRefStart()) ){
                    int overLapLength = twoSeqOfMsaResults[j-1].getRefEnd() - twoSeqOfMsaResults[j].getRefStart();
                    overLapLength++;
                    std::string newRef = twoSeqOfMsaResults[j].getRefSeq();
                    int corrected = 0;
                    for( int t=0; t<newRef.length(); t++ ){
                        if(corrected<overLapLength && newRef[t]!='-'){
                            newRef=newRef.substr(0, t)+"-"+newRef.substr(t+1, newRef.length()-t-1);
                            corrected++;
                        }
                    }
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
                    twoSeqOfMsaResults[j].setRefSeq(newRef);
                    twoSeqOfMsaResults[j].setRefStart(twoSeqOfMsaResults[j].getRefStart()+overLapLength);
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout <<"1824 to be complete , This should never appear. If you find it, please contact me." << std::endl;
                }else if( twoSeqOfMsaResults[j].getResultStart() < twoSeqOfMsaResults[j-1].getResultStart() ){
//                        std::cout << "j: " << j <<std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << " " << twoSeqOfMsaResults[j-1].getRefStart() <<" "<<  twoSeqOfMsaResults[j-1].getRefEnd() << std::endl;
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << " " << twoSeqOfMsaResults[j].getRefStart() << " " << twoSeqOfMsaResults[j].getRefEnd() << std::endl;
                    std::cout <<"1987 to be complete , This should never appear. If you find it, please contact me." << std::endl; // this should never run
                } else {
//                        std::cout << "j: " << j <<std::endl;
//                        std::cout << "j - 1: " << twoSeqOfMsaResults[j-1].getResultStart() << " " << twoSeqOfMsaResults[j-1].getResultEnd() << std::endl;
//                        std::cout << "j: " << twoSeqOfMsaResults[j].getResultStart() << " " << twoSeqOfMsaResults[j].getResultEnd() << std::endl;
//                        std::cout <<"1781 to be complete. Check it there any problem." << std::endl;
                }
                //}
            } // merge overlapped neighbor window end
            std::sort(twoSeqOfMsaResults.begin(), twoSeqOfMsaResults.end());

            for ( int iIndex = 0; iIndex < twoSeqOfMsaResults.size(); ++iIndex) {

                std::string thisResultSeq2 = twoSeqOfMsaResults[iIndex].getResultSeq();
                thisResultSeq2.erase(std::remove(thisResultSeq2.begin(), thisResultSeq2.end(), '-'), thisResultSeq2.end());
                assert( thisResultSeq2.size() == (twoSeqOfMsaResults[iIndex].getResultEnd()-twoSeqOfMsaResults[iIndex].getResultStart()+1) );

                std::string thisRefSeq = twoSeqOfMsaResults[iIndex].getRefSeq();
                std::string thisResultSeq = twoSeqOfMsaResults[iIndex].getResultSeq();

                thisRefSeq.erase(std::remove(thisRefSeq.begin(), thisRefSeq.end(), '-'), thisRefSeq.end());
                thisResultSeq.erase(std::remove(thisResultSeq.begin(), thisResultSeq.end(), '-'), thisResultSeq.end());

                int thisSubStart1 = twoSeqOfMsaResults[iIndex].getRefStart();
                int thisSubEnd1 = twoSeqOfMsaResults[iIndex].getRefEnd();

                int thisSubStart2 = twoSeqOfMsaResults[iIndex].getResultStart();
                int thisSubEnd2 = twoSeqOfMsaResults[iIndex].getResultEnd();

                std::string correctRefSeq = getSubsequence(referenceGenome, chrName, thisSubStart1, thisSubEnd1);
                std::string correctResultSeq = getSubsequence(targetGenome, chrName, thisSubStart2, thisSubEnd2);

                if( twoSeqOfMsaResults[iIndex].getResultStart() > twoSeqOfMsaResults[iIndex].getResultEnd() ){
                    correctResultSeq="";
                }
//                std::cout << "1986 ref: " << thisRefSeq << " correct: " << correctRefSeq << " " << twoSeqOfMsaResults[iIndex].getRefStart() << " " << twoSeqOfMsaResults[iIndex].getRefEnd() << std::endl;
//                std::cout << "1987 result:" << thisResultSeq << "correct: " << correctResultSeq << " " << twoSeqOfMsaResults[iIndex].getResultStart() << " " << twoSeqOfMsaResults[iIndex].getResultEnd() << std::endl;
                assert(thisRefSeq.compare(correctRefSeq) == 0);
                assert(thisResultSeq.compare(correctResultSeq) == 0);
            }

            std::cout << (*itName) << " begin: insert the SDI record before the first window" << std::endl;
            if( twoSeqOfMsaResults[0].getRefStart()>1){//begin: insert the SDI record before the first window
                int thisSubStart1 = 1;
                int thisSubEnd1 = twoSeqOfMsaResults[0].getRefStart()-1;
                std::string ori = getSubsequence(referenceGenome, chrName, thisSubStart1, thisSubEnd1);
                std::string result;
                if(twoSeqOfMsaResults[0].getResultStart()>1){
                    int thisSubStart2 = 1;
                    int thisSubEnd2 = twoSeqOfMsaResults[0].getResultStart()-1;
                    result = getSubsequence( targetGenome, chrName, thisSubStart2, thisSubEnd2);
                }else{
                    result = "-";
                }
                int position = 1;
                Variant mapSingleRecord = Variant(chrName, position, ori, result);

                if( ori.compare(result) != 0 ){
                    Data* data = new Data(mapSingleRecord);
                    sdiRecords.insertLast(data);
                }
            }else if(twoSeqOfMsaResults[0].getRefStart()==1 && twoSeqOfMsaResults[0].getResultStart()>1){
                int thisSubStart2 = 1;
                int thisSubEnd2 = twoSeqOfMsaResults[0].getResultStart()-1;
                std::string result = getSubsequence(targetGenome, chrName, thisSubStart2, thisSubEnd2);
                int position = 1;
                std::string ori = "-";
                Variant mapSingleRecord = Variant(chrName, position, ori, result);
                Data* data = new Data(mapSingleRecord);
                sdiRecords.insertLast(data);
            } //end: insert the first SDI record before the first windows

            std::cout << (*itName) << " begin: insert the windows region SDI record" << std::endl;
            //begin: insert the windows region SDI record
            int i=0;
            while(i<(twoSeqOfMsaResults.size()-1)) {
                std::string refSeq = twoSeqOfMsaResults[i].getRefSeq();
                std::string resultSeq = twoSeqOfMsaResults[i].getResultSeq();
                //std::cout << "ref:       " << refSeq << std::endl;
                //std::cout << "resultSeq: " << resultSeq << std::endl;
                int refLetterNumber = 0;
                for (int ai = 0; ai < refSeq.length(); ai++) {
                    if (refSeq[ai] != '-') {
                        ++refLetterNumber;
                    }
                    if (refSeq[ai] != resultSeq[ai]) {
                        if (refSeq[ai] != resultSeq[ai]) {
                            if (resultSeq[ai] == '-') {
//                            std::cout << "1754" << std::endl;
                                int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber - 1;
                                std::string ori(1, refSeq[ai]);
                                std::string result = "-";
                                Variant mapSingleRecord = Variant(chrName, position, ori, result);
                                Data* data = new Data(mapSingleRecord);
                                sdiRecords.insertLast(data);
//                            std::cout << "1761" << std::endl;
                            } else if (refSeq[ai] == '-') {
//                            std::cout << "1763" << std::endl;
                                if (sdiRecords.getLast() == NULL) {
//                                std::cout << "1755" << std::endl;
                                    int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                                    std::string ori = "-";
                                    std::string result(1, resultSeq[ai]);
                                    Variant mapSingleRecord = Variant(chrName, position, ori, result);
                                    Data* data = new Data(mapSingleRecord);
                                    sdiRecords.insertLast(data);
//                                std::cout << "1772" << std::endl;
                                } else {
//                                std::cout << "1774" << std::endl;
                                    if (NULL!=sdiRecords.getLast() && sdiRecords.getLast()->getMapSingleRecord().getPosition()==(twoSeqOfMsaResults[i].getRefStart() + refLetterNumber)
                                        && sdiRecords.getLast()->getMapSingleRecord().getChanginglength() > 0 &&
                                        sdiRecords.getLast()->getMapSingleRecord().getReference().compare("-") == 0) {

                                        int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                                        std::string ori = "-";
                                        std::string result = sdiRecords.getLast()->getMapSingleRecord().getAlternative() + resultSeq[ai];
                                        Variant mapSingleRecord = Variant(chrName, position, ori, result);
                                        Data* data = new Data(mapSingleRecord);
                                        sdiRecords.deleteLast();
                                        sdiRecords.insertLast(data);
                                    } else {
                                        int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                                        std::string ori = "-";
                                        std::string result(1, resultSeq[ai]);
                                        Variant mapSingleRecord = Variant(chrName, position, ori, result);
                                        Data* data = new Data(mapSingleRecord);
                                        sdiRecords.insertLast(data);
                                    }
//                                std::cout << "1793" << std::endl;
                                }
                            } else {
//                            std::cout << "1796" << std::endl;
                                int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber - 1;
                                std::string ori(1, refSeq[ai]);
                                std::string result(1, resultSeq[ai]);
                                Variant mapSingleRecord = Variant(chrName, position, ori, result);
                                Data* data = new Data(mapSingleRecord);
                                sdiRecords.insertLast(data);
//                            std::cout << "1803" << std::endl;
                            }
                        }
                    }
                }
                {
                    //begin: insert the inter window region SDI record
                    Variant mapSingleRecord;
                    std::string oriSeq;
//                    std::cout << "1811" << std::endl;
                    if ((twoSeqOfMsaResults[i].getRefEnd()) < (twoSeqOfMsaResults[i + 1].getRefStart() - 1)) {
                        int thisSubStart1 = twoSeqOfMsaResults[i].getRefEnd() + 1;
                        int thisSubEnd1 = twoSeqOfMsaResults[i + 1].getRefStart() - 1;
                        oriSeq = getSubsequence(referenceGenome, chrName, thisSubStart1,thisSubEnd1);
                    } else {
                        oriSeq = "-";
                    }
                    std::string resultSeq;
//                    int resultSeqLength = 0;
                    if ( (twoSeqOfMsaResults[i].getResultEnd() >= (twoSeqOfMsaResults[i].getResultStart()))  && (twoSeqOfMsaResults[i].getResultEnd() < (twoSeqOfMsaResults[i + 1].getResultStart() - 1)) ) {
                        int thisSubStart2 = twoSeqOfMsaResults[i].getResultEnd() + 1;
                        int thisSubEnd2 = twoSeqOfMsaResults[i + 1].getResultStart() - 1;
                        resultSeq = getSubsequence(targetGenome, chrName, thisSubStart2, thisSubEnd2);
//                        resultSeqLength = resultSeq.length();
                    } else if( (twoSeqOfMsaResults[i].getResultEnd() < (twoSeqOfMsaResults[i].getResultStart())) &&  (twoSeqOfMsaResults[i].getResultStart() < (twoSeqOfMsaResults[i + 1].getResultStart() - 1)) ){
                        int thisSubStart2 = twoSeqOfMsaResults[i].getResultStart() + 1;
                        int thisSubEnd2 = twoSeqOfMsaResults[i + 1].getResultStart() - 1;

                        resultSeq = getSubsequence(targetGenome, chrName, thisSubStart2,thisSubEnd2);
                    }// if the one is totally deleted
                    else {
                        resultSeq = "-";
                    }
                    if ( (oriSeq.compare(resultSeq)!=0) && (oriSeq.compare("-")==0) ){
                        int position = twoSeqOfMsaResults[i].getRefEnd() + 1;
                        mapSingleRecord = Variant(chrName, position, oriSeq, resultSeq);
//                        std::cout << "2140 " << chrName << " " << position << " " << oriSeq << " " << resultSeq << std::endl;
                        Data* data = new Data (mapSingleRecord);
                        sdiRecords.insertLast(data);
                    }else if (oriSeq.compare(resultSeq) != 0) {
//                        std::cout << "1828" << std::endl;
                        int position = twoSeqOfMsaResults[i].getRefEnd() + 1;
                        mapSingleRecord = Variant(chrName, position, oriSeq, resultSeq);
//                        std::cout << "1981 " << chrName << " " << position << " " << oriSeq << " " << resultSeq << std::endl;
                        Data* data = new Data (mapSingleRecord);
                        sdiRecords.insertLast(data);
                    }
                    i++;
                }
                refLetterNumber =0;
            }
            //end: insert the windows region SDI record
            std::cout << "begin: insert the last windows region SDI record" << std::endl;
            //begin: insert the last windows region SDI record
            if(i == (twoSeqOfMsaResults.size()-1)){
                //System.out.println("i="+i+" twoSeqOfMsaResults.size()-1:"+(twoSeqOfMsaResults.size()-1));
                int refLetterNumber = 0;
                std::string refSeq = twoSeqOfMsaResults[i].getRefSeq();
                std::string resultSeq = twoSeqOfMsaResults[i].getResultSeq();
                for(int ai=0; ai<refSeq.length(); ai++){
                    if (refSeq[ai] != '-') {
                        ++refLetterNumber;
                    }
                    if (refSeq[ai] != resultSeq[ai]) {
                        if (resultSeq[ai] == '-') {
//                            std::cout << "1754" << std::endl;
                            int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber - 1;
                            std::string ori(1, refSeq[ai]);
                            std::string result = "-";
                            Variant mapSingleRecord = Variant(chrName, position, ori, result);
                            Data* data = new Data(mapSingleRecord);
                            sdiRecords.insertLast(data);
//                            std::cout << "1761" << std::endl;
                        } else if (refSeq[ai] == '-') {
//                            std::cout << "1763" << std::endl;
                            if (sdiRecords.getLast() == NULL) {
//                                std::cout << "1755" << std::endl;
                                int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                                std::string ori = "-";
                                std::string result(1, resultSeq[ai]);
                                Variant mapSingleRecord = Variant(chrName, position, ori, result);
                                Data* data = new Data(mapSingleRecord);
                                sdiRecords.insertLast(data);
//                                std::cout << "1772" << std::endl;
                            } else {
//                                std::cout << "1774" << std::endl;
                                if (NULL!=sdiRecords.getLast() && sdiRecords.getLast()->getMapSingleRecord().getPosition()==(twoSeqOfMsaResults[i].getRefStart() + refLetterNumber)
                                    && sdiRecords.getLast()->getMapSingleRecord().getChanginglength() > 0 &&
                                    sdiRecords.getLast()->getMapSingleRecord().getReference().compare("-") == 0) {

                                    int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                                    std::string ori = "-";
                                    std::string result = sdiRecords.getLast()->getMapSingleRecord().getAlternative() + resultSeq[ai];
                                    Variant mapSingleRecord = Variant(chrName, position, ori, result);
                                    Data* data = new Data(mapSingleRecord);
                                    sdiRecords.deleteLast();
                                    sdiRecords.insertLast(data);
                                } else {
                                    int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber;
                                    std::string ori = "-";
                                    std::string result(1, resultSeq[ai]);
                                    Variant mapSingleRecord = Variant(chrName, position, ori, result);
                                    Data* data = new Data(mapSingleRecord);
                                    sdiRecords.insertLast(data);
                                }
//                                std::cout << "1793" << std::endl;
                            }
                        } else {
//                            std::cout << "1796" << std::endl;
                            int position = twoSeqOfMsaResults[i].getRefStart() + refLetterNumber - 1;
                            std::string ori(1, refSeq[ai]);
                            std::string result(1, resultSeq[ai]);
                            Variant mapSingleRecord = Variant(chrName, position, ori, result);
                            Data* data = new Data(mapSingleRecord);
                            sdiRecords.insertLast(data);
//                            std::cout << "1803" << std::endl;
                        }
                    }
                }
            }
            //end: insert the last windows region SDI record
            std::cout << "insert the SDI record after last window" << std::endl;
            //begin: insert the SDI record after last window
            int endIndex = twoSeqOfMsaResults.size()-1;
            std::string oriSeq;
            if((twoSeqOfMsaResults[endIndex].getRefEnd()) < referenceGenome[chrName].getSequence().length()){
                int thisSubStart1 = twoSeqOfMsaResults[endIndex].getRefEnd()+1;
                int thisSubEnd1 = targetGenome[chrName].getSequence().length();
                oriSeq = getSubsequence(referenceGenome, chrName, thisSubStart1, thisSubEnd1);
            }else{
                oriSeq="-";
            }
            if(twoSeqOfMsaResults[endIndex].getResultEnd() < (targetGenome[chrName].getSequence().length()-1)){
                int thisSubStart2 = twoSeqOfMsaResults[endIndex].getResultEnd()+1;
                int thisSubEnd2 = referenceGenome[chrName].getSequence().length();

                std::string resultSeq = getSubsequence(targetGenome, chrName, thisSubStart2, thisSubEnd2);
                if(oriSeq.compare(resultSeq) != 0){
                    int position = twoSeqOfMsaResults[endIndex].getRefEnd()+1;
//                    std::cout << "1917" << std::endl;
                    Variant mapSingleRecord = Variant(chrName, position, oriSeq, resultSeq);
                    Data* data = new Data(mapSingleRecord);
                    sdiRecords.insertLast(data);
                }
            }else if(twoSeqOfMsaResults[endIndex].getResultEnd() == (targetGenome[chrName].getSequence().length()-1)){
                int position = twoSeqOfMsaResults[endIndex].getRefEnd()+1;
                std::string result( 1, targetGenome[chrName].getSequence()[targetGenome[chrName].getSequence().length()-1] );
//                std::cout << "1925" << std::endl;
                Variant mapSingleRecord = Variant(chrName, position, oriSeq, result);
                Data* data = new Data(mapSingleRecord);
                sdiRecords.insertLast(data);
            }else if(twoSeqOfMsaResults[endIndex].getResultEnd() == (targetGenome[chrName].getSequence().length()) && oriSeq.length()>0 && oriSeq.compare("-")!=0 ){
                int position = twoSeqOfMsaResults[endIndex].getRefEnd()+1;
                std::string result = "-";
//                std::cout << "1932" << std::endl;
                Variant mapSingleRecord = Variant(chrName, position, oriSeq, result);
                Data* data = new Data(mapSingleRecord);
                sdiRecords.insertLast(data);
            }else if(oriSeq.length()>0){
                std::cout << "should never run here 2333" << std::endl;
            }
            //end: insert the SDI record after last window

            twoSeqOfMsaResults.clear();// for RAM saving
            //begin: merge link data structure
            std::cout << (*itName) << " link data structure start" << std::endl;

//            if( (sdiRecords[chrName].getFirst() != NULL)
//                && (sdiRecords[chrName].getFirst()->getNext() != NULL ) ){
            for( int runingCound = 0; runingCound<2; ++runingCound){
                if( (sdiRecords.getFirst() != NULL)
                    && (sdiRecords.getFirst()->getNext() != NULL ) ){
                    Data* prevOne = sdiRecords.getFirst();
                    Data* currOne = (sdiRecords.getFirst()->getNext());
//                    std::cout << "2107" << std::endl;
                    while( (currOne != NULL) && (NULL != currOne->getNext())){
//                        std::cout << "2109" << std::endl;
                        if( sdiRecords.getFirst() == currOne ){
                            prevOne = currOne;
                            currOne = prevOne->getNext();
                        }
                        if (currOne->getMapSingleRecord().getChanginglength()<0 && prevOne->getMapSingleRecord().getChanginglength()<0 &&
                            (prevOne->getMapSingleRecord().getPosition()+abs(prevOne->getMapSingleRecord().getChanginglength()))==currOne->getMapSingleRecord().getPosition()
                            && prevOne->getMapSingleRecord().getAlternative().compare("-")==0 && currOne->getMapSingleRecord().getAlternative().compare("-")==0) { // merge to deletions
//                            std::cout << "2117 delete prev" << std::endl;
                            int position = prevOne->getMapSingleRecord().getPosition();
                            std::string ori = prevOne->getMapSingleRecord().getReference() +
                                              currOne->getMapSingleRecord().getReference();
                            std::string result = "-";
                            Variant mapSingleRecord2(chrName, position, ori, result);
                            //delete prev one begin
                            if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                                sdiRecords.deleteFirst();
                            } else {
                                prevOne->getPrev()->setNext(currOne);
                                currOne->setPrev(prevOne->getPrev());
                                delete (prevOne);
                            } //delete prev one end
                            currOne->setMapSingleRecord(mapSingleRecord2);
                            prevOne = currOne->getPrev();
                            if( prevOne == NULL ){
                                prevOne = currOne;
                                currOne = prevOne->getNext();
                            }
//                            std::cout << "2133 delete prev" << std::endl;
                        } else if (currOne->getMapSingleRecord().getChanginglength()==0 && 0 == currOne->getMapSingleRecord().getReference().compare( currOne->getMapSingleRecord().getAlternative())){ // nonsense records
                            //delete current one
                            prevOne->setNext( currOne->getNext() );
                            currOne->getNext()->setPrev(prevOne);
                            delete(currOne);
                            currOne = prevOne->getNext();
//                            std::cout << "2140 delete prev" << std::endl;
                        } else if ( currOne->getMapSingleRecord().getChanginglength()<0 && prevOne->getMapSingleRecord().getChanginglength()>0
                                    && currOne->getMapSingleRecord().getReference().compare(prevOne->getMapSingleRecord().getAlternative())==0 &&
                                    currOne->getMapSingleRecord().getPosition() == prevOne->getMapSingleRecord().getPosition()){ //delete one insertion and next reverse sence deletion
                            //delete current one and prev
                            if( ((currOne->getPrev())) == (sdiRecords.getFirst()) ){
                                sdiRecords.deleteFirst();
                                sdiRecords.deleteFirst();
                            }else{
                                currOne->getPrev()->getPrev()->setNext(currOne->getNext());
                                currOne->getNext()->setPrev(currOne->getPrev()->getPrev());
                                Data *temp = currOne->getNext();
                                delete(currOne->getPrev());
                                delete(currOne);
                                currOne = temp;
                                prevOne = temp->getPrev();
                            }
                        } else if ( currOne->getMapSingleRecord().getChanginglength()>0 && prevOne->getMapSingleRecord().getChanginglength()<0
                                    && currOne->getMapSingleRecord().getAlternative().compare(prevOne->getMapSingleRecord().getReference())==0 &&
                                    currOne->getMapSingleRecord().getPosition() == prevOne->getMapSingleRecord().getPosition()){
                            //delete current one and prev
                            if( ((currOne->getPrev())) == (sdiRecords.getFirst()) ){
                                sdiRecords.deleteFirst();
                                sdiRecords.deleteFirst();
                            }else{
                                currOne->getPrev()->getPrev()->setNext(currOne->getNext());
                                currOne->getNext()->setPrev(currOne->getPrev()->getPrev());
                                Data *temp = currOne->getNext();
                                delete(currOne->getPrev());
                                delete(currOne);
                                currOne = temp;
                                prevOne = temp->getPrev();
                            }
                        } else {
                            prevOne=currOne;
                            currOne=prevOne->getNext();
                        }//std::cout <<  (*itName) << ": link data structure end " << currOne->getMapSingleRecord().getPosition() << std::endl;
                        //std::cout << "2009" << std::endl;
                    }
                }
            }
            //end: merge link data structure

            std::cout << (*itName) << " link data structure end" << std::endl;
            std::vector<Variant> sdiRecordsThisOne;
            if( sdiRecords.getFirst()!=NULL ){
                Data* thisone = sdiRecords.getFirst();
                while( thisone!=NULL ){
                    sdiRecordsThisOne.push_back(thisone->getMapSingleRecord());
                    thisone = (thisone->getNext());
                }
            }

            // clear RAM assigned by new Data() begin
            if( sdiRecords.getFirst()!=NULL ){
                Data* thisone = sdiRecords.getFirst();
                while( thisone!=NULL ){
                    Data *tempData = thisone;
                    thisone = (thisone->getNext());
                    delete(tempData);
                }
            }// clear RAM assigned by new Data() end

            std::cout << " begin to sort" << std::endl;
            // transform link to vector and sort and merge nearby records begin
            bool ifChanged = true;
            while(ifChanged){
                std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
                ifChanged = false;
                std::vector<int> sdiRecordsToRomove;
                int oldSize = sdiRecordsThisOne.size();
                for(int j=1; j<oldSize; j++){
                    if(sdiRecordsThisOne[j].getChanginglength()<0 && sdiRecordsThisOne[j-1].getChanginglength()<0 &&
                       (sdiRecordsThisOne[j-1].getPosition()+abs(sdiRecordsThisOne[j-1].getChanginglength()))==sdiRecordsThisOne[j].getPosition()
                       && sdiRecordsThisOne[j-1].getAlternative().compare("-")==0 && sdiRecordsThisOne[j].getAlternative().compare("-")==0 ){
                        int position = sdiRecordsThisOne[j-1].getPosition();
                        std::string ori = sdiRecordsThisOne[j-1].getReference()+sdiRecordsThisOne[j].getReference();
                        std::string result = "-";
                        Variant mapSingleRecord2(chrName, position, ori, result);
                        sdiRecordsToRomove.push_back(j-1);
                        sdiRecordsThisOne[j] = mapSingleRecord2;
                        j++;
                        ifChanged = true;
                    } else if( sdiRecordsThisOne[j].getReference().compare( sdiRecordsThisOne[j].getAlternative())==0){
                        sdiRecordsToRomove.push_back(j); // it does not affect sorting
                    }
                }
                for( int intTpRomoveIndex = sdiRecordsToRomove.size()-1; intTpRomoveIndex>=0 ;--intTpRomoveIndex ){
                    sdiRecordsThisOne.erase(sdiRecordsThisOne.begin()+sdiRecordsToRomove[intTpRomoveIndex]);
                }
            }
            // transform link to vector and sort and merge nearby records end

            std::ofstream ofile;
            ofile.open("newSdis/" + chrName + "/" + *itName + ".sdi");
            int itVariantNumebt = 0;
            for( std::vector<Variant>::iterator itVariant = sdiRecordsThisOne.begin();
                 itVariant!=sdiRecordsThisOne.end(); ++itVariant  ){
                if( itVariantNumebt >0 ){
                    ofile << std::endl << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" << itVariant->getChanginglength() << "\t" <<
                          itVariant->getReference() << "\t" << itVariant->getAlternative();
                }else{
                    ofile << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" << itVariant->getChanginglength() << "\t" <<
                          itVariant->getReference() << "\t" << itVariant->getAlternative();
                }
                ++itVariantNumebt;
            }
            ofile.close();
        }
    }
}

