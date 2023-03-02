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

