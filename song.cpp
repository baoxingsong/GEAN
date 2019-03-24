/*
 * =====================================================================================
 *
 *       Filename:  song.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/23/2017 21:51:57
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

#include "src/controlLayer.h"
#include "./googletest/googletest/include/gtest/gtest.h"


using namespace std;


int main(int argc, char** argv){
//testing::InitGoogleTest(&argc, argv);
//RUN_ALL_TESTS();
//return 0;
    if( argc<=1 ){
        usage();
        return 1;
    }
    std::string program = argv[1];
    if( program.compare("-h") == 0 || program.compare("--help") == 0 ){
        usage();
        exit(1);
    }

    InputParser inputParser (argc, argv);
    string parameterFile;
    std::string exepath = getexepath(argv);
    if( inputParser.cmdOptionExists("-parameter")){
        parameterFile = inputParser.getCmdOption("-parameter");
    }else{
        parameterFile = exepath + "/configure";
    }
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, exepath);
    if( program.compare("pseudogeno") == 0 ) {
        return getGenomeSequence(--argc, ++argv, parameters);
    }else if( program.compare("lift") == 0 ) {
        return getChangedFromReference(--argc, ++argv, parameters);
    }else if( program.compare("revlift") == 0 ) {
        return getReferenceFromChanged(--argc, ++argv, parameters);
    }else if( program.compare("liftgff") == 0 ) {
        return gffCoordinateLiftOver(--argc, ++argv, parameters);
    }else if( program.compare("revliftgff") == 0 ) {
        return revGffCoordinateLiftOver(--argc, ++argv, parameters);
    }else if( program.compare("gff2seq") == 0 ) {
        return getSequences(--argc, ++argv, parameters);
    }else if( program.compare("premsa") == 0 ) {
        return myPrepareForMsa(--argc, ++argv, parameters);
    }else if( program.compare("msatosdi") == 0 ){
        return myConstructSdiFromMsa(--argc, ++argv, parameters);
    }else if( program.compare("annowgr") == 0 ) {
        return reAnnotationAndExonerateAndNovo(--argc, ++argv, parameters);
    }else if( program.compare("transgff") == 0 ) {
        return TransGff(--argc, ++argv, parameters);
    }else if( program.compare("spltogff") == 0 ) {
        return spliceAlignmentToGff(--argc, ++argv, parameters);
    } else if( program.compare("purifygff") == 0 ) {
        return PurifyGff(--argc, ++argv, parameters);
    } else if( program.compare("sinsyn") == 0 ) {
        return syntenicSingleCopy(--argc, ++argv, parameters);
    } else if( program.compare("sinsyn2") == 0 ) {
        return syntenicSingleCopy2(--argc, ++argv, parameters);
    } else if( program.compare("quotasyn") == 0 ) {
        return syntenicSingleCopy3(--argc, ++argv, parameters);
    } else if( program.compare("orf") == 0 ) {
        return outPutORFConserveredTranscripts(--argc, ++argv, parameters);
    } else if( program.compare("reanva") == 0 ) {
        return Reanva(--argc, ++argv, parameters);
    } else if( program.compare("varcall") == 0 ) {
        return DenoveAssemblyVariantCalling(--argc, ++argv, parameters);
    }
    // The following functions are for testing purpose
    else if( program.compare("countNumberOfTwoneighborSNP") == 0 ) {
        return myCountNumberOfTwoneighborSNP(--argc, ++argv, parameters);
    }else if( program.compare("countNumberSNPAndIndel") == 0 ) {
        return mycountNumberSNPAndIndel(--argc, ++argv, parameters);
    }else if ( program.compare("randomVar") == 0 ){
        return myGenerateRandomSdi(--argc, ++argv, parameters);
    }else{
        usage();
    }
    return 0;
}
