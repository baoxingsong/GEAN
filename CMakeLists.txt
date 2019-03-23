cmake_minimum_required(VERSION 3.0)
project(GEAN)
#=SET (CMAKE_C_COMPILER "/usr/bin/clang")
#SET (CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
set(CMAKE_CXX_FLAGS "-Wall -msse2 -g3 -O3 -funroll-all-loops -mavx2 -I ./ -std=c++11 -pthread")
set(CMAKE_C_FLAGS "-Wall -msse2 -g3 -O3 -funroll-all-loops -mavx2 -I ./ -pthread")

#for google unit test
add_subdirectory(./googletest/googletest)
include_directories(googletest/googletest/include googletest/googletest)


add_executable(gean
        googletest
        googletest/googletest
        src/myImportandFunction/alignNeedlemanForTranscript.cpp
#       ksw2_extz2_sse.c
#       ksw2.h
#       src/myImportandFunction/alignNeedlemanForTranscript_simd.cpp
        src/myImportandFunction/alignNeedlemanForTranscript.h
        src/myImportandFunction/alignNeedlemanWunsch_v2.cpp
        src/myImportandFunction/alignNeedlemanWunsch.h
        src/model/GenomeBasicFeature.cpp
        src/model/GenomeBasicFeature.h
        src/controlLayer.cpp
        src/controlLayer.h
        src/InputParser.cpp
        src/InputParser.h
        src/model/model.h
        src/util/myutil.cpp
        src/util/myutil.h
        src/util/nucleotideCodeSubstitutionMatrix.cpp
        src/util/nucleotideCodeSubstitutionMatrix.h
        src/util/parameters.cpp
        src/util/parameters.h
        src/sole.h
        song.cpp
        src/model/Fasta.cpp
        src/model/Fasta.h
        src/model/variant.cpp
        src/model/variant.h
        src/model/Transcript.cpp
        src/model/Transcript.h
        src/model/Gene.cpp
        src/model/Gene.h
        src/model/MsaWindow.cpp
        src/model/MsaWindow.h
        src/model/MsaSingleRecord.cpp
        src/model/MsaSingleRecord.h
        src/model/MsaFileRecord.cpp
        src/model/MsaFileRecord.h
        src/model/TwoSeqOfMsaResult.cpp
        src/model/TwoSeqOfMsaResult.h
        src/model/FirstLastList.cpp
        src/model/FirstLastList.h
        src/model/STRAND.h
        src/model/REGION.h
        src/impl/impl.h
        src/impl/GetPseudoGenomeSequence.cpp
        src/impl/GetPseudoGenomeSequence.h
        src/impl/WriteFasta.cpp
        src/impl/WriteFasta.h
        src/impl/GetReverseComplementary.cpp
        src/impl/GetReverseComplementary.h
        src/impl/readFastaFile.cpp
        src/impl/readFastaFile.h
        src/impl/annotationLiftOver.cpp
        src/impl/annotationLiftOver.h
        src/impl/checkOrfState.cpp
        src/impl/checkOrfState.h
        src/util/util.h
        src/impl/nA2AA.cpp
        src/impl/nA2AA.h
        src/impl/getSubsequence.cpp
        src/impl/getSubsequence.h
        src/impl/readSdiFile.cpp
        src/impl/readSdiFile.h
        src/impl/readGffFile.cpp
        src/impl/readGffFile.h
        src/impl/coordinateLiftOver.cpp
        src/impl/coordinateLiftOver.h
        src/service/reannotation/TranscriptsTogenes.cpp
        src/service/reannotation/TranscriptsTogenes.h
        src/service/reannotation/reAnnotationAndExonerateAndNovo.cpp
        src/service/reannotation/reAnnotationAndExonerateAndNovo.h
        src/service/reannotation/transcriptRealignmentAndExonerate.cpp
        src/service/reannotation/transcriptRealignmentAndExonerate.h
        src/service/reannotation/runExonerate.cpp
        src/service/reannotation/runExonerate.h
        src/service/MSA/cutWindow.cpp
        src/service/MSA/cutWindow.h
        src/service/MSA/transcriptRealignmentPversion.cpp
        src/service/MSA/transcriptRealignmentPversion.h
        src/service/MSA/reSdiFromMsa.cpp
        src/service/MSA/reSdiFromMsa.h
        src/impl/getSequencesFromGff.cpp
        src/impl/getSequencesFromGff.h
        src/control/myCoordinateLiftOver.cpp
        src/control/myCoordinateLiftOver.h
        src/control/myReAnnotationAndExonerateAndNovo.cpp
        src/control/myReAnnotationAndExonerateAndNovo.h
        src/control/VariantsCountAsPhenotypeForAssociation.cpp
        src/control/VariantsCountAsPhenotypeForAssociation.h
        src/control/myControl.h
        src/service/service.h
        src/cns.cpp
        src/myImportandFunction/myImportantFunction.h
        src/impl/TranscriptUpdateInformation.cpp
        src/impl/TranscriptUpdateInformation.h
        src/impl/FileWriteCatch.cpp
        src/impl/FileWriteCatch.h
        src/myImportandFunction/alignNeedlemanWunsch_simd.cpp
        src/myImportandFunction/alignNeedlemanWunsch_simd.h
        src/tests/util/nucleotideCodeSubstitutionMatrix_test.cpp
        src/tests/myImportandFunction/alignNeedlemanWunsch_simd_test.cpp
        src/myImportandFunction/alignNeedlemanForTranscript_simd_avx2int32.cpp
        src/myImportandFunction/alignNeedlemanForTranscript_simd_avx2int32.h
        src/tests/control/myReAnnotationAndExonerateAndNovo_test.cpp
        src/tests/myImportandFunction/alignNeedlemanForTranscript_simd_avx2int32.cpp
        src/model/Range.cpp src/model/Range.h src/model/AlignmentMatch.cpp
        src/model/AlignmentMatch.h src/impl/NucmerRead.cpp
        src/impl/NucmerRead.h
        src/service/TransferGffWithNucmerResult.cpp
        src/service/TransferGffWithNucmerResult.h
        src/tests/service/TransferGffWithNucmerResult_test.cpp
        src/myImportandFunction/alignSlidingWindow.cpp
        src/myImportandFunction/alignSlidingWindow.h
        src/tests/myImportandFunction/alignSlidingWindow_test.cpp
        src/impl/readGffFileWithEverything.cpp
        src/impl/readGffFileWithEverything.h
        src/tests/impl/readGffFileWithEverything.cpp
        src/tests/std/bitsetTest.cpp
        src/impl/CheckAndUpdateTranscriptsEnds.cpp
        src/impl/CheckAndUpdateTranscriptsEnds.h
        src/service/UpdateVariants.cpp
        src/service/UpdateVariants.h
        src/tests/service/UpdateVariants_test.cpp
        src/myImportandFunction/alignNeedlemanForTranscript_simd_avx2int16.cpp
        src/myImportandFunction/alignNeedlemanForTranscript_simd_avx2int16.h
        src/myImportandFunction/alignTranscript.cpp
        src/myImportandFunction/alignTranscript.h
        src/tests/myImportandFunction/alignTranscript.cpp
        src/tests/control/myGffCoordinateLiftOver_test.cpp
        src/tests/impl/readFastaFile_test.cpp
        src/tests/impl/deNovoGenomeVariantCalling_test.cpp
        src/tests/impl/longest_test.cpp
        src/tests/myImportandFunction/hengAlign_test.cpp
        src/tests/impl/orf_test.cpp
        src/myImportandFunction/geneAnnotationAlignment.cpp
        src/myImportandFunction/geneAnnotationAlignment.h
        src/service/denovoGffLiftOverToUniqGffOutput.cpp
        src/service/denovoGffLiftOverToUniqGffOutput.h
        src/impl/organizeGffRecords.cpp
        src/impl/organizeGffRecords.h
        src/impl/deNovoGenomeVariantCalling.cpp
        src/impl/deNovoGenomeVariantCalling.h
        src/service/deNovoGenomeVariantCallingService.cpp
        src/service/deNovoGenomeVariantCallingService.h
        src/myImportandFunction/hengAlign.cpp
        src/myImportandFunction/hengAlign.h
        src/cns.cpp src/song_cns.cpp
        src/song_cns.h
        src/tests/service/TransferAllExonWithSpliceAlignmentResult_test.cpp
        src/model/OrthologPair.cpp
        src/model/OrthologPair.h
        src/impl/geneSyntenic.cpp
        src/impl/geneSyntenic.h
        src/impl/outPutORFConserveredTranscripts.cpp
        src/impl/outPutORFConserveredTranscripts.h)

target_link_libraries(gean gtest gtest_main)
