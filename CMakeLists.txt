cmake_minimum_required(VERSION 4.3)
project(GEAN)

set(CMAKE_CXX_STANDARD 14)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
else()
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
endif()



set(CMAKE_CXX_FLAGS "-Wall -g3 -O3 -funroll-all-loops -I ./ -std=gnu++14 -pthread")
set(CMAKE_C_FLAGS "-Wall -g3 -O3 -funroll-all-loops -I ./ -pthread")



add_executable(gean
        src/myImportandFunction/alignNeedlemanForTranscript.cpp
        src/myImportandFunction/alignNeedlemanForTranscript.h
        src/myImportandFunction/gotohAlignment.cpp
        src/myImportandFunction/gotohAlignment.h
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
        src/model/Matrix.cpp
        src/model/Matrix.h
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
        src/service/MSA/reSdiFromMsa.cpp
        src/service/MSA/reSdiFromMsa.h
        src/impl/getSequencesFromGff.cpp
        src/impl/getSequencesFromGff.h
        src/control/myCoordinateLiftOver.cpp
        src/control/myCoordinateLiftOver.h
        src/control/myReAnnotationAndExonerateAndNovo.cpp
        src/control/myReAnnotationAndExonerateAndNovo.h
        src/control/generateRandomSdi.cpp
        src/control/generateRandomSdi.h
        src/control/myControl.h
        src/service/service.h
        src/myImportandFunction/myImportantFunction.h
        src/impl/TranscriptUpdateInformation.cpp
        src/impl/TranscriptUpdateInformation.h
        src/impl/FileWriteCatch.cpp
        src/impl/FileWriteCatch.h
        src/myImportandFunction/alignNeedlemanForTranscript_simd_avx2int32.cpp
        src/myImportandFunction/alignNeedlemanForTranscript_simd_avx2int32.h
        src/model/Range.cpp src/model/Range.h src/model/AlignmentMatch.cpp
        src/model/AlignmentMatch.h src/impl/NucmerRead.cpp
        src/impl/NucmerRead.h
        src/service/TransferGffWithNucmerResult.cpp
        src/service/TransferGffWithNucmerResult.h
        src/myImportandFunction/alignSlidingWindow.cpp
        src/myImportandFunction/alignSlidingWindow.h
        src/impl/readGffFileWithEverything.cpp
        src/impl/readGffFileWithEverything.h
        src/impl/CheckAndUpdateTranscriptsEnds.cpp
        src/impl/CheckAndUpdateTranscriptsEnds.h
        src/service/UpdateVariants.cpp
        src/service/UpdateVariants.h
        src/myImportandFunction/alignNeedlemanForTranscript_simd_avx2int16.cpp
        src/myImportandFunction/alignNeedlemanForTranscript_simd_avx2int16.h
        src/myImportandFunction/alignTranscript.cpp
        src/myImportandFunction/alignTranscript.h
        src/myImportandFunction/geneAnnotationAlignment.cpp
        src/myImportandFunction/geneAnnotationAlignment.h
        src/service/denovoGffLiftOverToUniqGffOutput.cpp
        src/service/denovoGffLiftOverToUniqGffOutput.h
        src/impl/organizeGffRecords.cpp
        src/impl/organizeGffRecords.h
        src/impl/deNovoGenomeVariantCalling.cpp
        src/impl/deNovoGenomeVariantCalling.h
        src/model/OrthologPair.cpp
        src/model/OrthologPair.h
        src/impl/geneSyntenic.cpp
        src/impl/geneSyntenic.h
        src/impl/outPutORFConserveredTranscripts.cpp
        src/impl/outPutORFConserveredTranscripts.h)
