//
// Created by song on 8/4/18.
//

#ifndef ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H
#define ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H
#include "../impl/impl.h"
#include "../myImportandFunction/myImportantFunction.h"
#include <map>
#include <regex>
/*
void TransferAllExonWithSpliceAlignmentResult( const std::string & gffFilePath, const std::string & databaseFastaFilePath,
                                               const std::string & queryFastaFilePath, const std::string & samFile,
                                               std::map<std::string, std::string>& parameters, const std::string & outPutFilePath,
                                               const size_t & minIntron , const bool & slowMode, const int & slidingWindowSize, const size_t & maxLengthForStructureAlignment, int outputTag);
*/

void TransferAllExonWithNucmerResult( const std::string & gffFilePath, const std::string & databaseFastaFilePath,
                                      const std::string & queryFastaFilePath, const std::string & nucmerFilePath,
                                      std::map<std::string, std::string>& parameters, const std::string & outPutFilePath,
                                      const size_t & minIntron , const bool & slowMode, const int & slidingWindowSize,
                                      const size_t & maxLengthForStructureAlignment, const int & alignmentApproach);

void TransferGffWithNucmerResult( const std::string & gffFilePath, const std::string & databaseFastaFilePath,
                                  const std::string & queryFastaFilePath, const std::string & nucmerFilePath, std::map<std::string, std::string>& parameters,
                                  const std::string & outPutFilePath, const size_t & maxLengthForStructureAlignment, const int & alignmentApproach);
void TransferAllExonWithSpliceAlignmentResult( const std::string & gffFilePath, const std::string & databaseFastaFilePath,
                                               const std::string & queryFastaFilePath, const std::string & samFile,
                                               std::map<std::string, std::string>& parameters, const std::string & outPutFilePath,
                                               const size_t & minIntron, const int & slidingWindowSize, const size_t & maxLengthForStructureAlignment, int outputTag,
                                               const int & alignmentApproach);

void TransferAllExonWithNucmerResult(  std::map<std::string, std::vector<std::string> > & geneNameMap,
                                       std::map<std::string, Gene> & geneHashMap, std::map<std::string, Transcript> &transcriptHashMap,
                                       const std::string & databaseFastaFilePath,
                                       const std::string & queryFastaFilePath, std::map<std::string, std::vector<AlignmentMatch>> & alignmentMatchsMap,
                                       std::map<std::string, std::string>& parameters, const std::string & outPutFilePath,
                                       const size_t & minIntron , const bool & slowMode, const int & slidingWindowSize, const size_t & maxLengthForStructureAlignment,
                                       const std::string source, int outputTag, const int & alignmentApproach);

#endif //ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H
