/*
 * =====================================================================================
 *
 *       Filename:  controlLayer.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 09:38:17
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/

#ifndef _CONTROLLAYER_H
#define _CONTROLLAYER_H

#include <iostream>
#include "InputParser.h"
#include <sstream>
#include "./util/util.h"
#include "./model/model.h"
#include "./service/service.h"
#include "./control/myControl.h"
#include "./impl/impl.h"
#include "./myImportandFunction/myImportantFunction.h"

int getGenomeSequence(int argc, char** argv, std::map<std::string, std::string>& parameters);


int getChangedFromReference( int argc, char** argv, std::map<std::string, std::string>& parameters );
int getReferenceFromChanged( int argc, char** argv, std::map<std::string, std::string>& parameters );
int gffCoordinateLiftOver( int argc, char** argv, std::map<std::string, std::string>& parameters );
int revGffCoordinateLiftOver( int argc, char** argv, std::map<std::string, std::string>& parameters );
//int gffExonLiftOverCheckIntronSize( int argc, char** argv, std::map<std::string, std::string>& parameters );
int Reanva( int argc, char** argv, std::map<std::string, std::string>& parameters );
int getSequences(int argc, char** argv, std::map<std::string, std::string>& parameters);
int TransGff( int argc, char** argv, std::map<std::string, std::string>& parameters );
int spliceAlignmentToGff( int argc, char** argv, std::map<std::string, std::string>& parameters );
int PurifyGff( int argc, char** argv, std::map<std::string, std::string>& parameters );
int DenoveAssemblyVariantCalling( int argc, char** argv, std::map<std::string, std::string>& parameters );
int annotationLiftOver( int argc, char** argv );
int annotationLiftOverAndOrth( int argc, char** argv);
int reAnnotationAndExonerateAndNovo( int argc, char** argv, std::map<std::string, std::string>& parameters);
int myCountNumberOfTwoneighborSNP( int argc, char** argv, std::map<std::string, std::string>& parameters);
int mycountNumberSNPAndIndel( int argc, char** argv, std::map<std::string, std::string>& parameters);
int myGenerateRandomSdi( int argc, char** argv, std::map<std::string, std::string>& parameters);
int myPrepareForMsa( int argc, char** argv, std::map<std::string, std::string>& parameters);
int myConstructSdiFromMsa( int argc, char** argv, std::map<std::string, std::string>& parameters);
#endif
