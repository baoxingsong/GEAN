//
// Created by song on 8/4/18.
//

#ifndef ANNOTATIONLIFTOVER_NUCMERREAD_H
#define ANNOTATIONLIFTOVER_NUCMERREAD_H

#include "../model/model.h"
#include "../util/util.h"

void nucmerRead(const std::string & filePath, std::map<std::string, std::vector<AlignmentMatch>> & alignmentMatchsMap);


#endif //ANNOTATIONLIFTOVER_NUCMERREAD_H
