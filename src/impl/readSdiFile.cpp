//
// Created by baoxing on 10/10/17.
//

#include "readSdiFile.h"

void readSdiFile(const std::string& filePath, std::map<std::string, std::vector<Variant> >& variantsMap, const std::string& vcfFix, std::map<std::string, Fasta>& referenceGenome){
    readSdiFile(filePath, variantsMap, "", vcfFix, referenceGenome);
}
void readSdiFile (const std::string& filePath, std::map<std::string, std::vector<Variant> >& variantsMap, const std::string& chromosome, const std::string& vcfFix, std::map<std::string, Fasta>& referenceGenome){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening variants file " << filePath << std::endl;
        exit (1);
    }
//    std::regex reg("^(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)");
    std::string line;
//    std::cout << "143" << std::endl;
    while (std::getline(infile, line)){
        if( line.compare(0, 1, "#")==0 ){
            //std::cout << line << std::endl;
            continue;
        }
        if( chromosome.length()>0) {
            if( line.compare(0, chromosome.size(), chromosome)!=0 ){
                continue;
            }
        }
        std::vector<std::string> splits;
        char sep = '\t';
        split(line, sep, splits);
        if( splits.size() >=5 ){
            std::string thisChromosome;
            if( vcfFix.length()>0 ){
                thisChromosome = vcfFix + splits[0];
            }else{
                thisChromosome = splits[0];
            }
            int position = std::stoi(splits[1]);
            std::string reference = splits[3];
            std::string alternative = splits[4];
            transform(reference.begin(), reference.end(), reference.begin(),::toupper);
            transform(alternative.begin(), alternative.end(), alternative.begin(),::toupper);

            if (reference.compare("-") != 0) {
                if( referenceGenome.find(thisChromosome) == referenceGenome.end() ){
                } else if ( referenceGenome[thisChromosome].getSequence().size() < (position+reference.size()-1) ){
                    std::cerr << "the variant position is not in the reference genome range: " << thisChromosome << " " << (position+reference.size()-1) << ". The variant record will be ignored" << std::endl;
                }else if (referenceGenome[thisChromosome].getSequence().substr(position - 1, reference.size()).compare(reference) !=0) {
                    std::cerr << "the record does not confirm with reference genome sequence at: " << thisChromosome << "\t"
                              << position << std::endl<<
                              "reference genome: " << referenceGenome[thisChromosome].getSequence().substr(position - 1, reference.size()) << std::endl <<
                              "sdi: " << reference << ". The variant record will be ignored" << std::endl;
                }else{
                    Variant variant(thisChromosome, position, reference, alternative);
                    if( variantsMap.find(thisChromosome) == variantsMap.end() ){
                        variantsMap[thisChromosome]=std::vector<Variant>();
                    }
                    variantsMap[thisChromosome].push_back(variant);
                }
            }else{
                Variant variant(thisChromosome, position, reference, alternative);
                if( variantsMap.find(thisChromosome) == variantsMap.end() ){
                    variantsMap[thisChromosome]=std::vector<Variant>();
                }
                variantsMap[thisChromosome].push_back(variant);
            }
        }
    }

    for(std::map<std::string, std::vector<Variant> >::iterator it=variantsMap.begin(); it!=variantsMap.end(); ++it) {
        std::sort(it->second.begin(), it->second.end());
    }

    for(std::map<std::string, std::vector<Variant> >::iterator it=variantsMap.begin(); it!=variantsMap.end(); ++it){
        int lastTotalChanged=0;
        for( std::vector<Variant>::iterator it2=(*it).second.begin(); it2!=(*it).second.end(); it2++ ){
            (*it2).setLastTotalChanged(lastTotalChanged);

            int changedPoint= (*it2).getPosition() + lastTotalChanged;
            (*it2).setChangedPoint(changedPoint);
            lastTotalChanged +=(*it2).getChanginglength();
        }
    }
}