//
// Created by baoxing on 10/10/17.
//

#include "readGffFile.h"

#include <array>
#include <cctype>
#include <set>
#include <stdexcept>

namespace {

std::string trim(const std::string & value){
    size_t begin = 0;
    while( begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin])) ){
        ++begin;
    }
    size_t end = value.size();
    while( end > begin && std::isspace(static_cast<unsigned char>(value[end-1])) ){
        --end;
    }
    return value.substr(begin, end-begin);
}

std::string lowerCase(const std::string & value){
    std::string result = value;
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char character){
        return static_cast<char>(std::tolower(character));
    });
    return result;
}

bool splitGffColumns(const std::string & line, std::array<std::string, 9> & columns){
    size_t begin = 0;
    for( size_t column = 0; column < 8; ++column ){
        size_t tab = line.find('\t', begin);
        if( tab == std::string::npos ){
            return false;
        }
        columns[column] = line.substr(begin, tab-begin);
        begin = tab+1;
    }
    columns[8] = line.substr(begin);
    return true;
}

std::vector<std::string> splitAttributes(const std::string & attributes){
    std::vector<std::string> entries;
    size_t begin = 0;
    char quote = '\0';
    for( size_t index = 0; index < attributes.size(); ++index ){
        char character = attributes[index];
        if( quote == '\0' && (character == '\"' || character == '\'') ){
            quote = character;
        }else if( quote != '\0' && character == quote ){
            quote = '\0';
        }else if( quote == '\0' && character == ';' ){
            entries.push_back(attributes.substr(begin, index-begin));
            begin = index+1;
        }
    }
    entries.push_back(attributes.substr(begin));
    return entries;
}

void appendParentIds(const std::string & rawValue, std::vector<std::string> & parentIds,
                     std::set<std::string> & observedParentIds){
    std::string value = trim(rawValue);
    if( value.size() >= 2 &&
        ((value.front() == '\"' && value.back() == '\"') ||
         (value.front() == '\'' && value.back() == '\'')) ){
        value = value.substr(1, value.size()-2);
    }

    size_t begin = 0;
    while( begin <= value.size() ){
        size_t comma = value.find(',', begin);
        std::string parentId = trim(value.substr(begin, comma == std::string::npos ? std::string::npos : comma-begin));
        if( !parentId.empty() && observedParentIds.insert(parentId).second ){
            parentIds.push_back(parentId);
        }
        if( comma == std::string::npos ){
            break;
        }
        begin = comma+1;
    }
}

std::vector<std::string> getParentIds(const std::string & attributes, const std::regex & customParentRegex){
    std::vector<std::string> parentIds;
    std::set<std::string> observedParentIds;

    for( const std::string & rawEntry : splitAttributes(attributes) ){
        std::string entry = trim(rawEntry);
        if( entry.empty() ){
            continue;
        }

        size_t separator = entry.find('=');
        std::string key;
        std::string value;
        if( separator != std::string::npos ){
            key = trim(entry.substr(0, separator));
            value = entry.substr(separator+1);
        }else{
            separator = entry.find_first_of(" \t");
            if( separator == std::string::npos ){
                continue;
            }
            key = trim(entry.substr(0, separator));
            value = entry.substr(separator+1);
        }

        std::string normalizedKey = lowerCase(key);
        if( normalizedKey == "parent" || normalizedKey == "transcript_id" || normalizedKey == "transcriptid" ){
            appendParentIds(value, parentIds, observedParentIds);
        }
    }

    if( !parentIds.empty() ){
        return parentIds;
    }

    std::smatch match;
    if( std::regex_search(attributes, match, customParentRegex) && match.size() > 1 ){
        for( size_t group = match.size()-1; group > 0; --group ){
            if( match[group].matched && !match[group].str().empty() ){
                appendParentIds(match[group].str(), parentIds, observedParentIds);
                break;
            }
        }
    }
    return parentIds;
}

}

void readGffFile (const std::string& filePath, std::map<std::string, std::vector<Transcript> >& transcriptHashSet){
    std::string cdsParentRegex="([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$";
    readGffFile (filePath, transcriptHashSet, cdsParentRegex);
}

void get_transcript_to_gene_map_from_gff (const std::string& filePath, std::map<std::string, std::string >& transcript_to_gene_map, const std::string& regex){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit (1);
    }
    std::regex reg(regex);
    std::string line;
    while (std::getline(infile, line)){
        if( line.empty() ){
            continue;
        }
        std::smatch match;
        regex_search(line, match, reg);
        if( match.empty() || line[0]=='#'){
        }else{
            std::string transcript_id = match[1];
            std::string gene_id = match[2];
            transcript_to_gene_map[transcript_id] = gene_id;
        }
    }
}

void readGffFile (const std::string& filePath, std::map<std::string, std::vector<Transcript> >& transcriptHashSet, const std::string& cdsParentRegex){
    std::map<std::string, Transcript> transcriptHashMap;
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit (1);
    }
    std::regex customParentRegex(cdsParentRegex);
    std::string line;
    size_t lineNumber = 0;
    while (std::getline(infile, line)){
        ++lineNumber;
        if( line.empty() || line[0]=='#' ){
            continue;
        }

        std::array<std::string, 9> columns;
        if( !splitGffColumns(line, columns) || columns[2] != "CDS" ){
            continue;
        }

        int start;
        int end;
        try{
            size_t parsedStart = 0;
            size_t parsedEnd = 0;
            start = std::stoi(columns[3], &parsedStart);
            end = std::stoi(columns[4], &parsedEnd);
            if( parsedStart != columns[3].size() || parsedEnd != columns[4].size() ){
                throw std::invalid_argument("invalid CDS coordinate");
            }
        }catch( const std::exception & ){
            std::cerr << "invalid CDS coordinates at " << filePath << ":" << lineNumber << std::endl;
            continue;
        }
        if(start>end){
            std::swap(start, end);
        }

        std::vector<std::string> parentIds = getParentIds(columns[8], customParentRegex);
        if( parentIds.empty() ){
            std::cerr << "could not find a CDS parent at " << filePath << ":" << lineNumber << std::endl;
            continue;
        }

        STRAND strand = columns[6] == "-" ? NEGATIVE : POSITIVE;
        for( const std::string & parentId : parentIds ){
            if( transcriptHashMap.find(parentId) == transcriptHashMap.end() ){
                Transcript transcript(parentId, columns[0], strand);
                transcript.setSource(columns[1]);
                transcriptHashMap[parentId] = transcript;
            }
            GenomeBasicFeature cds(start, end, columns[5], columns[7], columns[8]);
            cds.setType(columns[2]);
            transcriptHashMap[parentId].addCds(cds);
        }
    }
    for (std::map<std::string, Transcript>::iterator it=transcriptHashMap.begin(); it!=transcriptHashMap.end(); ++it){
        if(transcriptHashSet.find(it->second.getChromeSomeName()) == transcriptHashSet.end() ){
            transcriptHashSet[it->second.getChromeSomeName()]=std::vector<Transcript>();
        }
        it->second.updateInforCds();
        transcriptHashSet[it->second.getChromeSomeName()].push_back(it->second);
    }
    for (std::map<std::string, std::vector<Transcript>>::iterator it=transcriptHashSet.begin();
            it!=transcriptHashSet.end(); ++it){
        std::sort(it->second.begin(), it->second.end(), [](Transcript a, Transcript b) {
            return a.getPStart() < b.getPStart();
        });
    }
}
