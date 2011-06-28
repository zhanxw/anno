#ifndef _GENOMESEQUENCE_H_
#define _GENOMESEQUENCE_H_

class GenomeSequence{
public:
    /**
     * @return true: if loads successful
     */
    bool open(const char* fileName){
        LineReader lr(fileName);
        std::string line;
        std::string chr = "";
        while(lr.readLine(&line) > 0) {
            if (line.size() > 0) {
                if (line[0] == '>') {
                    // new chromosome
                    unsigned int idx = line.find_first_of(' ');
                    chr = line.substr(1, idx - 1);
                    this->data[chr] = "";
                } else {
                    stringStrip(&line);
                    if (this->data.find(chr) != this->data.end())
                        this->data[chr] += line;
                }
            }
        }
        return true;
    };
    /**
     * @return total number of chromosome
     */
    int size() const {
        return this->data.size();
    };
    /**
     * @return total number of chromosome
     */
    long int getGenomeLength() const {
        long int l = 0;
        std::map<std::string, std::string>::const_iterator it;
        for (it = this->data.begin(); it != this->data.end() ; it++) {
            l += (it->second).size();
        }
        return l;
    };
    std::string& getChromosome(const char* c){
        return this->data[c];
    };
    const std::string& operator[] (const std::string& c) {
        return this->data[c];
    }
    bool exists(const std::string& c){
        if (this->data.find(c) != this->data.end())
         return true;
        return false;
    }
public:
    std::map<std::string, std::string> data;
};



#endif /* _GENOMESEQUENCE_H_ */
