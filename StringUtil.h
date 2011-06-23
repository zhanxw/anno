#ifndef _STRINGUTIL_H_
#define _STRINGUTIL_H_

/** tokenize the string
    @return number of tokens we obtained
    e.g. For empty input string, we will return 1, and result will have
only 1 element (the empty string)
*/
int stringTokenize(const std::string& str, const std::string& delim, std::vector<std::string>* result){
    assert(result);
    result->clear();
    std::string s;
    unsigned int l = str.size();
    unsigned int i = 0;
    while (i < l) {
        if (delim.find(str[i]) != std::string::npos) { // it's a delimeter
            result->push_back(s);
            s.clear();
        } else {
            s.push_back(str[i]);
        }
        ++i;
    };
    result->push_back(s);    
    return result->size();
};
int stringTokenize(const std::string& str, const char delim, std::vector<std::string>* result){
    std::string d;
    d.push_back(delim);
    return (stringTokenize(str, d, result));
};

//remove leading and trailing characters
void stringStrip(std::string* input, const char* characters = " ") {
    size_t beg = input->find_first_not_of(characters);
    size_t end = input->find_last_not_of(characters);
    input->assign( input->substr(beg, end-beg) );
};

//extract piece of string from @param beg(inclusive) to @param end(exclusive)
//NOTE: beg and end can be negative meaning count from right hand side. 
void stringSlice(std::string* input, int beg, int end) {
    assert(input);
    unsigned int len = input->size();
    if (beg < 0) beg += len;
    if (end < 0) end += len;
    assert (beg >= 0 && end >= 0);
    input -> assign ( input->substr(beg, end- beg)) ;
};

std::string stringJoin(const std::vector<std::string>& input, const std::string& delim) {
    std::string s;
    if (input.size() == 0) {
        return s;
    }
    s = input[0];
    for (unsigned int i = 1; i < input.size(); i++) {
        s+= delim;
        s+= input[i];
    }
    return s;
};

std::string toUpper(const std::string& s) {
    std::string r;
    for (unsigned int i = 0; i < s.size(); i++) {
        r.push_back(toupper(s[i]));
    }
    return r;
};

std::string toLower(const std::string& s) {
    std::string r;
    for (unsigned int i = 0; i < s.size(); i++) {
        r.push_back(toupper(s[i]));
    }
    return r;
};

#endif /* _STRINGUTIL_H_ */
