#ifndef _STRINGTEMPLATE_H_
#define _STRINGTEMPLATE_H_

/**
 * Specify a template format
 * e.g.
 * int n;
 * std::string result;
 * StringTemplate s("Hello $(NAME).");
 * s.add("NAME", "world");
 * s.translate(&result);
 * printf("%s\n", result.c_str());
 *
 * Output:
 * Hello world.
 */
class StringTemplate{
public:
    StringTemplate(const char* t){
        this->temp = t;
    };
    StringTemplate(const std::string& t){
        this->temp = t;
    };
    void add(const char* key, const char* value){
        std::string k = key;
        this->data[k] = value;
    };
    /**
     *@param s: translation results is s
     *@return: how many key words are translateds
     */
    int translate(std::string* s){
        assert(s);
        s->assign(this->temp);
        int nTrans = 0; // count how many keywords are translated
        std::map<std::string, std::string>::const_iterator it;
        for (it = this->data.begin(); it != this->data.end(); it++){
            const std::string& key = it->first;
            const std::string& value = it->second;
            int pos = s->find(key);
            int kLen= key.size();
            // check if the keyword in the template looks like $(key)
            if (pos == std::string::npos || 
                pos - 2 < 0 || (*s)[pos-2] != '$' ||
                (*s)[pos-1] != '(' ||
                pos + kLen >= s->size() || (*s)[pos+kLen] != ')') 
                continue;
            nTrans++;
            s->replace(pos - 2 , kLen + 3, value);
        }
        return nTrans;
    };
private:
    std::map<std::string, std::string> data;
    std::string temp; // temp meaning template, as template is a keyword.
};

#endif /* _STRINGTEMPLATE_H_ */
