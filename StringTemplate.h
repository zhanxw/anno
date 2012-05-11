#ifndef _STRINGTEMPLATE_H_
#define _STRINGTEMPLATE_H_

#include <string>
#include <vector>
#include <map>
#include <stdio.h>

/**
 * For ANNO, support keywords are:
 * TOP_GENE_NAME1
 * TOP_GENE_NAME
 * TOP_TYPE
 * GENE_NAME       can be vector
 * GENE_STRAND     ..
 * TYPE
 * TYPE1
 */

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
  struct KEY;
  struct VALUE;
  class Array {
 public:
    int translate(std::string* str, const std::map<std::string, VALUE>& dict);
    /**
     * Parse s[beg...end]  beg: inclusive end: exclusive
     */
    int parse(const char* s, int param_beg, int param_end);
    void clear() {
      this->data.clear();
      this->dict.clear();
      this->delim.clear();
    };
 private:
    std::vector<KEY> data;
    std::map<std::string, VALUE> dict; // store key->value pair
    std::string delim;
  };

  typedef enum {
    UNDEFINED_KEY = 0,
    TEXT = 1,
    KEYWORD = 2,
    ARRAY = 3            // will be expanded recursively
  } KEY_TYPE;

  typedef enum {
    UNDEFINED_TYPE = 0,
    STRING = 1,
    STRING_ARRAY = 2
  } VALUE_TYPE;

  struct VALUE {
    VALUE_TYPE type;
    std::string string;
    std::vector<std::string> array;
  } ;

  struct KEY{
    KEY_TYPE type;
    std::string text;
    std::string keyword;
    Array array;
    void clear() {
      type = UNDEFINED_KEY;
      text.clear();
      keyword.clear();
      array.clear();
    }

  };

  explicit StringTemplate() {};
  explicit StringTemplate(const char* s) {
    this->parse(s);
  };
  explicit StringTemplate(const std::string& s) {
    this->parse(s);
  };
private: // forbid two default copy-constructor
  StringTemplate(const StringTemplate& s) {
    if (this != &s){
      this->data = s.data;
      this->dict = s.dict;
    }
  };
  StringTemplate& operator= (const StringTemplate& s) {
    this->clear();
    this->data = s.data;
    this->dict = s.dict;
  };
public:
  StringTemplate& operator= (const char* s) {
    this->clear();
    this->parse(s);
    return *this;
  };
  StringTemplate& operator= (const std::string& s) {
    this->clear();
    this->parse(s);
    return *this;
  };

  std::vector<KEY> data;
  std::map<std::string, VALUE> dict; // store key->value pair

  static bool isValidKeyword(const char* s){
    while ( *s != '\0') {
      if (!isalnum(*s))
        return false;
      s++;
    }
    return true;
  }
  void add(const char* key, const char* value) {
    VALUE v;
    v.type = STRING;
    v.string = value;
    this->dict[key] = v;
  };
  void add(const char* key, const std::string& value) {
    this->add(key, value.c_str());
  };
  void add(const char* key, const std::vector<std::string>& value) {
    VALUE v;
    v.type = STRING_ARRAY;
    v.array = value;
    this->dict[key] = v;
  };

  /**
   * @param str: store translated result in str
   * @return 0: if success
   */
  int translate(std::string* str) {
    std::string& s = *str;
    s.clear();
    for (unsigned int i = 0; i < data.size(); i++) {
      KEY& k = data[i];
      switch(k.type) {
        case TEXT:
          s += k.text;
          break;
        case KEYWORD:
          if (dict.find(k.keyword) == dict.end()) {
            fprintf(stderr, "Data translation error for key %s!\n", k.keyword.c_str());
          } else {
            VALUE& v = dict[k.keyword];
            if (v.type == STRING) {
              s += v.string;
            } else {
              if (v.array.size() > 0) {
                s += v.array[0];
              }
            };
          };
          break;
        case ARRAY:
          int ret = k.array.translate(str, this->dict);
          if (ret) {
            return ret;
          };
          break;
      };
    }
    return 0;
  };
  /**
   *  @return 0: when parse succeeds
   */
  int parse(const std::string& s) {
    return this->parse(s.c_str());
  };
  int parse(const char* s) {
    int beg = 0;   // mark the begin of array
    int end = 0;

    KEY key;
    key.type = TEXT;
    while (true) {
      if (s[end] == '\0') {
        this->data.push_back(key);
        key.clear();
        /* if (mode != 0) { */
        /*   fprintf(stderr, "Malformated data %s!", s); */
        /*   return -1; */
        /* } */
        break;
      }

      switch (key.type) {
        case TEXT:
          if (s[end] == '$' && s[end-1] != '\\') { // change mode
            key.type = TEXT;
            this->data.push_back(key);
            key.clear();
            if (s[end + 1] == '[') {
              beg = end + 2;
              key.type = ARRAY;
            } else if (s[end + 1] == '(') {
              beg = end + 2;
              key.type = KEYWORD;
            } else {
              fprintf(stderr, "Wrong usage of $ in %s!\n", s);
              return -1;
            };
            end ++ ;
          } else {
            key.text.push_back(s[end]);
          };
          break;
        case KEYWORD:
          if (s[end] == ')' && s[end-1] != '\\'){  // change mode
            key.type = KEYWORD;
            key.keyword = std::string(s+beg, s+end);
            this->data.push_back(key);
            if (!isValidKeyword(key.keyword.c_str())) {
              fprintf(stderr, "Wrong keyword %s!\n", s);
            }
            key.clear();
            key.type = TEXT;
          } else {
            if ( !isalnum(s[end]) ) {
              fprintf(stderr, "Wrong keyword name: %c in %s", s[end], s);
              return -1;
            }
          };
          break;
        case ARRAY:
          if (s[end] == ']' && s[end-1] != '\\') { // change mode
            if (key.array.parse(s, beg, end) != 0) {
              fprintf(stderr, "Cannot parse this part: %s!\n", s+beg);
            }
            this->data.push_back(key);
            key.clear();
            key.type = TEXT;
          } else {
          }
          break;
      }
      end ++ ;
    };
    return 0;
  };
  void clear() {
    this->data.clear();
    this->dict.clear();
  }


  /**
   */
#if 0
public:
  StringTemplate(const char* t){
    this->temp = t;
  };
  StringTemplate(const std::string& t){
    this->temp = t;
  };
  /// totally reset StringTemplate
  void clear() {

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
#endif
};

int StringTemplate::Array::translate(std::string* str, const std::map<std::string, VALUE>& dict) {
  std::string& s = *str;
  std::map<std::string, VALUE>::const_iterator iter;

  unsigned int maxStringArraySize = 0;
  unsigned int idx = 0;
  bool reachValueEnd = true;
  do {
    reachValueEnd = true;
    if (idx)
      s += this->delim;
    for (unsigned int i = 0; i < data.size(); i++) {
      KEY& k = data[i];
      switch(k.type) {
        case TEXT:
          s += k.text;
          break;
        case KEYWORD:
          iter = dict.find(k.keyword);
          if (iter == dict.end()) {
            fprintf(stderr, "Data translation error: not found key %s!\n", k.keyword.c_str());
          } else {
            const VALUE& v = iter->second;
            if (v.type == STRING) {
              s += v.string;
            } else {
              if (idx < v.array.size() - 1) {
                reachValueEnd = false;
              } else {
                reachValueEnd = true;
              }
              if (idx == 0) { // make sure all string array have the same size
                if (maxStringArraySize == 0) {
                  if (v.array.size() > 0) {
                    maxStringArraySize = v.array.size();
                  }
                } else {
                  if (maxStringArraySize != v.array.size()) {
                    fprintf(stderr, "Unbalanced vector size. Stopped when tranlating %s!\n", k.keyword.c_str());
                    return -1;
                  };
                }
              }
              if (v.array.size() > 0) {
                s += v.array[idx];
              }
            };
          };
          break;
        case ARRAY:
          int ret = k.array.translate(str, this->dict);
          if (ret) {
            return ret;
          };
          break;
      };
    }
    idx ++;
  } while (! reachValueEnd);
  return 0;
};
int StringTemplate::Array::parse(const char* s, int param_beg, int param_end){
  int real_beg = param_beg;
  int real_end = param_end;
  for (int i = param_end; i >= param_beg; i--){
    if (s[i] == ' ') {
      this->delim = std::string(s + i + 1, s + param_end);
      real_end = i;
      break;
    }
  };
  int beg = real_beg;
  int end = real_beg;
  KEY* pKey = new KEY;
  KEY& key = *pKey;
  key.type = TEXT;
  while (true) {
    if (end == real_end) {
      this->data.push_back(key);
      key.clear();
      break;
    }

    switch (key.type) {
      case TEXT:
        if (s[end] == '$' && s[end-1] != '\\') { // change mode
          key.type = TEXT;
          this->data.push_back(key);
          key.clear();
          if (s[end + 1] == '[') {
            beg = end + 2;
            key.type = ARRAY;
          } else if (s[end + 1] == '(') {
            beg = end + 2;
            key.type = KEYWORD;
          } else {
            fprintf(stderr, "Wrong usage of $ in %s!\n", s);
            return -1;
          };
          end ++ ;
        } else {
          key.text.push_back(s[end]);
        };
        break;
      case KEYWORD:
        if (s[end] == ')' && s[end-1] != '\\'){  // change mode
          key.type = KEYWORD;
          key.keyword = std::string(s+beg, s+end);
          this->data.push_back(key);
          key.clear();
          key.type = TEXT;
          if (!isValidKeyword(key.keyword.c_str())) {
            fprintf(stderr, "Wrong keyword %s!\n", s);
          }
        } else {
          if ( !isalnum(s[end]) ) {
            fprintf(stderr, "Wrong keyword name: %c in %s", s[end], s);
            return -1;
          }
        };
        break;
      case ARRAY:
        if (s[end] == ']' && s[end-1] != '\\') { // change mode
          this->data.push_back(key);
          key.clear();
          key.type = TEXT;
          if (key.array.parse(s, beg, end) != 0) {
            fprintf(stderr, "Cannot parse this part: %s!\n", s+beg);
          }
        } else {
        }
        break;
    }
    end ++ ;
  };
  return 0;
};


#endif /* _STRINGTEMPLATE_H_ */
