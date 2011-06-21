#ifndef _ORDEREDMAP_H_
#define _ORDEREDMAP_H_

template <class KEY, class TYPE>
class OrderedMap{
  public:
    bool find(const KEY& key) {
        if (this->keyTypeMap.find(key) == this->keyTypeMap.end()){
            return false;
        }
        return true;
    }
    TYPE& operator[] (const KEY& key) {
        if (!this->find(key)){
            this->keyVec.push_back(key);
        }
        return this->keyTypeMap[key];
    }
    void front(KEY* k, TYPE* v) {
        *k = this->keyVec.front();
        *v = this->keyTypeMap[(*k)];
    }
    bool at(unsigned int idx, KEY* k, TYPE* v) {
        if (idx >= this->size()) false;
        *k = this->keyVec[idx];
        *v = this->keyTypeMap[(*k)];
        return true;
    }
    unsigned int size() const { return this->keyVec.size();} ;
  private:
    std::vector < KEY > keyVec;
    std::map < KEY, TYPE > keyTypeMap;
};

#endif /* _ORDEREDMAP_H_ */
