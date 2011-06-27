#ifndef _FREQTABLE_H_
#define _FREQTABLE_H_


template <class T>
class FreqTable{
public:
    void add(const T& t) {
        if (this->data.find(t) == this->data.end()) {
            this->data[t] = 1;
        } else {
            this->data[t] ++; 
        }
        this->isSorted = false;
    };
    void remove(const T& t) {
        if (this->data.find(t) == this->data.end()) {
            return false;
        }
        this->data[t] -- ;
        this->isSorted = false;
    };
    unsigned int size() const{ return this->data.size();}; 
    // return the frequency in ascending order
    void at(const unsigned int idx, T* t, int* v) {
        if (!this->isSorted) 
            this->sort();
        *v = this->orderedData[idx].first;
        *t = *(this->orderedData[idx].second);
    };
    void clear() {
        this->data.clear();
        this->orderedData.clear();
        this->isSorted = false;
    };
private:
    void sort() {
        this->orderedData.clear();
        typename std::map<T, int >::iterator it;
        for (it = this->data.begin(); 
             it != this->data.end() ; it++) {
            this->orderedData.push_back(std::make_pair( (*it).second, &((*it).first)) );
        }
        std::sort(this->orderedData.begin(), this->orderedData.end());
        this->isSorted = true;
    };
    std::map< T, int > data;
    std::vector< std::pair<int, const T*> > orderedData;
    bool isSorted;
 };

//////////////////////////////////////////////////////////////////////
// Test case 1 //
/*
    FreqTable<std::string> codonFreq;
    codonFreq.add("A");
    codonFreq.add("B");
    codonFreq.add("C");
    codonFreq.add("D");
    codonFreq.add("A");
    codonFreq.add("B");
    codonFreq.add("C");
    codonFreq.add("A");
    codonFreq.add("B");
    codonFreq.add("A");
    
    std::string s;
    int f;
    for (unsigned int i = 0 ; i < codonFreq.size(); i++) {
        codonFreq.at(i, &s, &f);
        printf("freq of %s is %d\n", s.c_str(), f);
    };
    return 0;
*/

// Test case 2 //
/*
    FreqTable<char> codonFreq;
    codonFreq.add('A');
    codonFreq.add('B');
    codonFreq.add('C');
    codonFreq.add('D');
    codonFreq.add('A');
    codonFreq.add('B');
    codonFreq.add('C');
    codonFreq.add('A');
    codonFreq.add('B');
    codonFreq.add('A');
    
    char s;
    int f;
    for (unsigned int i = 0 ; i < codonFreq.size(); i++) {
        codonFreq.at(i, &s, &f);
        printf("freq of %c is %d\n", s, f);
    };
    return 0;

*/

#endif /* _FREQTABLE_H_ */
