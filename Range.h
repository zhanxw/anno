#ifndef _RANGE_H_
#define _RANGE_H_

// all are 1-based index, inclusive on boundaries.
// NOTE: UCSC use 0-based start and 1-based end, this is convenient to calculate length
//       but for complex case, this causes confusion.
struct Range{
    int start;
    int end;
    Range(int s, int e): start(s), end(e) {};
    Range(): start(-1), end(-1) {};
    int length() const { 
        int l = end - start + 1;
        if (l < 0) {
            printf("getLength() < 0 for start(%d) and end(%d)\n", start, end);
        }
        return l;
    };
    bool inRange(int pos) {
        if (this->start <= pos && pos <= this->end) 
            return true;
        return false;
    };
};
#endif /* _RANGE_H_ */
