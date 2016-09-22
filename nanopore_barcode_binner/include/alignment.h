
#ifndef ALIGNMENT_H
#define ALIGNMENT_H


#include <string>
#include <seqan/basic.h>
#include <seqan/align.h>

using namespace seqan;


class ScoredAlignment {
public:
    ScoredAlignment(Align<Dna5String, ArrayGaps> & alignment,
                    int readLength, int barcodeLength, int score);
    std::string getString();

    int m_readLength;
    int m_barcodeLength;
    int m_readStartPos;
    int m_readEndPos;
    int m_barcodeStartPos;
    int m_barcodeEndPos;
    int m_rawScore;
};

#endif // ALIGNMENT_H
