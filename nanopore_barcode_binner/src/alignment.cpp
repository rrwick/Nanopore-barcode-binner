
#include "alignment.h"

#include <iostream>

ScoredAlignment::ScoredAlignment(Align<Dna5String, ArrayGaps> & alignment,
                                 int readLength, int barcodeLength, int score):
    m_readLength(readLength), m_barcodeLength(barcodeLength),
    m_readStartPos(-1), m_barcodeStartPos(-1), m_rawScore(score)
{
    // Extract the alignment sequences into C++ strings for constant time random access.
    std::ostringstream stream1;
    stream1 << row(alignment, 0);
    std::string readAlignment =  stream1.str();
    std::ostringstream stream2;
    stream2 << row(alignment, 1);
    std::string barcodeAlignment =  stream2.str();

    int alignmentLength = std::max(readAlignment.size(), barcodeAlignment.size());
    if (alignmentLength == 0)
        return;

    int alignmentStartPos = -1;
    int alignmentEndPos = -1;

    // We consider the alignment to have started when we've encountered a base in both
    // sequences (though not necessarily at the same time).
    bool readStarted = false;
    bool barcodeStarted = false;
    for (int i = 0; i < alignmentLength; ++i) {
        if (readAlignment[i] != '-')
            readStarted = true;
        if (barcodeAlignment[i] != '-')
            barcodeStarted = true;
        if (readStarted && barcodeStarted) {
            alignmentStartPos = i;
            break;
        }
    }

    // We use the same logic to see when the alignment has ended.
    bool readEnded = false;
    bool barcodeEnded = false;
    for (int i = alignmentLength - 1; i >= 0; --i) {
        if (readAlignment[i] != '-')
            readEnded = true;
        if (barcodeAlignment[i] != '-')
            barcodeEnded = true;
        if (readEnded && barcodeEnded) {
            alignmentEndPos = i;
            break;
        }
    }

    if (alignmentStartPos == -1 || alignmentEndPos == -1)
        return;

    int readBases = 0, barcodeBases = 0;
    for (int i = 0; i < alignmentLength; ++i) {
        char base1 = readAlignment[i];
        char base2 = barcodeAlignment[i];

        if (i == alignmentStartPos) {
            m_readStartPos = readBases;
            m_barcodeStartPos = barcodeBases;
        }
        if (i == alignmentEndPos) {
            m_readEndPos = readBases;
            m_barcodeEndPos = barcodeBases;
        }

        if (base1 != '-')
            ++readBases;
        if (base2 != '-')
            ++barcodeBases;
    }
}

std::string ScoredAlignment::getString() {
    return std::to_string(m_readStartPos) + "," + 
           std::to_string(m_readEndPos) + "," + 
           std::to_string(m_barcodeStartPos) + "," + 
           std::to_string(m_barcodeEndPos) + "," + 
           std::to_string(m_rawScore);
}
