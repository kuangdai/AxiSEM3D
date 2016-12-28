// NullSource.h
// created by Kuangdai on 1-Nov-2016 
// null source, for scaling test only

#pragma once
#include "Source.h"

class NullSource: public Source {
public:
    NullSource(): Source() {};
    std::string verbose() const;
    
protected:    
    void computeSourceFourier(const Quad &myQuad, const RDColP &interpFactZ,
        arPP_CMatX3 &fouriers) const;
};

