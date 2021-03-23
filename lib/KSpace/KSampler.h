//
// Created by shaune on 2021/3/22.
//

#ifndef KSPACESOLVER_KSAMPLER_H
#define KSPACESOLVER_KSAMPLER_H

#include <vector>

class KSampler {
public:
    virtual ~KSampler();

    KSampler();

private:
    using KSamplerContainerType = std::vector<std::vector<float>>;
};

#endif //KSPACESOLVER_KSAMPLER_H
