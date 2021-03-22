



#ifndef __KSAMPLE_H_INCLUDE__
#define __KSAMPLE_H_INCLUDE__

#include <Containers/OutputStreamContainer.h>

using KSamplesBase = OutputStreamContainer;


class KSamples: public KSamplesBase 
{
public:
	KSamples() = default;
	~KSamples() {};
};

#endif /* ifndef __KSAMPLE_H_INCLUDE__ */




