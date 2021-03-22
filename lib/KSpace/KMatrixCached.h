



#ifndef __KCACHED_MATRIX_H_INCLUDE__
#define __KCACHED_MATRIX_H_INCLUDE__ 


#include <KSpace/KMatrix.h>
#include <unordered_map>

using KMatrixCached = std::unordered_map<std::string, 
	                                     std::shared_ptr<KBaseMatrix> >;

#endif /* ifndef __KCACHED_MATRIX_H_INCLUDE__ */




