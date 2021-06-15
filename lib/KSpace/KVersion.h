//
// Created by shaune on 2021/3/23.
//

#ifndef KSPACESOLVER_KVERSION_H
#define KSPACESOLVER_KVERSION_H

#include <string>


inline std::string getVersion()
{
    const std::string version = "Version";
    const std::string majorNo = "2";
    const std::string minorNo = "20";

    return version + '-' + majorNo + '.' + minorNo;
}

#endif //KSPACESOLVER_KVERSION_H
