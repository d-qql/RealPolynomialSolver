//
// Created by d-qql on 24.10.2020.
//

#ifndef SPARCE_MATRICES_TABS_H
#define SPARCE_MATRICES_TABS_H
template<typename T>
T Tabs(T num){
    if(num<T(0)) return -num;
    else return num;
}
#endif //SPARCE_MATRICES_TABS_H
