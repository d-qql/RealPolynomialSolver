//
// Created by d-qql on 24.10.2020.
//

#ifndef SPARCE_MATRICES_POLYNOMIAL_H
#define SPARCE_MATRICES_POLYNOMIAL_H
#include <map>
#include <cmath>
#include <vector>
#include "Tabs.h"
template<typename T>
struct Trip{
    std::pair<T, T> Segment;
    int count;
};
using namespace std;
template<typename T>
class Polynomial {
public:
    map<int, T> poly;
    Polynomial() = default;
    explicit Polynomial(map<int, T> input){
        this->poly = input;
    }
        //static vector<double> Solution(map<int, double> poly);
};
template<typename T>
pair<T, T> Lagrange(const Polynomial<T>& polynomial){
    T c = polynomial.poly.rbegin()->second;
    for(auto i : polynomial.poly){
        polynomial.poly[i.first]/=c;
        // cout<<i.second<<"x^"<<i.first<<" ";
    }
    int k=polynomial.poly.rbegin()->first;
    for(auto i = polynomial.poly.rbegin(); i != polynomial.poly.rend(); i++){
        if(i->second<0){
            k = polynomial.poly.rbegin()->first-i->first;
            break;
        }
    }
    T A = 0;
    for(auto i : polynomial.poly){
        if(i.second<0 && i.second<A) A = i.second;
    }
    //cout<<endl<<k<<" "<<A<<endl;
    double max = 1 + pow(-A, 1./k);

    for(auto i : polynomial.poly){
        if(i.first % 2 == 1){
            polynomial.poly[i.first] *= -1;
        }
    }
    for(auto i : polynomial.poly){
        //cout<<i.second<<"x^"<<i.first<<" ";
    }
    c = polynomial.poly.rbegin()->second;
    //cout<<c<<endl;
    for(auto i : polynomial.poly){
        polynomial.poly[i.first]/=c;
    }
    for(auto i : polynomial.poly){
        //cout<<i.second<<"x^"<<i.first<<" ";
    }
    k=polynomial.poly.rbegin()->first;
    for(auto i = polynomial.poly.rbegin(); i != polynomial.poly.rend(); i++){
        if(i->second<0){
            k = polynomial.poly.rbegin()->first-i->first;
            break;
        }
    }
    A = 0;
    for(auto i : polynomial.poly){
        if(i.second < 0 && i.second < A) A = i.second;
    }
    //cout<<endl<<k<<" "<<A;
    T min = -1 - pow(-A, T(1)/k);
    //cout<<endl<<pow(-A, 1./k)<<endl;
    return {min, max};
}

template<typename T>
int baseNumber(const vector<Polynomial<T>>& ShturmChain, pair<T, T> Segment){
    int changeA = 0, changeB = 0;
    T prevA = 0, prevB = 0, currentA = 0, currentB = 0;
    for(auto i : (ShturmChain.begin())->poly){
        prevA += pow(Segment.first, i.first)*i.second;
        prevB += pow(Segment.second, i.first)*i.second;
    }
    for(auto &k : ShturmChain){
        for(auto &i : k.poly){
            currentA += pow(Segment.first, i.first)*i.second;
            currentB += pow(Segment.second, i.first)*i.second;
        }
        if(currentA*prevA<T(0)) changeA++;
        if(currentB*prevB<T(0)) changeB++;
        prevA = currentA;
        prevB = currentB;
        currentA = T(0);
        currentB = T(0);
    }
    cout<<Tabs(changeA-changeB)<<endl;
    return Tabs(changeA-changeB);
}

template<typename T>
vector<Polynomial<T>> ShturmChain(const Polynomial<T>& polynomial, const Polynomial<T>& derivative){
    vector<Polynomial<T>> result;
    result.push_back(polynomial);
    result.push_back(derivative);
    for(int i = 2; i < polynomial.poly.rbegin()->first; i++){
        if(result[i-1].poly.empty()) break;
        Polynomial<T> temp = Division(result[i-2], result[i-1]);
        for(auto j : (temp.poly)){
            temp.poly[j.first]*=-1;
        }
        result.push_back(temp);
    }
    return result;
}

template<typename T>
int Budan_Fourier(const vector<Polynomial<T>>& derivatives, pair<T, T> Segment){
    int chA = 0, chB = 0;
    T poliA_prev = 0, poliB_prev = 0;
    for(auto j : derivatives.begin().poly) {
        poliA_prev += j.second * pow(Segment.first, j.first);
        poliB_prev += j.second * pow(Segment.second, j.first);
    }
    for(auto i = next(derivatives.begin()); i != derivatives.end(); i++){
        T poliA = 0, poliB = 0;
        for(auto &j : *i.poly){
            poliA += j.second*pow(Segment.first, j.first);
            poliB += j.second*pow(Segment.second, j.first);
        }
        if(poliA*poliA_prev<T(0)) chA++;
        if(poliB*poliB_prev<T(0)) chB++;
        poliA_prev = poliA;
        poliB_prev = poliB;
    }
    return chA-chB;
}

template<typename T>
Polynomial<T> Division(const Polynomial<T>& up, const Polynomial<T>& down){
    T c;
    while(up.poly.rbegin()!=up.poly.rend() && up.poly.rbegin()->first >= down.poly.rbegin()->first){
        c = T(up.poly.rbegin()->second)/T(down.poly.rbegin()->second);
        int deltaPOW = up.poly.rbegin()->first-down.poly.rbegin()->first;
        for(auto i:down.poly){
            if(up.poly.find(i.first+deltaPOW) != up.poly.end()){
                up.poly.find(i.first+deltaPOW)->second-=c*i.second;
                if(Tabs(up.poly.find(i.first+deltaPOW)->second) < T(1e-10)){
                    up.poly.erase(up.poly.find(i.first+deltaPOW));
                }
            }else{
                up.poly.insert({i.first+deltaPOW, -c*i.second});
            }
        }
    }
    return up;
}
template<typename T>
vector<Polynomial<T>> Derivatives(const Polynomial<T>& polynomial){
    vector<Polynomial<T>> result;
    result.push_back(polynomial);
    for(int i = 1; i <= polynomial.poly.rbegin()->first; i++){
        Polynomial<T> temp;
        for(auto & j : result[i-1]){
            if(j.first >= 1) {
                temp.poly.insert({j.first - 1, j.second * j.first});
            }
        }
        result.push_back(temp);
    }
    return result;
}

template<typename T>
Polynomial<T> Derivative(const Polynomial<T>& polynomial){
    Polynomial<T> result = Polynomial<T>();
    for(auto & j : polynomial.poly){
        if(j.first >= 1) {
            result.poly.insert({j.first - 1, j.second * j.first});
        }
    }
    return result;
}
template<typename T>
vector<T> oneRadicalSegments(const vector<Polynomial<T>>& ShturmChain, pair<T, T> Segment){
    vector<Trip<T>> result;
    vector<T> R;
    result.push_back({Segment, baseNumber(ShturmChain, Segment)});
    bool running = true;
    while (running) {
        running = false;
        for (auto &i : result) {
            if (i.count > 1) {
                running = true;
            }
        }
        if(!running) break;
        for (int i = 0; i < result.size(); i++) {
            if (result[i].count > 1) {
                int leftCount = baseNumber(ShturmChain, {result[i].Segment.first, (result[i].Segment.first + result[i].Segment.second) / 2});
                int rightCount = baseNumber(ShturmChain, {(result[i].Segment.first + result[i].Segment.second) / 2, result[i].Segment.second});
                if(leftCount > 0 || rightCount > 0) {
                    if (leftCount > 0)
                        result.push_back({{result[i].Segment.first, (result[i].Segment.first + result[i].Segment.second) / 2}, leftCount}); //Как же он сука был прав(пуш бэк изменяет итератор)
                    if (rightCount > 0)
                        result.push_back({{(result[i].Segment.first + result[i].Segment.second) / 2, result[i].Segment.second}, rightCount});
                    result.erase(result.begin()+i);
                }
            }
        }
    }
    running = true;
    while(running){
        running = false;
        for(auto &i: result){
            if (i.Segment.second-i.Segment.first>1e-6){
                running = true;
            }
        }
        if(!running) break;
        for(int i = 0; i < result.size(); i++){
            if(result[i].Segment.second-result[i].Segment.first>T(1e-6)){
                //cout<<result[i].Segment.second-result[i].Segment.first<<endl;
                int leftCount = baseNumber(ShturmChain, {result[i].Segment.first, (result[i].Segment.first + result[i].Segment.second) / 2});
                int rightCount = baseNumber(ShturmChain, {(result[i].Segment.first + result[i].Segment.second) / 2, result[i].Segment.second});
                if(leftCount == 1){
                    result.push_back({{result[i].Segment.first, (result[i].Segment.first + result[i].Segment.second) / 2}, leftCount});
                    result.erase(result.begin()+i);
                    continue;
                }
                if(rightCount == 1){
                    result.push_back({{(result[i].Segment.first + result[i].Segment.second) / 2, result[i].Segment.second}, rightCount});
                    result.erase(result.begin()+i);
                    continue;
                }
            }
        }
    }
    for(auto &i : result){
        R.push_back(i.Segment.first);
    }
    return R;
}
#endif //SPARCE_MATRICES_POLYNOMIAL_H
