/*
    This file is part of HMMLab.

    HMMLab is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HMMLab is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HMMLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <map>
#include <vector>
#include <algorithm>
using namespace std;

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

template <typename T>
class List: public vector<T>
{
public:
    typedef typename List<T>::iterator list_iterator;

    List(): vector<T>() {};

    ~List() {}

    List(int size, T value): vector<T>() {
        for(int i = 0; i < size; i++) {
            append(value);
        }
    };

    List<T>& operator+=(List<T>& list) {
        this->insert(this->end(), list.begin(), list.end());
        return *this;
    };

    List<T>& operator+(List<T>& list) {
        List<T>* result = new List<T>();
        *result = *this;
        return (*result += list);
    };

    void append(T elem) {
        this->push_back(elem);
    };

    void remove(int i) {
        this->erase(List<T>::begin() + i);
    };

    void remove_value(T value) {
        for(list_iterator it = List<T>::begin(); it < List<T>::end(); it++) {
            if((*it) == value) {
                this->erase(it);
            }
        }
    };

    int index(T value) {
        int i = 0;
        for(list_iterator it = List<T>::begin(); it < List<T>::end(); it++, i++) {
            if((*it) == value) {
                return i;
            }
        }
        return -1;
    };

    T* c_array() {
        T* result = new T[List<T>::size()];
        int i = 0;
        for(list_iterator it = List<T>::begin(); it < List<T>::end(); it++) {
            result[i] = *it;
            i++;
        }
        return result;
    };

    void listsort() {
        sort(this->begin(), this->end());
    };
};

template <typename T, typename U>
class Dict: public map<T, U>
{
public:
    typedef typename Dict<T, U>::iterator dict_iterator;

    Dict(): map<T, U>() {};

    ~Dict() {};

    void update(Dict<T, U>& dict) {
        for(dict_iterator it = dict.begin(); it != dict.end(); it++) {
            (*this)[it->first] = it->second;
        }
    };

    List<T> keys() {
        List<T> result;
        for(dict_iterator it = Dict<T, U>::begin(); it != Dict<T, U>::end(); it++) {
            result.append(it->first);
        }
        return result;
    };
};

template <typename T>
T sum(List<T> list)
{
    T res = 0;
    for(unsigned int i = 0; i < list.size(); i++) {
        res += list[i];
    }
    return res;
};

template <typename T>
List<T> set2List(set<T> s)
{
    List<T> result;
    copy(s.begin(), s.end(), back_inserter(result));
    return result;
};

#endif
