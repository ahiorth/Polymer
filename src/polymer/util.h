#ifndef UTIL_H
#define UTIL_H

#include <vector>
template < typename T>
std::pair<bool, int > findInVector(const std::vector<T>  & vecOfElements, const T  & element)
{
    std::pair<bool, int > result;
    // Find given element in vector
    auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);
    if (it != vecOfElements.end())
    {
        result.second = distance(vecOfElements.begin(), it);
        result.first = true;
    }
    else
    {
        result.first = false;
        result.second = -1;
    }
    return result;
}

// Example of Use
//if (result.first)
//std::cout << "Element Found at index : " << result.second << std::endl;
//else
//std::cout << "Element Not Found" << std::endl;


template < typename T>
std::pair<int, T > findInRange(const std::vector<T>& vecOfElements, const T& element)
{
    std::pair<int, T > result;
    std::vector<T>::const_iterator it;
    T prev_element = 0;
    // Find where the element is located in a list of increasing elements
    // e.g. if vecOfElements = [0,1,2,3,6,10]
    // findInRange(vecOfElements, 5) will return true, 4
    for(it=vecOfElements.begin();it!=vecOfElements.end();++it)
    {
        if (*it >= element)
        {
            result.first = std::distance(vecOfElements.begin(), it);
            result.second = *it-prev_element;
            return result;
        }
        prev_element = *it;
    }
    result.second = -1;
    result.first  = -1;
    return result;
}




#endif
