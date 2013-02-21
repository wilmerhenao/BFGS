#ifndef _CONTAINER_HPP_
#define _CONTAINER_HPP_

#include<iostream>
#include<string>
#include<map>

class StringComparerForMap{
public:
    bool operator()(const std::string x, const std::string y){
	if(x.compare(y) > 0)
	    return true;
	return false;
    }
}

// Declares a map with keys as strings that contains pointers to functions    
typedef std::map<std::string, void(*)(T*, T*, T*, int), StringComparerForMap> myMap;

#endif // _CONTAINER_HPP_
