// =====================================================================================
//
//       Filename:  myutil.h
//
//    Description:
//
//        Version:  1.0
//        Created:  02/27/2018 04:53:23 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
//        Company:  unemployed
//
// =====================================================================================
/*************************************************************************
 *
 *
 * several frequently used functions
 *
 *
 *
************************************************************************/

#ifndef _MYUTIL_H
#define _MYUTIL_H
#include <string>
#include <cassert>
#include <vector>
#include "../model/My_Vector.h"
#include <fstream>
#include <regex>
#include <sstream>
#include <iostream>
#include <sys/time.h>



void split(const std::string &s, std::vector<std::string> &elements);
//std::vector<std::string> split(const std::string &s, char delim);
void clear_ram(double **a, const int &m);
void clear_ram(int **a, const int &m);
const std::string currentDateTime();
struct greater {
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};

template <typename T>
const int argmax(const My_Vector<T> & a){
    int max_i = 0;
    double max = a.get_array()[max_i];
    for ( int i=1; i<a.get_length(); ++i ){
        if ( a.get_array()[i] > max ){
            max_i=i;
            max = a.get_array()[max_i];
        }
    }
    return max_i;
}

size_t getFileLineNumber ( std::string & filePath );
#endif

