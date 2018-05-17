// =====================================================================================
//
//       Filename:  myutil.cpp
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


#include "myutil.h"

void split(const std::string &s, std::vector<std::string> &elements){
    std::istringstream iss(s);
    do {
        std::string subs;
        iss >> subs;
        if (subs.length() > 0) {
            elements.push_back(subs);
        }
    } while (iss);
}

void clear_ram(double **a, const int &m){
    for( int i=0; i<m; ++i ){
        delete [] a[i];
    }
    delete[] a;
}

const std::string currentDateTime() {
    char timeBuf  [256];
    struct timeval tv;
    struct timezone tz;
    struct tm *tm;
    gettimeofday(&tv, &tz);
    tm=localtime(&tv.tv_sec);
    sprintf (timeBuf, "%02d:%02d:%02d:%03d",tm->tm_hour, tm->tm_min, tm->tm_sec, (tv.tv_usec/1000) );
    return timeBuf;
}

size_t getFileLineNumber ( std::string & filePath ){
    // check file begin
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening file " << filePath << std::endl;
        exit (1);
    }
    // check file end
    infile.clear();
    // count variants and individuals begin
    size_t number_of_line = 0;
    std::string line;
    while (std::getline(infile, line)){
        if( line.length() > 0 ){
            ++number_of_line;
        }
    }
    infile.close();
    return number_of_line;
}
