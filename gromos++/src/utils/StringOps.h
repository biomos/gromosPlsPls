/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   StringOps.h
 * Author: bschroed
 *
 * Created on December 18, 2017, 3:46 PM
 */



#ifndef STRINGOPS_H
#define STRINGOPS_H

#include <string>
#include <list>
#include <iostream>
#include <sstream>

using namespace std;

namespace utils {
/**
     * @class StrigOps
     * 
     * contains some internal string operations
     * @param
     * @return 
     */
    class StringOps {
    public:
        //StringOps();
        //StringOps(const StringOps& orig);
        //virtual ~StringOps();

        //function replaces all replacepattern occurences with insertpatterns
        static string replaceall(std::string input, std::string replacepattern, std::string insertpattern);

    private:

    };
}
#endif /* STRINGOPS_H */
