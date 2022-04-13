/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   StringOps.cc
 * Author: bschroed
 * 
 * Created on December 18, 2017, 4:04 PM
 */

#include "StringOps.h"

namespace utils{
    string StringOps::replaceall(string input, string replacepattern, string insertpattern){

        int replacepattternlength= replacepattern.length();
        int startpos = 0;
        int found_pos = input.find( replacepattern, startpos); //initial search
        ostringstream finalstr;

        //get all start positions
        while(found_pos != string::npos){
            finalstr << input.substr(startpos, found_pos-startpos) << insertpattern; //found a position, build up new string in stream final
            startpos = found_pos + replacepattternlength; // set new star position for search
            found_pos = input.find( replacepattern, startpos); // get new position
        }

        finalstr << input.substr(startpos, input.length()-startpos);    //append suffix
        return finalstr.str();
    }
}
