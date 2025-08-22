#define PROFILE  // turns on the reporting of timing results

#include "openfhe.h"
#include "utils.h"
#include "testcode.h"
#include "algorithms.h"


#include "schemerns/rns-leveledshe.h"

#include <iostream>
#include <getopt.h>
// #include <ctime>


using namespace lbcrypto;
using namespace std;
using namespace ckkseif;

int main(int argc, char **argv) {

    size_t size = (argc > 1) ? (size_t) std::stoul(argv[1]) : 0;
    size_t iteration = (argc > 2) ? (size_t) std::stoi(argv[2]) : 8;

    bool assertion = assert_pow2(size, 8, 256);
    if(assertion == false)return 0;

    if(size != 0){
        cout << "size: " << size <<  ", iteration: " << iteration << endl;
        SortTest(45,size,128, 0, 16, iteration); // sf, size, domain, bootlevel, iter
    }else{

        cout << "Operating experiments over all sizes: iteration: " << iteration << endl;
        for(usint i=3;i<9;i++){ 
            SortTest(45,1 << i,128, 0, 16, iteration); // sf, size, domain, bootlevel, iter
        } 

    }

 

    return 0;

}         