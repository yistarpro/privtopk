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

    int32_t bootlv = 4; // 4 on default, 0 for no bootstrapping

    size_t numthread = (argc > 1) ? (size_t) std::stoul(argv[1]): 0;
    size_t initctnum = (argc > 2) ? (size_t) std::stoi(argv[2]) : 0;
    size_t iteration = (argc > 3) ? (size_t) std::stoi(argv[3]) : 8;



    if(numthread != 0){
        if(initctnum != 0){
            cout << "init number of ciphertexts: " << initctnum << ", number of threads: " << numthread <<  ", iteration: " << iteration << endl;
        	TopkTest_MT_opt(49, 1 << 17, 128, 2, bootlv, iteration, numthread, 64, initctnum, false, 2);
        }else{
            cout << "Operating experiments over all number of init ciphertexts, with numthread: " << numthread <<  ", iteration: " << iteration << endl;
            for(usint j=0;j<3;j++){ 
        	    TopkTest_MT_opt(49, 1 << 17, 128, 2, bootlv, iteration, numthread, 64, 8 << j, false, 2);
            }
        }
    }else{

        if(initctnum != 0){
            cout << "Operating experiments over all numthread, with init number of ciphertexts: " << initctnum <<  ", iteration: " << iteration << endl;
            for(usint i=0;i<8;i++){ 
        	    // TopkTest_MT_opt(49, 1 << 17, 128, 2, bootlv, iteration, 128 >> i, 64, initctnum, false, 2);
            }
        }else{
            cout << "Operating experiments over all numthread and init number of ciphertexts, iteration: " << iteration << endl;
            for(usint i=0;i<8;i++){ //varying the number of simultaneous operation
                TopkTest_MT_opt(49, 1 << 17, 128, 2, bootlv, iteration, 128 >> i, 64, 8, false, 2);
            }

            for(usint i=0;i<8;i++){ //varying the number of simultaneous operation
                TopkTest_MT_opt(49, 1 << 17, 128, 2, bootlv, iteration, 128 >> i, 64, 16, false, 2);
            }

            for(usint i=0;i<8;i++){ //varying the number of simultaneous operation
                TopkTest_MT_opt(49, 1 << 17, 128, 2, bootlv, iteration, 128 >> i, 64, 32, false, 2);
            }

            for(usint i=0;i<8;i++){ //varying the number of simultaneous operation
                TopkTest_MT_opt(49, 1 << 17, 128, 2, bootlv, iteration, 128 >> i, 64, 64, false, 4);
            }
        }

    }

 


    return 0;

}   