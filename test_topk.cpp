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

    size_t size = (argc > 1) ? (size_t) std::stoul(argv[1]) : 0;
    size_t numk = (argc > 2) ? (size_t) std::stoi(argv[2]) : 0;
    size_t iteration = (argc > 3) ? (size_t) std::stoi(argv[3]) : 8;

    bool assertion = assert_pow2(size, 512, (1<<17));
    if(assertion == false)return 0;
    assertion = assert_pow2(numk, 1, 16);
    if(assertion == false)return 0;


    if(size != 0){
        if(numk != 0){
            cout << "numk: " << numk << "size: " << size <<  "iteration: " << iteration << endl;
            TopkTest_MT_opt(49, size, 128, numk, 4, iteration, 128, 64, 128, true);
        }else{
            cout << "Operating experiments over all k, with size: " << size <<  "iteration: " << iteration << endl;
            for(usint j=0;j<5;j++){ 
                TopkTest_MT_opt(49, size, 128, 1 << j, bootlv, iteration, 128, 64, 128, true);
            }
        }
    }else{

        if(numk != 0){
            cout << "Operating experiments over all size, with k: " << numk <<  "iteration: " << iteration << endl;
            for(usint i=11;i<18;i++){ 
                TopkTest_MT_opt(49, 1 << i, 128, numk, bootlv, iteration, 128, 64, 128, true);
            }
        }else{
            cout << "Operating experiments over all size and k, iteration: " << iteration << endl;
            for(usint j=0;j<5;j++){ 
                for(usint i=11;i<18;i++){ 
                    TopkTest_MT_opt(49, 1 << i, 128, 1 << j, bootlv, iteration, 128, 64, 128, true);
                }
            }
            cout << "Additional experiments for optimized param, iteration: " << iteration << endl;
            // Exceptions:
            TopkTest_MT_opt(49, 1 << 17, 128, 1, bootlv, iteration, 128, 64, 16, false, 4);
            TopkTest_MT_opt(49, 1 << 16, 128, 2, bootlv, iteration, 128, 64, 16, false, 4);
        }

    }

 

    return 0;

}   