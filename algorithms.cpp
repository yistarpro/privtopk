#include "openfhe.h"
#include "utils.h"
#include "algorithms.h"
#include <iostream>
#include <vector>
#include <cmath>
// #include "openfhecore.h"
// #include "math/chebyshev.h"
#include <map>
#include <omp.h>

using namespace lbcrypto;
using namespace std;

namespace ckkseif {

std::vector<int32_t> GenIdxForRotsum(const int32_t from, const int32_t to, const bool bidirectional){
    if(bidirectional){
        const int32_t round = log2(from/to);
        std::vector<int32_t> arr(2*round);

        int32_t rotnum=from;
        for(int32_t s=1 ; s < round+1 ; s++){
            arr[s-1] = rotnum >> s;
            arr[round+s-1] = - arr[s-1];
        }
        return arr;
    }else{
        const int32_t round = log2(from/to);
        std::vector<int32_t> arr(round);

        int32_t rotnum=from;
        for(int32_t s=1 ; s < round+1 ; s++){
            arr[s-1] = rotnum >> s;
        }
        return arr;
    }
    
}

std::vector<int32_t> GenIdxForMultiples(const int32_t base, const int32_t number, const bool bidirectional){
    if(bidirectional){
        std::vector<int32_t> arr(2*number);

        for(int32_t s=1 ; s < number+1 ; s++){
            arr[s-1] = base * s;
            arr[number+s-1] = - arr[s-1];
        }
        return arr;
    }else{
        std::vector<int32_t> arr(number);

        for(int32_t s=1 ; s < number+1 ; s++){
            arr[s-1] = base * s;
        }
        return arr;
    }
    
}



void AddRotKeyForPo2(const PrivateKey<DCRTPoly> privateKey, CryptoContext<DCRTPoly> cc, const int32_t batchSize){
    int32_t copy=batchSize;
    std::vector<int32_t> arr(2*log2(batchSize));
    for(long i = 0 ; i < log2(batchSize) ; i++){
        copy >>= 1;
        arr[i]=(copy);
        arr[log2(batchSize)+i]=-(copy);
    }
    cc->EvalRotateKeyGen(privateKey, arr);

}

void AddRotKeyForSort(const PrivateKey<DCRTPoly> privateKey, CryptoContext<DCRTPoly> cc, const int32_t arrsize){
    int32_t arrsizesquare=arrsize*arrsize;
    int32_t arrsizesquare2=arrsize*arrsize-arrsize;
    int32_t interval = 2*log2(arrsize);
    
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize();
    int32_t numct = arrsizesquare / batchSize;

    std::vector<int32_t> arr(interval*3+numct);
    for(int32_t i = 0 ; i < interval ; i++){
        arrsizesquare >>= 1;
        arr[i]=(arrsizesquare);
        arr[i+interval]= -(arrsizesquare);
    }
    for(usint i = 0 ; i < log2(arrsize) ; i++){
        arrsizesquare2 >>= 1;
        arr[i+2*interval]= (arrsizesquare2);
        arr[i+2*interval+log2(arrsize)]= -(arrsizesquare2);
    }

    //For sort_multi:
    if(numct!=0){
        int32_t slicelength = arrsize / numct;
        for(int32_t i=0;i<numct;i++){
            arr[3*interval+i] = (slicelength*i);
        }
    }

    cc->EvalRotateKeyGen(privateKey, arr);
}

void AddRotKeyForTopk_multithread_opt(const PrivateKey<DCRTPoly> privateKey, CryptoContext<DCRTPoly> cc, const int32_t size, const int32_t numk, const int32_t arrsize_init, const int32_t mergecrit){
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize();
    bool verbose = false;


    // int32_t arrsize = 2*numk;
    int32_t arrsize = arrsize_init;
    int32_t currentsize = size;

    // int32_t numct = 2*numk*size/ batchSize;
    int32_t numct = size * arrsize_init / batchSize;

    // if(batchSize * threads_actual > size * numk){ // threads_actual - numbered ciphertext input, if initial arrsize >  numk
    //     arrsize = batchSize * threads_actual / size;
    //     if(arrsize > (1 << 6)){ //arrsize is too large for taking topk_large as initial phase
    //         threads_actual = size / (batchSize >> 6);
    //         if(threads_actual ==0)threads_actual = 1; 
    //         arrsize = batchSize * threads_actual / size;
    //         cout << "threads altered: " << threads_actual << endl;
    //     }
    //     numct = threads_actual;
    // }

    if(verbose)cout << "Initial arrsize: " << arrsize << ", numct: " << numct  << endl;

    //rankselect shares po2 keys for various arrsize, and we generate it here at once.
    std::vector<int32_t> arrpo2 = GenIdxForRotsum(batchSize, 1, true);
    if(verbose)cout << arrpo2 << endl;
    cc->EvalRotateKeyGen(privateKey, arrpo2);


    while(numct > 1){

        if(verbose)cout << "numct: " << numct << ", arrsize: " << arrsize << ", num element: " << currentsize << endl;
        //other rotsums require po2 keys, and transpose keys only matters
        //transpose keys for rankselect
        std::vector<int32_t> arr = GenIdxForRotsum(arrsize*arrsize-arrsize, arrsize-1, true);
        cc->EvalRotateKeyGen(privateKey, arr);
        if(verbose)cout << arr << endl;

        int32_t next_size = currentsize * numk / arrsize;
        int32_t next_arrsize = 64;

        // //condition check
        // if(numk * next_size < batchSize){ //merge to one
        //     multiplier = numct;
        // }else{
        //     int32_t target_output_bound = min(256*next_size / batchSize, boot_threshold); // number of output: (condition 2, empty sub-matrix occurence) (condition 3, initial merge for threads_actual)
        //     if( target_output_bound < numct ){
        //         multiplier = numct / target_output_bound; 
        //         next_arrsize = target_output_bound * batchSize / next_size;
        //     }
        // }
        // if(binarymerge){
        //     multiplier = 2;
        //     next_arrsize = arrsize;
        // }

        
   

        // ================================New: rotkeys for merge and divide
        int32_t numpacking = max(numct / arrsize, 1); //number of packing ciphertexts for bootstrapping
        int32_t numpacked = numct/numpacking; //number of ciphertexts packed in a intermediate ciphertext

        //Rankvec
        std::vector<int32_t> arrm1  = GenIdxForMultiples(arrsize, numpacked, true);
        cc->EvalRotateKeyGen(privateKey, arrm1);
        if(verbose)cout << "merge and divide after comp" << arrm1 << endl;
        // ================================New: rotkeys for merge and divide

        
        if(numk * next_size * mergecrit < batchSize){
            next_arrsize = min(next_size, batchSize / next_size);
            int32_t arrsize_square = arrsize * arrsize ;
            int32_t next_arrsizesquare = next_arrsize * next_arrsize ;

            if(verbose)cout << "merge to one, arrsize: " << arrsize << ", next_arrsize: "<< next_arrsize <<  endl;;
            //Merge Keys                    


            //Variables for merge
            int32_t var0 = max(numk / max(next_arrsizesquare / arrsize, 1),1);
            int32_t var1 = next_arrsize * var0 / numk; 
            int32_t var2 = var0 * max(next_arrsizesquare / arrsize ,1) * arrsize;
            int32_t var3 = max(var0, next_arrsize / arrsize);

            std::vector<int32_t> arr2(numct-1);
            for(int32_t idx = 1 ; idx < numct ; idx++){
                int32_t idx0 = idx % var1;
                int32_t idx1 = idx / var1;
                arr2[idx-1]=  -(idx0 + var2 * idx1);
            }
            if(verbose)cout << "topk_large_merge Key 1 " << arr2 << endl;
            cc->EvalRotateKeyGen(privateKey, arr2);        

            int32_t rotvar0 = next_arrsize / (numk * numct);
            int32_t rotvar1 = numk / var3;
            if(rotvar0 > 1){
                std::vector<int32_t> arr3 = GenIdxForRotsum((arrsize_square - numct) * rotvar0, (arrsize_square - numct));
                if(verbose)cout << "topk_large_merge Key 2 " << arr3 << endl;
                cc->EvalRotateKeyGen(privateKey, arr3);   
            }

            if(rotvar1 > 1){
                std::vector<int32_t> arr4 = GenIdxForRotsum((max(arrsize, next_arrsize) - var1) * rotvar1, (max(arrsize, next_arrsize) - var1));
                if(verbose)cout << "topk_large_merge Key 3 " << arr4 << endl;
                cc->EvalRotateKeyGen(privateKey, arr4);   
            }


            //update
            numct = 1;
        }else{
            //gathering keys for blocktopks
            // std::vector<int32_t> arr2 = GenIdxForRotsum((arrsize*arrsize - numk ) * next_arrsize/numk, (arrsize*arrsize - numk ));
            // if(verbose)cout << "Without merge " << arr2 << endl;
            // cc->EvalRotateKeyGen(privateKey, arr2);

            // ================================New: rotkeys for merge and divide
            next_arrsize = min(arrsize * arrsize / numk, 256);

            numpacking = max(numct * numk / (arrsize * arrsize), 1); //number of packing ciphertexts for bootstrapping
            numpacked = numct / numpacking; //number of ciphertexts packed in an intermediate ciphertext
            int32_t var1 = arrsize /  min(arrsize, numct); //required number of blocks to fill a row
            int32_t var2 = arrsize /  (numk * max(numpacked / arrsize, 1)); //required number of blocks to fill columns

            //Gather
            std::vector<int32_t> arrm2(numpacked);
            for(int32_t idx=1;idx< numpacked;idx++){
                arrm2[idx-1] = - (idx % arrsize + numk * arrsize * (idx / arrsize));
            }
            cc->EvalRotateKeyGen(privateKey, arrm2);
            if(verbose)cout << "merge and divide blocktopkS _ 1" << arrm2 << endl;
        
            //Rearrange
            if(var1 > 1){
                std::vector<int32_t> arrm3 = GenIdxForRotsum(((next_arrsize*next_arrsize/var1 - numct)*var1), (next_arrsize*next_arrsize/var1 - numct), true); //gather to construct full rows
                cc->EvalRotateKeyGen(privateKey, arrm3);
                if(verbose)cout << arrm3 << endl;
            }
            if(var2 > 1){
                std::vector<int32_t> arrm3 = GenIdxForRotsum(((next_arrsize*next_arrsize/(var1*var2) - numk*arrsize* max(numpacked / arrsize, 1)))*var2, (next_arrsize*next_arrsize/(var1*var2) - numk*arrsize* max(numpacked / arrsize, 1)), true);  //gather to construct full columns
                cc->EvalRotateKeyGen(privateKey, arrm3);
                if(verbose)cout << arrm3 << endl;
            }
            
            //Extraction
            int32_t resultnumct = numct * numk * next_arrsize / (arrsize * arrsize);
            int32_t resultnumpacked = resultnumct / numpacking; //number of result ciphertexts to be divdided per am intermediate ciphertext

            std::vector<int32_t> arrm4 = GenIdxForMultiples(next_arrsize, resultnumpacked, false);
            cc->EvalRotateKeyGen(privateKey, arrm4);
            if(verbose)cout << "merge and divide blocktopkS _ 4" << arrm4 << endl;

            numct = resultnumct;
            //New: rotkeys for merge and divide

        }
        
        currentsize = next_size;
        arrsize = next_arrsize;

    }


    // Single Ciphertext Case
    ///////////*********************///////////
    if(currentsize > (1 << 12) && numk==1){
        //Gen rotkey for gathering after iterative Max, for size = 2^13 ~ 2^15.
        std::vector<int32_t> arr = GenIdxForRotsum(240, 30);// v2 * (v2 * v1 * v1 - v1);, for v1 = batchSize / size = initial arrsize, and v2 = size / 2^3
        if(verbose)cout << "Iter Max for: " << currentsize << " : "<< arr << endl;
        cc->EvalRotateKeyGen(privateKey, arr);

        arrsize = 16;
        currentsize = 1 << 12;
    }
    ///////////*********************///////////
    
    //Keys for blocktopks
    while(currentsize > 256){
        if(verbose)cout << "arrsize: " << arrsize << ", num element: " << currentsize << endl;
        //transpose keys for rankselect
        std::vector<int32_t> arr = GenIdxForRotsum(arrsize*arrsize-arrsize, arrsize-1, true);
        cc->EvalRotateKeyGen(privateKey, arr);
        if(verbose)cout << arr << endl;

        int32_t next_arrsize = arrsize*arrsize/numk;
        if(next_arrsize * next_arrsize > batchSize) next_arrsize = batchSize / next_arrsize;
        currentsize = currentsize * numk / arrsize;

        //gathering keys
        std::vector<int32_t> arr2 = GenIdxForRotsum((arrsize*arrsize - numk ) * next_arrsize/numk, (arrsize*arrsize - numk ));
        if(verbose)cout << "true " << arr2 << endl;
        cc->EvalRotateKeyGen(privateKey, arr2);

        //update
        arrsize = next_arrsize;
    }

    if(verbose)cout << "Final arrsize: " << arrsize << endl;

    //transpose keys for rankselect
    std::vector<int32_t> arr2 = GenIdxForRotsum(arrsize*arrsize-arrsize, arrsize-1, true);
    cc->EvalRotateKeyGen(privateKey, arr2);
    if(verbose)cout << "Final key" << arr2 << endl;

}



void AddRotKeyForBoot(const PrivateKey<DCRTPoly> privateKey, CryptoContext<DCRTPoly> cc, const int32_t size){
    int32_t sizesquare=size*size;
    int32_t sizesquare2=size*size-size;
    int32_t interval = 2*log2(size);
    
    std::vector<int32_t> arr(interval*3);
    for(int32_t i = 0 ; i < interval ; i++){
        sizesquare >>= 1;
        arr[i]=(sizesquare);
        arr[i+interval]= -(sizesquare);
    }
    for(usint i = 0 ; i < log2(size) ; i++){
        sizesquare2 >>= 1;
        arr[i+2*interval]= (sizesquare2);
        arr[i+2*interval+log2(size)]= -(sizesquare2);
    }
    cc->EvalRotateKeyGen(privateKey, arr);
}

Ciphertext<DCRTPoly> BootAuto(Ciphertext<DCRTPoly> ciphertext){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    const auto cryptoParamsCKKS =
    std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(
          cc->GetCryptoParameters());

	usint RingDim = log2(cc->GetRingDimension());
    usint levelBudgetElmt= (RingDim >15 ) ? 1 << (RingDim-14) : 2 ;
    std::vector<uint32_t> levelBudget = {levelBudgetElmt, levelBudgetElmt};

    Ciphertext<DCRTPoly> result;
    usint depth = FHECKKSRNS::GetBootstrapDepth(levelBudget, cryptoParamsCKKS->GetSecretKeyDist());
    usint current = ciphertext->GetLevel();
    if(depth >= current)
        result = ciphertext->Clone();
    else{
        result = cc->EvalBootstrap(ciphertext);
    }
    return result;
}

Ciphertext<DCRTPoly> BootWithPrec(Ciphertext<DCRTPoly> ciphertext, const usint verbose, const KeyPair<DCRTPoly> keys){
    Ciphertext<DCRTPoly> result;
    Plaintext pt;
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    const int32_t batchSize = cc->GetEncodingParams()->GetBatchSize();
    vector<double> vals1(batchSize);

    // cout << "Boot: " <<  ciphertext->GetLevel() << endl;
    if(verbose == 1){
        cc->Decrypt(keys.secretKey, ciphertext, &pt);
        vals1 = pt->GetRealPackedValue();
    }
    if(verbose == 2){
        std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
        pt->SetLength(256);
        cout << pt << endl;
    }
    result = cc->EvalBootstrap(ciphertext,2,10);
    if(verbose == 1){
        cc->Decrypt(keys.secretKey, result, &pt);
    }
    if(verbose == 2){
        cc->Decrypt(keys.secretKey, result, &pt);
        std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
        pt->SetLength(256);
        cout << pt << endl;

    }
    if(verbose == 1){
        precision(pt, vals1, batchSize);
    }

    return result;
}


Ciphertext<DCRTPoly> Product(const vector<Ciphertext<DCRTPoly>> ciphertext){
    usint num = ciphertext.size();
    usint phase = ceil(log2(num));
    usint res = 0;
    const CryptoContext<DCRTPoly> cc = ciphertext[0]->GetCryptoContext();

    vector<Ciphertext<DCRTPoly>> result(num);
    for(usint i=0; i<num ; i++){
        result[i] = ciphertext[i]->Clone();
    }
    for(usint i=0; i< phase; i++){
        res = num%2;
        num /=2;
        for(usint j=0; j< num; j++ ){
            result[j]=cc->EvalMult(result[2*j],result[2*j+1]);
            cc->ModReduceInPlace(result[j]);
        }
        if(res==1){
            result[num]=result[2*num]->Clone();
            num+=res;
        }
    }
    return result[0];
}



Ciphertext<DCRTPoly> EvalLog(const Ciphertext<DCRTPoly> ciphertext, const double bound, const double base, const usint degree){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    Ciphertext<DCRTPoly> result = cc->EvalMult(ciphertext, 1/bound);
    cc->ModReduceInPlace(result);
    result = cc->EvalSub(result, 1); // y = x/bound -1

    vector<double> coeff(degree+1);
    double logbase = 1;
    if(base>1)logbase = log2(base);
    double coeff0 = 1/logbase;
    coeff[0]=0;
    for(usint i=1; i<degree+1; i++){
        coeff[i]=coeff0/(double)i;
        coeff0=-coeff0;
    }
    //log (1+y) +log(bound) = log x
    result = cc->EvalPoly(result, coeff);

    // Ciphertext<DCRTPoly> pert = EvalInverseAuto(ciphertext, bound);
    // pert = cc->EvalSub(pert, 1);
    // for(usint i=0; i< log2(degree); i++){
    //     cc->EvalMultInPlace(pert,pert);
    //     cc->ModReduceInPlace(pert);
    // }
    // cc->EvalMultInPlace(pert,log(bound)/logbase);
    // cc->ModReduceInPlace(pert);


    // result = cc->EvalAdd(result, log(bound)/logbase);
    // result = cc->EvalAdd(result, pert);

    return result;
}


Ciphertext<DCRTPoly> EvalLogLike(const Ciphertext<DCRTPoly> ciphertext, const double bound){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    Ciphertext<DCRTPoly> result = cc->EvalMult(ciphertext, 1/bound);
    cc->ModReduceInPlace(result);

    vector<double> coeff={0, 1, -0.5};
    result = cc->EvalPoly(result, coeff);

    return result;
}


Ciphertext<DCRTPoly> EvalInverse(const Ciphertext<DCRTPoly> ciphertext, const double bound, const usint degree){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    Ciphertext<DCRTPoly> sq = cc->EvalMult(ciphertext, 1/(2*bound)); 
    cc->ModReduceInPlace(sq);
    cc->EvalNegateInPlace(sq);
    Ciphertext<DCRTPoly> result = cc->EvalAdd(sq, 2); //1+y
    sq = cc->EvalAdd(sq, 1); // y = 1- x/2*bound.  [0, bound] -> [0, 0.5] -> [0.5, 1]

    for(usint i=1; i<degree; i++){
        sq = cc->EvalMult(sq, sq);
        cc->ModReduceInPlace(sq);
        Ciphertext<DCRTPoly> tmp=cc->EvalAdd(sq, 1);
        result = cc->EvalMult(result, tmp);
        cc->ModReduceInPlace(result);
    }

    result = cc->EvalMult(result, 0.5/bound);
    cc->ModReduceInPlace(result);

    return result;
}



vector<double> GetCoeff(const usint bound){
    vector<double> coeff(bound);

    for(usint i=0; i < bound; i++){
        coeff[i]=0;
	}
    coeff[0]=1;
    for(usint i=1; i < bound; i++){
        for(usint j=i; j!=0; j--){
            coeff[j]=coeff[j-1]-coeff[j]*((double)i/(double)bound);
            
        }
        coeff[0]*=(-(double)i/(double)bound);    
	}


    // for(usint i=1; i < bound; i++){
    //     coeff[i]/=coeff[0];
    // }
    // coeff[0]=1;
    return coeff;
}


Plaintext GenIndicatorChecker(const usint bound, const CryptoContext<DCRTPoly> cc){
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize();
    std::vector<double> num(bound);
    for(usint i=0 ; i < bound; i ++){
        num[i]=i;
    }

    std::vector<double> nums(batchSize);
    nums = fullCopy(num, batchSize, bound);
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(nums);
    return ptxt;
}

Plaintext GenIndicatorCheckerInterval(const usint from, const usint size, const CryptoContext<DCRTPoly> cc){
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize();

    std::vector<double> num(batchSize/size);
    for(usint i=0 ; i < batchSize/size; i ++){
        num[i]=from+i;
    }

    std::vector<double> nums(batchSize);
    nums = repeat(num, batchSize, size);
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(nums);
    return ptxt;
}


Plaintext GenIndicatorCheckerIntervalRecursive(const usint from, const usint to, const usint size, const CryptoContext<DCRTPoly> cc){
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize();

    std::vector<double> num(to-from);
    for(usint i=0 ; i < to-from; i ++){
        num[i]=from+i;
    }

    std::vector<double> nums(batchSize/size);
    nums = fullCopy(num, batchSize/size, to-from);

    std::vector<double> numss(batchSize);
    numss = repeat(nums, batchSize, size);
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(numss);
    return ptxt;
}

// Multi-version included
Plaintext GenIndicatorCheckerForSort(const int32_t arrsize, const CryptoContext<DCRTPoly> cc, const int32_t numk){
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize();
    int32_t numofrows = batchSize / arrsize;
    std::vector<double> num(numofrows);

    

    for(int32_t i=0 ; i < numofrows; i ++){
        num[i]=(double)(i % arrsize);
        if(num[i]>=numk && numk!=0){
            num[i]=arrsize; //padding
        }
    }
    
    std::vector<double> nums(batchSize);
    nums = repeat(num, batchSize, arrsize);
 
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(nums);
    return ptxt;
    
}



vector<Plaintext> GenIndicatorCheckerForSIMDCOUNT(const usint base, const usint size, const usint paral,const CryptoContext<DCRTPoly> cc){
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize();

    std::vector<double> num(base);
    for(usint i=0 ; i < base; i ++){
        num[i]=(double)i;
    }

    vector<Plaintext> ptxts(paral);

    for(usint i=0;i<paral;i++){
        usint pp=pow(base,i);
        std::vector<double> nums(pp*size*base);
        nums = repeat(num, pp*size*base, pp*size);

        std::vector<double> numss(batchSize);
        numss = fullCopy(nums, batchSize, pp*size*base);
        
        ptxts[i] = cc->MakeCKKSPackedPlaintext(numss);
    }

    return ptxts;
}


Plaintext GenIndicatorCheckerPartialArray(const usint from, const usint size, const CryptoContext<DCRTPoly> cc, const vector<usint> list){
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize();

    std::vector<double> num(batchSize/size);
    for(usint i=0 ; i < batchSize/size; i ++){
        num[i]=list[from+i];
    }

    std::vector<double> nums(batchSize);
    nums = repeat(num, batchSize, size);
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(nums);
    return ptxt;
}


vector<usint> GenIndicatorRounds(const usint bound, const usint scaleModSize){
	vector<usint> rounds(2);
    const usint boundbits=log2(bound);
    rounds[0]=2+boundbits*2;
    if(boundbits > 4)rounds[0]-=1;
    if(boundbits > 6)rounds[0]-=1;

    rounds[1]=1;
	if(boundbits > 3)rounds[1]+=1;
    if(boundbits > 5)rounds[1]+=1;

	if(scaleModSize >= 40){
        rounds[0]+=1;
        if(boundbits<=4)rounds[0]-=1;
        if(boundbits==7)rounds[0]+=1;
	}
    if(scaleModSize >= 45){
        if(boundbits<=4)rounds[0]+=1;
        if(boundbits==7)rounds[0]+=1;
        if(boundbits==9)rounds[0]+=1;
        if(boundbits==7)rounds[1]-=1;
        if(boundbits==6)rounds[1]-=1;
        if(boundbits==9)rounds[1]-=1;

	}
    if(scaleModSize >= 50){
        if(boundbits>6)rounds[0]+=1;
        if(boundbits==9)rounds[0]-=1;
        if(boundbits==13)rounds[0]-=1;
        if(boundbits==4)rounds[1]-=1;
        if(boundbits>6)rounds[1]-=1;
        if(boundbits==7){rounds[0]-=2; rounds[1]+=1;}
        if(boundbits==9)rounds[1]+=1;
        if(boundbits>=12)rounds[1]+=1;
	}
    if(scaleModSize >=55){
        if(boundbits>4 && boundbits < 9)rounds[0]+=1;
        if(boundbits>4 && boundbits < 9)rounds[1]-=1;

        if(boundbits>12)rounds[0]+=1;
        if(boundbits>11)rounds[1]-=1;


        // if(boundbits==12)rounds[1]-=1;
    }

    if(scaleModSize >= 59){
        if(boundbits==9)rounds[0]+=1;
        if(boundbits==9)rounds[1]-=1;

	}
    if(boundbits==1)rounds[0]=0;
    // if(boundbits==1 && scaleModSize >= 40)rounds[0]=1;
    if(boundbits==1 && scaleModSize >= 59)rounds[1]+=1;


    return rounds;
}


vector<usint> GenZeroTestRounds(const usint bound, const usint scaleModSize){
	vector<usint> rounds(2);
    const usint boundbits=log2(bound);
    rounds[0]=2+boundbits;

    rounds[1]=1;
	if(boundbits > 3)rounds[1]+=1;

	if(scaleModSize >= 40){
        rounds[0]+=1;
        if(boundbits==4)rounds[1]-=1;
        if(boundbits==5)rounds[1]-=1;
	}
    if(scaleModSize >= 45){
        if(boundbits==6)rounds[1]-=1;
        if(boundbits==7)rounds[1]-=1;
	}
    if(scaleModSize >= 50){
        if(boundbits>2)rounds[0]+=1;
        if(boundbits>6)rounds[0]-=1;

        if(boundbits==8)rounds[1]-=1;
        if(boundbits==9)rounds[1]-=1;
        if(boundbits==10)rounds[1]-=1;

	}
    if(scaleModSize >= 59){
        if(boundbits==2)rounds[0]+=1;
        if(boundbits>6)rounds[0]+=1;

        if(boundbits>10)rounds[1]-=1;
	}
    if(boundbits==1)rounds[0]=0;
    if(boundbits==1)rounds[1]=0;

    return rounds;
}

Ciphertext<DCRTPoly> RotAndSum(const Ciphertext<DCRTPoly> ciphertext, const int32_t from, const int32_t to, const bool verbose) {
    Ciphertext<DCRTPoly> result = ciphertext->Clone();
    Ciphertext<DCRTPoly> tmp; 
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();

    if(verbose)cout << "RotAndSum - from: " << from << ", to: " << to << endl;

    const int32_t round = log2(from/to);
	int32_t rotnum=from;
	for(int32_t s=1 ; s < round+1 ; s++){
        if(verbose)cout << "Rotating: " << (rotnum >> s) << endl;

        tmp = cc->EvalRotate(result, rotnum >> s);
        result= cc->EvalAdd(result, tmp);
	}

    return result;
}




Ciphertext<DCRTPoly> Cleanse(Ciphertext<DCRTPoly> ciphertext, const usint round) {
	const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    Ciphertext<DCRTPoly> result = ciphertext->Clone();

	for(usint i=0; i < round; i++){
        Ciphertext<DCRTPoly> power = cc->EvalSquare(result);
        cc->ModReduceInPlace(power);
        Ciphertext<DCRTPoly> power2 = cc->EvalAdd(power,power);

        //cc->LevelReduceInPlace(result, nullptr, 1);

        result = cc->EvalMult(result, power2);
        cc->ModReduceInPlace(result);

        cc->EvalAddInPlace(power,power2);
        //cc->LevelReduceInPlace(power, nullptr, 1);
        result = cc->EvalSub(power, result);
	}
    return result;
}

Ciphertext<DCRTPoly> Indicator(const Ciphertext<DCRTPoly> ciphertext, const usint bound, const vector<usint> rounds, const double numtocheck, const usint levlimit){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    const double div =1 / (double) bound;


    Ciphertext<DCRTPoly> result;
    if(bound > 2){
        result = cc->EvalSub(ciphertext, numtocheck);
        cc->EvalMultInPlace(result, div);
        cc->ModReduceInPlace(result);    
        result = cc->EvalSquare(result);
        cc->ModReduceInPlace(result);
        cc->EvalAddInPlace(result, result);
        cc->EvalAddInPlace(result, -1);
    }
    if(bound == 2 && numtocheck==1)result = ciphertext->Clone();
    if(bound == 2 && numtocheck==0){
        result = cc->EvalAdd(ciphertext, -1);
        cc->EvalNegateInPlace(result);
    }

	for(usint i=0; i<rounds[0]; i++){
        result = cc->EvalSquare(result);
        cc->ModReduceInPlace(result);
	}
	for(usint i=0; i < rounds[1]; i++){
        Ciphertext<DCRTPoly> power = cc->EvalSquare(result);
        cc->ModReduceInPlace(power);
        Ciphertext<DCRTPoly> power2 = cc->EvalAdd(power,power);

        //cc->LevelReduceInPlace(result, nullptr, 1);

        result = cc->EvalMult(result, power2);
        cc->ModReduceInPlace(result);

        cc->EvalAddInPlace(power,power2);
        //cc->LevelReduceInPlace(power, nullptr, 1);
        result = cc->EvalSub(power, result);
	}

    return result;
}

Ciphertext<DCRTPoly> IndicatorBinary(const Ciphertext<DCRTPoly> ciphertext, const vector<usint> rounds){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();

	
    Ciphertext<DCRTPoly> result = cc->EvalSquare(ciphertext);
    cc->ModReduceInPlace(result);
    cc->EvalAddInPlace(result, -1);
    cc->EvalNegateInPlace(result);
	for(usint i=0; i<rounds[0]; i++){
        result = cc->EvalSquare(result);
        cc->ModReduceInPlace(result);
	}
	for(usint i=0; i < rounds[1]; i++){
        Ciphertext<DCRTPoly> power = cc->EvalSquare(result);
        cc->ModReduceInPlace(power);
        Ciphertext<DCRTPoly> power2 = cc->EvalAdd(power,power);

        //cc->LevelReduceInPlace(result, nullptr, 1);

        result = cc->EvalMult(result, power2);
        cc->ModReduceInPlace(result);

        cc->EvalAddInPlace(power,power2);
        //cc->LevelReduceInPlace(power, nullptr, 1);
        result = cc->EvalSub(power, result);
	}

    return result;
}


Ciphertext<DCRTPoly> ZeroTest(const Ciphertext<DCRTPoly> ciphertext, const usint bound, const vector<usint> rounds){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    const double div =1 / (double) bound;

	
    Ciphertext<DCRTPoly> result = ciphertext->Clone();
    if(bound!=1.0){
        result = cc->EvalMult(result, div);
        cc->ModReduceInPlace(result);
    }
    cc->EvalAddInPlace(result, result);
    cc->EvalAddInPlace(result, -1);

	for(usint i=0; i<rounds[0]; i++){
        result = cc->EvalSquare(result);
        cc->ModReduceInPlace(result);
	}
	for(usint i=0; i < rounds[1]; i++){
        Ciphertext<DCRTPoly> power = cc->EvalSquare(result);
        cc->ModReduceInPlace(power);
        Ciphertext<DCRTPoly> power2 = cc->EvalAdd(power,power);

        //cc->LevelReduceInPlace(result, nullptr, 1);

        result = cc->EvalMult(result, power2);
        cc->ModReduceInPlace(result);

        cc->EvalAddInPlace(power,power2);
        //cc->LevelReduceInPlace(power, nullptr, 1);
        result = cc->EvalSub(power, result);
	}

    return result;
}

Ciphertext<DCRTPoly> IndicatorSIMD(const Ciphertext<DCRTPoly> ciphertext, const usint bound, const vector<usint> rounds, const Plaintext numtocheck, const usint levlimit){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    // usint budget = cc->GetMultiplicativeDepth();
    const double div =1 / (double) bound;
	
    Ciphertext<DCRTPoly> result;
    
    if(bound > 2){
        result = cc->EvalSub(ciphertext, numtocheck);
        cc->EvalMultInPlace(result, div); 
        cc->ModReduceInPlace(result);    
        result = cc->EvalSquare(result);
        cc->ModReduceInPlace(result);
        cc->EvalAddInPlace(result, result);
        cc->EvalAddInPlace(result, -1);
    }

    if(bound == 2){
        vector<double> checker = numtocheck->GetRealPackedValue();
        vector<double> checkeradd(checker.size());
        vector<double> checkermult(checker.size());
        for(usint i=0; i<checker.size(); i++){
            if(checker[i]>0.5){
                checkeradd[i]=0.0;
                checkermult[i]=1.0;
            }else{
                checkeradd[i]=1.0;
                checkermult[i]=-1.0;
            }
        }
        Plaintext ptxtadd = cc->MakeCKKSPackedPlaintext(checkeradd);
        Plaintext ptxtmult = cc->MakeCKKSPackedPlaintext(checkermult);
        result = cc->EvalMult(ciphertext, ptxtmult);
        result = cc->EvalAdd(result, ptxtadd);
        // cc->EvalNegateInPlace(result);
    }

    //     result = cc->EvalSquare(result);
    //     cc->ModReduceInPlace(result);
    //     cc->EvalAddInPlace(result, -1);
    //     cc->EvalNegateInPlace(result);
    // }

	for(usint i=0; i<rounds[0]; i++){
        if(levlimit < (result->GetLevel()))result = cc->EvalBootstrap(result,2,10);
        result = cc->EvalSquare(result);
        cc->ModReduceInPlace(result);
	}
	for(usint i=0; i < rounds[1]; i++){
        if(levlimit < (result->GetLevel()))result = cc->EvalBootstrap(result,2,10);
        Ciphertext<DCRTPoly> power = cc->EvalSquare(result);
        cc->ModReduceInPlace(power);
        Ciphertext<DCRTPoly> power2 = cc->EvalAdd(power,power);

        //cc->LevelReduceInPlace(result, nullptr, 1);

        result = cc->EvalMult(result, power2);
        cc->ModReduceInPlace(result);

        cc->EvalAddInPlace(power,power2);
        //cc->LevelReduceInPlace(power, nullptr, 1);
        result = cc->EvalSub(power, result);
	}

    return result;
}




//----------------------------------------------------------------------------------
//   Comparison
//----------------------------------------------------------------------------------

Ciphertext<DCRTPoly> normalize(const Ciphertext<DCRTPoly> ciphertext, double bound){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    
    const double div =1 / (double) bound;
    Ciphertext<DCRTPoly> result = cc->EvalMult(ciphertext, div);
    cc->ModReduceInPlace(result);    
    
    return result;
}

// Ciphertext<DCRTPoly> comparison(const Ciphertext<DCRTPoly> ciphertext, const usint degf, const usint degg, double bound, const usint ver){
//     const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
//     auto result = ciphertext->Clone();

//     std::vector<double> fcoeff;
//     std::vector<double> gcoeff;

//     if(ver==1){
//         fcoeff= {0, 1.5, 0, -0.5};
//         gcoeff= {0, 2.076171875, 0, -1.3271484375};        
//     }

//     if(ver==2){
//         fcoeff = {0, 1.875, 0, -1.25, 0, 0.375};
//         gcoeff = {0, 3.255859375, 0, -5.96484375, 0, 3.70703125};
//     }

//     if(ver==3){
//         fcoeff = {0, 2.1875, 0, -2.1875, 0, 1.3125, 0, -0.3125};
//         gcoeff = {0, 4.4814453125, 0, -16.1884765625, 0, 25.013671875, 0, -12.55859375};
//     }

//     if(ver==4){
//         fcoeff= {0, 2.4609375, 0, -3.28125, 0, 2.953125, 0, -1.40625, 0 , 0.2734375};
//         gcoeff= {0, 5.712890625, 0, -34.154296875, 0, 94.7412109375, 0, -110.83203125, 0 , 45.5302734375};
//     }
//     for(usint i=0;i<degg;i++){
//         result = cc->EvalPoly(result, gcoeff);
//     }
//     for(usint i=0;i<degf;i++){
//         result = cc->EvalPoly(result, fcoeff);
//     }
//     return result;
// }

Ciphertext<DCRTPoly> comparison(const Ciphertext<DCRTPoly> ciphertext, const usint degf, const usint degg, double bound, const usint ver){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    Ciphertext<DCRTPoly> result = ciphertext->Clone();

    std::vector<double> fcoeff;
    std::vector<double> gcoeff;

    if(ver==1){
        fcoeff= {0, 1.5, 0, -0.5};
        gcoeff= {0, 2.076171875, 0, -1.3271484375};        
    }

    if(ver==2){
        fcoeff = {0, 1.875, 0, -1.25, 0, 0.375};
        gcoeff = {0, 3.255859375, 0, -5.96484375, 0, 3.70703125};
    }

    if(ver==3){
        fcoeff = {0, 2.1875, 0, -2.1875, 0, 1.3125, 0, -0.3125};
        gcoeff = {0, 4.4814453125, 0, -16.1884765625, 0, 25.013671875, 0, -12.55859375};
    }

    if(ver==4){
        fcoeff= {0, 2.4609375, 0, -3.28125, 0, 2.953125, 0, -1.40625, 0 , 0.2734375};
        gcoeff= {0, 5.712890625, 0, -34.154296875, 0, 94.7412109375, 0, -110.83203125, 0 , 45.5302734375};
    }
    for(usint i=0;i<degg;i++){
        vector<Ciphertext<DCRTPoly>> powers(3);
        powers[0]=cc->EvalMult(result, result);
        cc->ModReduceInPlace(powers[0]);
        powers[1]=cc->EvalMult(powers[0], powers[0]);
        cc->ModReduceInPlace(powers[1]);
        powers[0]=cc->EvalMult(powers[0], result);
        cc->ModReduceInPlace(powers[0]);
        powers[2]=cc->EvalMult(powers[1], powers[0]);
        cc->ModReduceInPlace(powers[2]);
        powers[1]=cc->EvalMult(powers[1], result);
        cc->ModReduceInPlace(powers[1]);

        result=cc->EvalMult(result, 4.4814453125);
        powers[0]=cc->EvalMult(powers[0], -16.1884765625);
        powers[1]=cc->EvalMult(powers[1], 25.013671875);
        powers[2]=cc->EvalMult(powers[2], -12.55859375);
        result = cc-> EvalAdd(result,powers[0]);
        result = cc-> EvalAdd(result,powers[1]);
        result = cc-> EvalAdd(result,powers[2]);
        cc->ModReduceInPlace(result);
    }
    for(usint i=0;i<degf;i++){
        vector<Ciphertext<DCRTPoly>> powers(3);
        powers[0]=cc->EvalMult(result, result);
        cc->ModReduceInPlace(powers[0]);
        powers[1]=cc->EvalMult(powers[0], powers[0]);
        cc->ModReduceInPlace(powers[1]);
        powers[0]=cc->EvalMult(powers[0], result);
        cc->ModReduceInPlace(powers[0]);
        powers[2]=cc->EvalMult(powers[1], powers[0]);
        cc->ModReduceInPlace(powers[2]);
        powers[1]=cc->EvalMult(powers[1], result);
        cc->ModReduceInPlace(powers[1]);

        result=cc->EvalMult(result, 2.1875);
        powers[0]=cc->EvalMult(powers[0], -2.1875);
        powers[1]=cc->EvalMult(powers[1], 1.3125);
        powers[2]=cc->EvalMult(powers[2], -0.3125);
        result = cc-> EvalAdd(result,powers[0]);
        result = cc-> EvalAdd(result,powers[1]);
        result = cc-> EvalAdd(result,powers[2]);
        cc->ModReduceInPlace(result);    
    }
    return result;
}

Ciphertext<DCRTPoly> comp(const Ciphertext<DCRTPoly> ciphertext, const double bound, const bool lastmod, const usint levlimit){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    Ciphertext<DCRTPoly> result = ciphertext->Clone();
    
    usint degg=2;
    usint logbound=log2(bound);
    // if(logbound>5)degg+=1;
    // if(logbound>7)degg+=1;
    // if(logbound>9)degg+=1;
    if(logbound>3)degg=logbound/2;



    usint degf=2;

    vector<Ciphertext<DCRTPoly>> powers(3);


    for(usint i=0;i<degg;i++){
        if(levlimit < (result->GetLevel()))result = cc->EvalBootstrap(result,2,10);

        powers[0]=cc->EvalMult(result, result);
        cc->ModReduceInPlace(powers[0]);
        powers[1]=cc->EvalMult(powers[0], powers[0]);
        cc->ModReduceInPlace(powers[1]);
        powers[0]=cc->EvalMult(powers[0], result);
        cc->ModReduceInPlace(powers[0]);
        powers[2]=cc->EvalMult(powers[1], powers[0]);
        cc->ModReduceInPlace(powers[2]);
        powers[1]=cc->EvalMult(powers[1], result);
        cc->ModReduceInPlace(powers[1]);

        result=cc->EvalMult(result, 4.4814453125);
        powers[0]=cc->EvalMult(powers[0], -16.1884765625);
        powers[1]=cc->EvalMult(powers[1], 25.013671875);
        powers[2]=cc->EvalMult(powers[2], -12.55859375);
        result = cc-> EvalAdd(result,powers[0]);
        result = cc-> EvalAdd(result,powers[1]);
        result = cc-> EvalAdd(result,powers[2]);
        cc->ModReduceInPlace(result);
    }

    // if(levlimit < (result->GetLevel()))result = cc->EvalBootstrap(result,2,10);

    for(usint i=0;i<degf;i++){
        if(levlimit < (result->GetLevel()))result = cc->EvalBootstrap(result,2,10);

        powers[0]=cc->EvalMult(result, result);
        cc->ModReduceInPlace(powers[0]);
        powers[1]=cc->EvalMult(powers[0], powers[0]);
        cc->ModReduceInPlace(powers[1]);
        powers[0]=cc->EvalMult(powers[0], result);
        cc->ModReduceInPlace(powers[0]);
        powers[2]=cc->EvalMult(powers[1], powers[0]);
        cc->ModReduceInPlace(powers[2]);
        powers[1]=cc->EvalMult(powers[1], result);
        cc->ModReduceInPlace(powers[1]);

        if(i==degf-1 && lastmod == true){
            result=cc->EvalMult(result, 1.09375);
            powers[0]=cc->EvalMult(powers[0], -1.09375);
            powers[1]=cc->EvalMult(powers[1], 0.65625);
            powers[2]=cc->EvalMult(powers[2], -0.15625);
            result = cc-> EvalAdd(result,powers[0]);
            result = cc-> EvalAdd(result,powers[1]);
            result = cc-> EvalAdd(result,powers[2]);
            cc->ModReduceInPlace(result);
            result = cc-> EvalAdd(result, 0.5);
        }else{
            result=cc->EvalMult(result, 2.1875);
            powers[0]=cc->EvalMult(powers[0], -2.1875);
            powers[1]=cc->EvalMult(powers[1], 1.3125);
            powers[2]=cc->EvalMult(powers[2], -0.3125);
            result = cc-> EvalAdd(result,powers[0]);
            result = cc-> EvalAdd(result,powers[1]);
            result = cc-> EvalAdd(result,powers[2]);
            cc->ModReduceInPlace(result);
        }

    }


    return result;
}

Ciphertext<DCRTPoly> compandUp(const Ciphertext<DCRTPoly> ciphertext, const double bound, const bool boot,  const usint up){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    Ciphertext<DCRTPoly> result = ciphertext->Clone();
    
    usint degg=2;
    usint logbound=log2(bound);
    if(logbound>5)degg+=1;
    if(logbound>7)degg+=1;
    if(logbound>9)degg+=1;

    usint degf=2;

    vector<Ciphertext<DCRTPoly>> powers(3);


    for(usint i=0;i<degg;i++){
        powers[0]=cc->EvalMult(result, result);
        cc->ModReduceInPlace(powers[0]);
        powers[1]=cc->EvalMult(powers[0], powers[0]);
        cc->ModReduceInPlace(powers[1]);
        powers[0]=cc->EvalMult(powers[0], result);
        cc->ModReduceInPlace(powers[0]);
        powers[2]=cc->EvalMult(powers[1], powers[0]);
        cc->ModReduceInPlace(powers[2]);
        powers[1]=cc->EvalMult(powers[1], result);
        cc->ModReduceInPlace(powers[1]);

        result=cc->EvalMult(result, 4.4814453125);
        powers[0]=cc->EvalMult(powers[0], -16.1884765625);
        powers[1]=cc->EvalMult(powers[1], 25.013671875);
        powers[2]=cc->EvalMult(powers[2], -12.55859375);
        result = cc-> EvalAdd(result,powers[0]);
        result = cc-> EvalAdd(result,powers[1]);
        result = cc-> EvalAdd(result,powers[2]);
        cc->ModReduceInPlace(result);
    }

    if(boot)result = cc->EvalBootstrap(result);

    for(usint i=0;i<degf;i++){
        powers[0]=cc->EvalMult(result, result);
        cc->ModReduceInPlace(powers[0]);
        powers[1]=cc->EvalMult(powers[0], powers[0]);
        cc->ModReduceInPlace(powers[1]);
        powers[0]=cc->EvalMult(powers[0], result);
        cc->ModReduceInPlace(powers[0]);
        powers[2]=cc->EvalMult(powers[1], powers[0]);
        cc->ModReduceInPlace(powers[2]);
        powers[1]=cc->EvalMult(powers[1], result);
        cc->ModReduceInPlace(powers[1]);

        if(i==degf-1 && up != 1){
            double upscale = (double)(1 << up);
            result=cc->EvalMult(result, 1.09375*upscale);
            powers[0]=cc->EvalMult(powers[0], -1.09375*upscale);
            powers[1]=cc->EvalMult(powers[1], 0.65625*upscale);
            powers[2]=cc->EvalMult(powers[2], -0.15625*upscale);
            result = cc-> EvalAdd(result,powers[0]);
            result = cc-> EvalAdd(result,powers[1]);
            result = cc-> EvalAdd(result,powers[2]);
            cc->ModReduceInPlace(result);
            result = cc-> EvalAdd(result, 0.5*upscale);
        }else{
            result=cc->EvalMult(result, 2.1875);
            powers[0]=cc->EvalMult(powers[0], -2.1875);
            powers[1]=cc->EvalMult(powers[1], 1.3125);
            powers[2]=cc->EvalMult(powers[2], -0.3125);
            result = cc-> EvalAdd(result,powers[0]);
            result = cc-> EvalAdd(result,powers[1]);
            result = cc-> EvalAdd(result,powers[2]);
            cc->ModReduceInPlace(result);
        }

    }


    return result;
}

Ciphertext<DCRTPoly> compDecrete(const Ciphertext<DCRTPoly> ciphertext, const int32_t bound){
    const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    Ciphertext<DCRTPoly> input=cc->EvalAdd(ciphertext, bound);
    Ciphertext<DCRTPoly> tmp;
    uint32_t scalingfactor = cc->GetEncodingParams()->GetPlaintextModulus();

    vector<usint> rounds=GenIndicatorRounds(2*bound,scalingfactor);
    Ciphertext<DCRTPoly> result=Indicator(input, 2*bound, rounds, (double)bound);
    result = cc->EvalMult(result,0.5);
    cc->ModReduceInPlace(result);

    for(int32_t j=bound+1; j< 2*bound-1; j++){
	    tmp = Indicator(input, 2*bound, rounds, (double)j);        
        result = cc->EvalAdd(result, tmp);
    }
	
    return result;
}

Ciphertext<DCRTPoly> fakeboot(const Ciphertext<DCRTPoly> ciphertext, const CryptoContext<DCRTPoly> cc, const KeyPair<DCRTPoly> keys, bool verbose){
    Plaintext result;
    usint depth = 19;
    if(verbose)cout << "Before: " << ciphertext->GetLevel() << endl;

    cc->Decrypt(keys.secretKey, ciphertext, &result);
    vector<double> x1 = result->GetRealPackedValue();
    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x1, 1, depth);
    Ciphertext<DCRTPoly> c1 = cc->Encrypt(keys.publicKey, ptxt1);

    // usint count = 65536;
    // for(usint i=0; i< 65536; i++){
    //     if(x1[i] < 0.00001)count-=1;
    // }
    // cout << "number of nonzeros: " << count << endl;


    if(verbose)cout << "After: " << c1->GetLevel() << endl;

    return c1;
}

Ciphertext<DCRTPoly> boot(const Ciphertext<DCRTPoly> ciphertext, const CryptoContext<DCRTPoly> cc, const KeyPair<DCRTPoly> keys){
    Plaintext result;
    cc->Decrypt(keys.secretKey, ciphertext, &result);
    vector<double> x1 = result->GetRealPackedValue();
    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x1);
    Ciphertext<DCRTPoly> c1 = cc->Encrypt(keys.publicKey, ptxt1);
    return c1;
}






vector<Plaintext> maskPrecompute(const int32_t arrsize, const int32_t batchSize, const int32_t skewnormalizer, const CryptoContext<DCRTPoly> cc, const bool skew, const int32_t numk){
    vector<Plaintext> result(5);
	vector<double> masking(batchSize);
    int32_t arrsizesquare = arrsize*arrsize;

    //vertical mask : Extract Column
	for(int32_t s=0; s<batchSize; s++){
		masking[s]=0.0;
	}
	for(int32_t i=0; i<(batchSize/arrsizesquare);i++){
		int32_t interval=i*(arrsizesquare);
		for(int32_t s=0; s<arrsize; s++){
			masking[interval+s*arrsize]=1.0;
		}
	}
	
    result[0] = cc->MakeCKKSPackedPlaintext(masking);

    //horizontal mask : Extract Row
	for(int32_t s=0; s<batchSize; s++){
		masking[s]=0.0;
	}
	for(int32_t i=0; i<(batchSize/arrsizesquare);i++){
		int32_t interval=i*(arrsizesquare);
		for(int32_t s=0; s<arrsize; s++){
			masking[interval+s]=1.0;
		}
	}
    
    result[1] = cc->MakeCKKSPackedPlaintext(masking);

    if(skew){
        //skew mask
        for(int32_t s=0; s<batchSize; s++){
            masking[s]=1/(double)skewnormalizer;
        }
        for(int32_t i=0; i<(batchSize/arrsizesquare);i++){
            int32_t interval=i*(arrsizesquare);
            // for(int32_t s=0; s<size; s++){
            // 	masking[interval+s*size+s]=-0.5; //diagonal entries
            // }
            for(int32_t s=0; s<arrsize; s++){
                for(int32_t t=0; t<s+1; t++){
                    masking[interval+s*arrsize+t]=-1/(double)skewnormalizer; //lowertriangle and diagonal
                }
            }
        }

        result[2] = cc->MakeCKKSPackedPlaintext(masking);
    }else{

        int32_t sublength = numk;
        if(numk==0)sublength = arrsize / 2;
        //vertical mask for topk : Extract Column
        for(int32_t s=0; s<batchSize; s++){
            masking[s]=0.0;
        }
        for(int32_t i=0; i<(batchSize/arrsizesquare);i++){
            int32_t interval=i*(arrsizesquare);
            for(int32_t s=0; s < sublength ; s++){
                masking[interval+s*arrsize]=1.0;
            }
        }
        
        result[2] = cc->MakeCKKSPackedPlaintext(masking);

    }


    //Bootstrapping mask - horizontal : Extract Row, divide by arrsize
	for(int32_t s=0; s<batchSize; s++){
		masking[s]=0.0;
	}
	for(int32_t i=0; i<(batchSize/arrsizesquare);i++){
		int32_t interval=i*(arrsizesquare);
		for(int32_t s=0; s<arrsize; s++){
			masking[interval+s]=1.0 / ((double)arrsize);
		}
	}
    
    result[3] = cc->MakeCKKSPackedPlaintext(masking);


    //Bootstrapping mask - horizontal : Extract Row, mult by arrsize
	for(int32_t s=0; s<batchSize; s++){
		masking[s]=0.0;
	}
	for(int32_t i=0; i<(batchSize/arrsizesquare);i++){
		int32_t interval=i*(arrsizesquare);
		for(int32_t s=0; s<arrsize; s++){
			masking[interval+s]=(double)arrsize;
		}
	}
    
    result[4] = cc->MakeCKKSPackedPlaintext(masking);


    return result;
}


vector<Plaintext> maskPrecompute_full(const int32_t size, const int32_t batchSize, const int32_t skewnormalizer,  const CryptoContext<DCRTPoly> cc){
	vector<double> masking(batchSize);
    int32_t sizesquare = size*size;
    int32_t numct = sizesquare / batchSize;
    int32_t slicelength = batchSize/size;
    if(slicelength > size)slicelength=size;
    vector<Plaintext> result(2+numct);

    //vertical mask : Extract Column
	for(int32_t s=0; s<batchSize; s++){
		masking[s]=0.0;
	}
    for(int32_t s=0; s<slicelength; s++){
		masking[s*size]=1.0;
	}
	for(int32_t i=1; i<(batchSize/sizesquare);i++){
		int32_t interval=i*(sizesquare);
		for(int32_t s=0; s<size; s++){
			masking[interval+s*size]=1.0;
		}
	}
	
    result[0] = cc->MakeCKKSPackedPlaintext(masking);

    //horizontal mask : Extract Row
	for(int32_t s=0; s<batchSize; s++){
		masking[s]=0.0;
	}
    for(int32_t s=0; s<size; s++){
		masking[s]=1.0;
	}
	for(int32_t i=1; i<(batchSize/sizesquare);i++){
		int32_t interval=i*(sizesquare);
		for(int32_t s=0; s<size; s++){
			masking[interval+s]=1.0;
		}
	}
    
    result[1] = cc->MakeCKKSPackedPlaintext(masking);

    //skew mask


    for(int32_t n=0; n< numct; n++){
        for(int32_t s=0; s<batchSize; s++){
            masking[s]=1/(double)skewnormalizer;
        }
        int32_t initpt = n*slicelength;
        // for(int32_t s=0; s<slicelength; s++){
        //     masking[initpt+s*size+s]=-1.0/(double)skewnormalizer; //diagonal entries
        // }
        for(int32_t s=0; s<slicelength; s++){
            for(int32_t t=0; t<s+initpt+1; t++){
                masking[s*size+t]=-1/(double)skewnormalizer; //lowertriangle and diagonal
            }
        }
        result[2+n] = cc->MakeCKKSPackedPlaintext(masking);
    }
	


    return result;
}




Ciphertext<DCRTPoly> sort(const Ciphertext<DCRTPoly> ciphertext, const int32_t arrsize, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys, const bool boot1, const bool boot2, const Plaintext targetranks, const int32_t padding, const bool verbose){
	const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize(); 

    int32_t arrsizesquare = arrsize*arrsize;

    TimeVar t;

    TIC(t);
    vector<Plaintext> masks = maskPrecompute(arrsize, batchSize, 2*bound, cc);

    Ciphertext<DCRTPoly> tmp, copy;
    Plaintext pt;
    int32_t bootthreshold = 34;
    if(arrsize == 2)bootthreshold = 33;


    int32_t cipherlevel = ciphertext->GetLevel();

    if(verbose)cout <<"input level " << cipherlevel << ", arrsize " << arrsize << endl;
    copy = ciphertext->Clone();

    if(boot1){
        if(cipherlevel > bootthreshold){
            if(verbose)cout << "Boot: " <<  copy->GetLevel() << endl;
            copy = cc->EvalBootstrap(copy,2,10);
        }
    }else{
        if(boot2){
            copy = fakeboot(copy, cc, keys, verbose); // if do not want to track the depth consumption, set boot2 = true.
        }else{
            if(cipherlevel > bootthreshold){
                copy = fakeboot(copy, cc, keys, verbose);
            }
        }
        
    }


    tmp = RotAndSum(copy, -(arrsizesquare-arrsize), -(arrsize-1)); //transpose row to col
    tmp = cc->EvalMult(tmp, masks[0]);
    cc->ModReduceInPlace(tmp);
    tmp =  RotAndSum(tmp, -arrsize, -1); //copy col = ct_cc
    copy =  RotAndSum(copy, -arrsizesquare, -arrsize); // copy row = ct_rc
    tmp = cc->EvalSub(tmp, copy); // ct_cc-ct_rc
    tmp = cc->EvalAdd(tmp, tmp);
    tmp = cc->EvalAdd(tmp, masks[2]);

    tmp = comp(tmp, 2*bound, true, 45); //obtain M_comp

 
    cipherlevel = tmp->GetLevel();

    if(verbose)cout << "After comp: " << cipherlevel << endl;
    
    if(boot1){
        if(cipherlevel > bootthreshold){
            if(verbose)cout << "Boot: " <<  tmp->GetLevel() << endl;
            tmp = cc->EvalBootstrap(tmp,2,10);
        }
    }else{
        if(boot2){
            tmp = fakeboot(tmp, cc, keys, verbose); // if do not want to track the depth consumption, set boot2 = true.
        }else{
            if(cipherlevel > bootthreshold){
                tmp = fakeboot(tmp, cc, keys, verbose);
            }
        }
        
    }


    if(arrsizesquare==batchSize){
        tmp =  RotAndSum(tmp, arrsizesquare, arrsize); //row sum & copy row
    }else{
        tmp =  RotAndSum(tmp, arrsizesquare, arrsize); //row sum, order vector

        tmp = cc->EvalMult(tmp, masks[1]);
        cc->ModReduceInPlace(tmp); 

        tmp = RotAndSum(tmp, -arrsizesquare, -arrsize); //copy row
    }
    
    vector<usint> rounds = GenIndicatorRounds(max(arrsize,4), scaleModSize);   
    tmp = IndicatorSIMD(tmp, max(arrsize+padding,4), rounds, targetranks, 49);  // if arrsize = 2, use indicator for 4 instead.


    tmp = cc->EvalMult(tmp, copy);
    cc->ModReduceInPlace(tmp);
    tmp = RotAndSum(tmp, arrsize, 1); // col sum

    if(verbose)cout << "Time " << TOC(t) <<", Comp: " << cipherlevel << ", Ind: " << tmp->GetLevel() <<  endl;

    return tmp;
}





vector<Ciphertext<DCRTPoly>> sort_SIMD(const vector<Ciphertext<DCRTPoly>> ciphertext, const int32_t arrsize, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys, const int32_t numct, const bool boot1, const Plaintext targetranks, const int32_t padding, const bool multithread, const bool verbose){
	const CryptoContext<DCRTPoly> cc = ciphertext[0]->GetCryptoContext();
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize(); 

    int32_t arrsizesquare = arrsize*arrsize;

    TimeVar t;

    TIC(t);
    vector<Plaintext> masks = maskPrecompute(arrsize, batchSize, 2*bound, cc);

    Plaintext pt;
    // int32_t bootthreshold = 34;
    // if(arrsize == 2)bootthreshold = 33;


    int32_t cipherlevel = ciphertext[0]->GetLevel();
    if(verbose)cout <<"input level " << cipherlevel << ", arrsize " << arrsize << endl;
    bool compboot = false;
    if(cipherlevel > 5)compboot = true; 

    vector<Ciphertext<DCRTPoly>> copy(numct);
    vector<Ciphertext<DCRTPoly>> sorted(numct);

    int32_t numpacking = max(numct / arrsize, 1); //number of packing ciphertexts for bootstrapping
    int32_t numpacked = numct/numpacking; //number of ciphertexts packed in a intermediate ciphertext

    #pragma omp parallel for if(multithread)
    for(int32_t i = 0; i < numct; i++){
        // TimeVar tt;    
        // TIC(tt);

        // cout << "Sort part 1 " << i << endl;
        copy[i] = ciphertext[i]->Clone();
        sorted[i] = RotAndSum(copy[i], -(arrsizesquare-arrsize), -(arrsize-1)); //transpose row to col
        sorted[i] = cc->EvalMult(sorted[i], masks[0]);
        cc->ModReduceInPlace(sorted[i]);
        sorted[i] =  RotAndSum(sorted[i], -arrsize, -1); //copy col = ct_cc
        copy[i] =  RotAndSum(copy[i], -arrsizesquare, -arrsize); // copy row = ct_rc
        sorted[i] = cc->EvalSub(sorted[i], copy[i]); // ct_cc-ct_rc
        sorted[i] = cc->EvalAdd(sorted[i], sorted[i]);
        sorted[i] = cc->EvalAdd(sorted[i], masks[2]);

        sorted[i] = comp(sorted[i], 2*bound, true, 45); //obtain M_comp

        sorted[i] =  RotAndSum(sorted[i], arrsizesquare, arrsize); //row sum, rank vector

        int32_t idx1 = i % numpacked;

        if(compboot){
            sorted[i] = cc->EvalMult(sorted[i], masks[3]); //row mask for bootstrapping, embed to [0,1]
        }else{
            sorted[i] = cc->EvalMult(sorted[i], masks[1]); //row mask
        }

        if(idx1!=0 && compboot)sorted[i]=cc->EvalRotate(sorted[i], -idx1*arrsize); 

        // cout << "Time " << i << " : " << TOC(tt) << " ,,,, " << TOC(t) << endl;
    }

    cipherlevel = sorted[0]->GetLevel();
    if(verbose)cout << "After comp: " << cipherlevel << endl;

    vector<Ciphertext<DCRTPoly>> gathered(numpacking);


    if(compboot){
        cout << "The number of Bootstrapping ciphertexts: " << numpacking << endl;

        // cc->Decrypt(keys.secretKey, sorted[0], &pt);
        // std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
        // printpt(pt, arrsize * 2, arrsize, 1, false, true);
        // printpt(pt, arrsize, 0, arrsize, false, true);

        // Aggregate
        for(int32_t i=0;i< log2(numpacked);i++){
            int32_t interval = 1 << i;
            #pragma omp parallel for if(multithread)
            for(int32_t j=0;j < numct / (2*interval) ; j++){
                // cout << j*interval*2 << ", " << j*interval*2 + interval << endl;
                sorted[j*interval*2]=cc->EvalAdd(sorted[j*interval*2], sorted[j*interval*2 + interval]);
            }
        }

        // Boot
        // #pragma omp parallel for if(multithread && (numpacking > 2))
        for(int32_t i=0;i< numpacking; i++){
            cc->ModReduceInPlace(sorted[i*numpacked]); 
            if(boot1){
                if(verbose && i==0)cout << "Boot: " <<  sorted[i*numpacked]->GetLevel() << endl;
                gathered[i] = cc->EvalBootstrap(sorted[i*numpacked],2,10);
            }else{
                gathered[i] = fakeboot(sorted[i*numpacked], cc, keys, verbose);
            }
        }

        // Divide
        #pragma omp parallel for if(multithread)
        for(int32_t i=0;i< numct;i++){
            int32_t idx0 = i / numpacked;
            int32_t idx1 = i % numpacked;
            sorted[i] = gathered[idx0]->Clone();
            if(idx1 != 0)sorted[i] = cc->EvalRotate(sorted[i], idx1*arrsize);
            sorted[i] = cc->EvalMult(sorted[i], masks[4]); // rowmask for post-bootstrapping, un-embed from [0,1]
            cc->ModReduceInPlace(sorted[i]); 
        }

        // cc->Decrypt(keys.secretKey, sorted[0], &pt);
        // std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
        // printpt(pt, arrsize * 2, arrsize, 1, false, true);
        // printpt(pt, arrsize, 0, arrsize, false, true);


    }   
    
    vector<usint> rounds = GenIndicatorRounds(max(arrsize,4), scaleModSize);   
   
    #pragma omp parallel for if(multithread)
    for(int32_t i = 0; i < numct; i++){
        // TimeVar tt;    
        // TIC(tt);
        sorted[i] = RotAndSum(sorted[i], -arrsizesquare, -arrsize); //copy row
    

        sorted[i] = IndicatorSIMD(sorted[i], max(arrsize+padding,4), rounds, targetranks); // if arrsize = 2, use indicator for 4 instead.


        sorted[i] = cc->EvalMult(sorted[i], copy[i]);
        cc->ModReduceInPlace(sorted[i]);
        sorted[i] = RotAndSum(sorted[i], arrsize, 1); // col sum
        
        // cout << "Time " << i << " : " << TOC(tt) << " ,,,, " << TOC(t)<<  endl;

    }



    if(verbose)cout << "Time " << TOC(t) <<", Comp: " << cipherlevel << ", Ind: " << sorted[0]->GetLevel() <<  endl;

    // cc->Decrypt(keys.secretKey, sorted[0], &pt);
    // std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
    // printpt(pt, arrsize * 2, arrsize, 1, false, true);
	// printpt(pt, arrsize, 0, arrsize, false, true);


    return sorted;
}

vector<Ciphertext<DCRTPoly>> rearrange_SIMD(const vector<Ciphertext<DCRTPoly>> ciphertext, const int32_t arrsize, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys,  const int32_t numk, const int32_t numct, const bool boot1, const bool multithread, const bool verbose){
    const CryptoContext<DCRTPoly> cc = ciphertext[0]->GetCryptoContext();
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize(); 

    int32_t arrsizesquare = arrsize*arrsize;
    
    vector<Plaintext> masks = maskPrecompute(arrsize, batchSize, 1, cc, false, numk); // 0: colmask 1: rowmask 2: topkmask

    int32_t numpacking = max(numct * numk / arrsizesquare, 1); //number of packing ciphertexts for bootstrapping
    int32_t numpacked = numct / numpacking; //number of ciphertexts packed in an intermediate ciphertext
    int32_t var1 = arrsize /  min(arrsize, numct); //required number of blocks to fill a row
    int32_t var2 = arrsize /  (numk * max(numpacked / arrsize, 1)); //required number of blocks to fill columns

    vector<Ciphertext<DCRTPoly>> gathered(numpacking);
    vector<Ciphertext<DCRTPoly>> sorted(numct);

    int32_t next_arrsize = min(arrsizesquare / numk, 256);

    int32_t resultnumct = numct * numk * next_arrsize / arrsizesquare;
    int32_t resultnumpacked = resultnumct / numpacking; //number of result ciphertexts to be divdided per am intermediate ciphertext

    cout << "Merge and divide, num of ciphertext: " << numct << ", number of Bootstrapping ciphertexts: " << numpacking << ", current sub-array size: " << arrsize << ", next sub-array size: " << next_arrsize << ", The number of output ciphertexts: " << resultnumct  << endl;
    
    // cout << "var1" << var1 << "var2 " << var2 << endl;
    TimeVar t;
    // Plaintext pt;

    // check_slot_validity(ciphertext, arrsize, bound, keys, numk, false);

    TIC(t);
    
    // Rotate
    #pragma omp parallel for if(multithread)
    for(int32_t i=0;i< numct;i++){
        int32_t idx = i % numpacked;
        sorted[i] = cc->EvalMult(ciphertext[i], masks[2]); //colmask
        if(idx != 0){
            idx = idx % arrsize + numk * arrsize * (idx / arrsize);
            // cout << i << ": " << idx << endl;
            sorted[i]=cc->EvalRotate(sorted[i], -idx); 
        }
    }

   
    // Aggregate
    for(int32_t i=0;i< log2(numpacked);i++){
        int32_t interval = 1 << i;
        #pragma omp parallel for if(multithread)
        for(int32_t j=0;j < numct / (2*interval) ; j++){
            sorted[j*interval*2]=cc->EvalAdd(sorted[j*interval*2], sorted[j*interval*2 + interval]);
        }
    }

    // Rearrange
    #pragma omp parallel for if(multithread && (numpacking > 1))
    for(int32_t i=0;i< numpacking; i++){
        cc->ModReduceInPlace(sorted[i*numpacked]); 
        if(var1 > 1)sorted[i*numpacked] = RotAndSum(sorted[i*numpacked], ((next_arrsize*next_arrsize/var1 - numct)*var1), (next_arrsize*next_arrsize/var1 - numct)); //gather to construct full rows
        if(var2 > 1)sorted[i*numpacked] = RotAndSum(sorted[i*numpacked], ((next_arrsize*next_arrsize/(var1*var2) - numk*arrsize* max(numpacked / arrsize, 1)))*var2, (next_arrsize*next_arrsize/(var1*var2) - numk*arrsize* max(numpacked / arrsize, 1))); //gather to construct full columns
    }

    // Boot
    // #pragma omp parallel for if(multithread && (numpacking > 2))
    for(int32_t i=0;i< numpacking; i++){
        cc->ModReduceInPlace(sorted[i*numpacked]); 
        if(boot1){
            if(verbose && i==1)cout << "Boot: " <<  sorted[i*numpacked]->GetLevel() << endl;
            gathered[i] = cc->EvalBootstrap(sorted[i*numpacked],2,10);
        }else{
            gathered[i] = fakeboot(sorted[i*numpacked], cc, keys, verbose);
        }
    }

    vector<Plaintext> masks2 = maskPrecompute(next_arrsize, batchSize, 2*bound, cc);
    vector<Ciphertext<DCRTPoly>> result(resultnumct);


    #pragma omp parallel for if(multithread)
    for(int32_t i=0;i< resultnumct;i++){
        int32_t idx0 = i / resultnumpacked;
        int32_t idx1 = i % resultnumpacked;
        result[i] = gathered[idx0]->Clone();
        if(idx1 != 0)result[i] = cc->EvalRotate(result[i], idx1*next_arrsize); 
        result[i] = cc->EvalMult(result[i], masks2[1]);
        cc->ModReduceInPlace(result[i]); 
    }


    cout << "Time " << TOC(t) <<", level: " << result[0]->GetLevel() <<  endl;

    // check_slot_validity(result, next_arrsize, bound, keys, numk, true);
    // cc->Decrypt(keys.secretKey, result[0], &pt);
    // std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
	// printpt(pt, 256, 64, 1, false, true);


    return result;
}




vector<Ciphertext<DCRTPoly>> bootPacked(const vector<Ciphertext<DCRTPoly>> ciphertext, const int32_t arrsize, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys,  const int32_t numk, const int32_t numct, const bool bootstrapping_int, const bool boot1, const bool multithread, const bool verbose){
    const CryptoContext<DCRTPoly> cc = ciphertext[0]->GetCryptoContext();
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize(); 

    vector<Plaintext> masks = maskPrecompute(arrsize, batchSize, 1, cc, false, numk); // 0: colmask 1: rowmask 2: topkmask

    int32_t numpacking = max(numct/arrsize, 1); //number of packing ciphertexts for bootstrapping
    int32_t numpacked = numct / numpacking; //number of ciphertexts packed in an intermediate ciphertext

    vector<Ciphertext<DCRTPoly>> gathered(numpacking);
    vector<Ciphertext<DCRTPoly>> result(numct);

    cout << "Merge and divide, num of ciphertext: " << numct << ", number of Bootstrapping ciphertexts: " << numpacking << ", current sub-array size: " << arrsize << endl;
    
    TimeVar t;

    // Rotate
    #pragma omp parallel for if(multithread)
    for(int32_t i=0;i< numct;i++){
        int32_t idx = i % arrsize;
        if(bootstrapping_int){
            result[i] = cc->EvalMult(ciphertext[i], masks[3]); //rowmask, embedd to [0,1]
        }else{
            result[i] = cc->EvalMult(ciphertext[i], masks[1]); //rowmask
        }
        if(idx != 0){
            result[i]=cc->EvalRotate(result[i], -idx*arrsize); 
        }
    }

    // Aggregate
    for(int32_t i=0;i< log2(numpacked);i++){
        int32_t interval = 1 << i;
        #pragma omp parallel for if(multithread)
        for(int32_t j=0;j < numct / (2*interval) ; j++){
            // cout << j*interval*2 << ", " << j*interval*2 + interval << endl;
            result[j*interval*2]=cc->EvalAdd(result[j*interval*2], result[j*interval*2 + interval]);
        }
    }

    // Boot
    #pragma omp parallel for if(multithread && (numpacking > 2))
    for(int32_t i=0;i< numpacking; i++){
        cc->ModReduceInPlace(result[i*numpacked]); 
        if(boot1){
            if(verbose && i==1)cout << "Boot: " <<  result[i*numpacked]->GetLevel() << endl;
            gathered[i] = cc->EvalBootstrap(result[i*numpacked],2,10);
        }else{
            gathered[i] = fakeboot(result[i*numpacked], cc, keys, verbose);
        }
    }

    #pragma omp parallel for if(multithread)
    for(int32_t i=0;i< numct;i++){
        int32_t idx0 = i / arrsize;
        int32_t idx1 = i % arrsize;
        result[i] = gathered[idx0]->Clone();
        if(idx1 != 0)result[i] = cc->EvalRotate(result[i], idx1*arrsize);
        if(bootstrapping_int){
            result[i] = cc->EvalMult(ciphertext[i], masks[4]); //rowmask, un-embedd from [0,1]
        }else{
            result[i] = cc->EvalMult(ciphertext[i], masks[1]); //rowmask
        }
        cc->ModReduceInPlace(result[i]); 
    }
    
    cout << "Time " << TOC(t) <<", level: " << result[0]->GetLevel() <<  endl;

    return result;
}

// vector<Ciphertext<DCRTPoly>> sort_multi(const Ciphertext<DCRTPoly> ciphertext, const int32_t size, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys, const bool boot1, const bool boot2){
// 	const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
//     int32_t batchSize = cc->GetEncodingParams()->GetBatchSize(); 

//     usint sizesquare = size*size;
//     usint numct = sizesquare / batchSize;
//     usint slicelength = size / numct;

//     vector<Plaintext> masks = maskPrecompute_full(size, batchSize, 2*bound, cc);

//     Ciphertext<DCRTPoly> tmp, tmp2, copy, result;
//     // vector<Ciphertext<DCRTPoly>> tmps2(numct);

//     Plaintext pt, ptboot;

//     cout << "numct: " << numct << ", slicelength: " << slicelength << endl;

//     copy =  RotAndSum(ciphertext, -batchSize, -size); // copy row, Note: size*slicelength = batchSize

    
//     Ciphertext<DCRTPoly> transposed = RotAndSum(ciphertext, -(size-1)*slicelength, -(size-1));//transpose
//     for(usint i=0; i<numct; i++){
//         // cout << i << endl;
//         tmp = cc->EvalRotate(transposed, i*slicelength);
//         tmp = cc->EvalMult(tmp, masks[0]);
//         cc->ModReduceInPlace(tmp);
//         tmp =  RotAndSum(tmp, -size, -1); //copy col
//         tmp = cc->EvalSub(tmp, copy); // ct_cc-ct_rc
//         tmp = cc->EvalAdd(tmp, tmp);
//         tmp = cc->EvalAdd(tmp, masks[2+i]);
//         tmp = comp(tmp, 2*bound, true); //obtain M_comp

//         // auto tmplevel = tmp->GetLevel();
//         // cout << tmplevel << endl;
//         // if(boot1 && tmplevel > 20)tmp = cc->EvalBootstrap(tmp,2,10);

//         if(i==0){
//             // result = RotAndSum(tmp, batchSize, size); // row sum
//             result = tmp->Clone();
//         }else{
//             // tmp = RotAndSum(tmp, batchSize, size); // row sum
//             // result = cc->EvalAdd(tmp, result);
//             result = cc->EvalAdd(tmp, result);
//         }
//     }


 
//     // cc->Decrypt(keys.secretKey, result, &pt);
//     // std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
//     // pt->SetLength(size);
//     // cout << "step1: " << pt << endl;

 



//     // result = cc->EvalMult(result, masks[1]);
//     // cc->ModReduceInPlace(result);
//     // result = RotAndSum(result, -batchSize, -size); //copy row

//     result = RotAndSum(result, -batchSize, -size); //row sum & copy row


//     vector<usint> rounds = GenIndicatorRounds(size, scaleModSize);

//     vector<Ciphertext<DCRTPoly>> results(numct);
//     for(usint i=0; i<numct; i++){
//         Plaintext numstocheck = GenIndicatorCheckerForSort(size, cc, i);
//         tmp2 = IndicatorSIMD(result, size, rounds, numstocheck);
//         tmp2 = cc->EvalMult(tmp2, copy);
//         cc->ModReduceInPlace(tmp2);
//         results[i] = RotAndSum(tmp2, size, 1); // col sum
//     }


//     // const double div =1 / (double) size;
//     // tmp = cc->EvalSub(tmp, numstocheck);
//     // cc->EvalMultInPlace(tmp, div);
//     // if(boot2==true)tmp = cc->EvalBootstrap(tmp);
//     // usint cleanseiter=rounds[1];
//     // rounds[1]=0;
//     // tmp2 = Indicator(tmp, 1, rounds, 0);
//     // if(boot2==true)tmp2 = cc->EvalBootstrap(tmp2);
// 	// tmp2 = Cleanse(tmp2, cleanseiter);

//     // cc->Decrypt(keys.secretKey, tmp2, &pt);
//     // std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
//     // pt->SetLength(sizesquare);
//     // pt->SetLength(8);
//     // cout << "Indicator: " << pt << endl;


//     // cc->Decrypt(keys.secretKey, result, &pt);
//     // std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
//     // pt->SetLength(size);
//     // cout << "step3: " << pt << endl;


//     return results;
// }

vector<Ciphertext<DCRTPoly>> encryptForTopk(const vector<double> vals, const int32_t size, const int32_t numk, const PublicKey<DCRTPoly> publicKey, CryptoContext<DCRTPoly> cc, const int32_t arrsize_init){
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize(); 
    int32_t numct = 1;
    int32_t arrsize = 2*numk;
    
    if(arrsize_init == 0){
        if(log2(batchSize) < 1+ log2(numk*size)){ // Multiple Ciphertext Case
            numct = 2*numk*size/batchSize;
            arrsize = 2*numk;
        }else{
            arrsize = batchSize/size;
            if(log2(size)*2 <= log2(batchSize)){ // Single Ciphertext Case
                arrsize = size;
            }
        }
    }else{
        arrsize = arrsize_init;
        numct = arrsize_init * size / batchSize;
    }
    


    int32_t arrsizesquare = arrsize*arrsize;
    vector<Ciphertext<DCRTPoly>> c1(numct);
    int32_t numofsubarray = size/arrsize;
    cout << "!!!!!!!! Encryption: number of input ciphertexts: " << numct << ", sub array length: " << arrsize << ", number of sub array: "  << numofsubarray <<  " !!!!!!!!" <<endl;
    // vector<double> msg(batchSize);
    // for(usint j=0; j<batchSize; j++){
    //     msg[j]=0;
    // }



    for(int32_t i=0; i<numct; i++){
        vector<double> msg(batchSize);
        for(int32_t j=0; j<numofsubarray/numct; j++){
            for(int32_t k=0; k<arrsize; k++){
                msg[j*arrsizesquare+k]=vals[i*size/numct + j*arrsize+k];
            }
        }
        Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(msg);
        c1[i] = cc->Encrypt(publicKey, ptxt1);

        
    }
    return c1;
}


Ciphertext<DCRTPoly> blocktopkS(const Ciphertext<DCRTPoly> ciphertext, const int32_t arrsize, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys, const int32_t numk, const bool boot1, const bool boot2, const bool verbose){
	const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize(); 

    int32_t arrsizesquare = arrsize*arrsize;
    // vector<Plaintext> masks = maskPrecomputeForTopk(arrsize, batchSize, numk, cc); // 0: colmask 1: rowmask
    vector<Plaintext> masks = maskPrecompute(arrsize, batchSize, 1, cc, false); // 0: colmask 1: rowmask


    Plaintext targetranks = GenIndicatorCheckerForSort(arrsize, cc, numk);
    Ciphertext<DCRTPoly> c1 = sort(ciphertext, arrsize, bound, scaleModSize, keys,  boot1, boot2, targetranks, 1, verbose);

    if((c1->GetLevel()) > 47){
        if(boot1){
            if(verbose)cout << "Boot: " <<  c1->GetLevel() << endl;
            c1 = cc->EvalBootstrap(c1,2,10);
        }else{
            if(boot2)c1 = fakeboot(c1, cc, keys);
        }
    }


    c1 = cc->EvalMult(c1, masks[0]);
    cc->ModReduceInPlace(c1);
    c1 = RotAndSum(c1, (arrsize-1)*numk, (arrsize-1));//col transpose
    c1 = cc->EvalMult(c1, masks[1]);
    cc->ModReduceInPlace(c1);

    // int32_t newsize = arrsizesquare/numk;
    // if(2*log2(newsize) < log2(batchSize)){
    //     int32_t tmpsize = newsize/numk;
    //     c1 = RotAndSum(c1, (arrsizesquare-numk)*tmpsize, (arrsizesquare-numk));//gather sub-arrays
    //     vector<Plaintext> masks2 = maskPrecompute(newsize, batchSize, 1, cc, false); // 0: colmask 1: rowmask
    //     c1 = cc->EvalMult(c1, masks2[1]);
    //     cc->ModReduceInPlace(c1);
    // }else{
    //     int32_t tmpsize = batchSize/arrsizesquare;
    //     c1 = RotAndSum(c1, (arrsizesquare-numk)*tmpsize, (arrsizesquare-numk));//gather sub-arrays
    //     vector<Plaintext> masks2 = maskPrecompute(tmpsize*numk, batchSize, 1, cc, false); // 0: colmask 1: rowmask
    //     c1 = cc->EvalMult(c1, masks2[1]);
    //     cc->ModReduceInPlace(c1);
    // }

    int32_t next_arrsize = arrsizesquare/numk;
    if(next_arrsize * next_arrsize > batchSize)next_arrsize = batchSize / next_arrsize;
    c1 = RotAndSum(c1, (arrsizesquare-numk)*next_arrsize/numk, (arrsizesquare-numk));//gather sub-arrays
    vector<Plaintext> masks2 = maskPrecompute(next_arrsize, batchSize, 1, cc, false); // 0: colmask 1: rowmask
    c1 = cc->EvalMult(c1, masks2[1]);
    cc->ModReduceInPlace(c1);


    return c1;
}


vector<Ciphertext<DCRTPoly>> blocktopkM(const vector<Ciphertext<DCRTPoly>> ciphertexts, const int32_t size, const int32_t arrsize, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys, const int32_t numk, const bool boot1, const bool multithread, const int32_t mergecrit){
	const CryptoContext<DCRTPoly> cc = ciphertexts[0]->GetCryptoContext();
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize();

    int32_t numct = size * arrsize / batchSize;
    int32_t nextsize = size * numk / arrsize;

    //Merge Keys                    
    
    vector<Ciphertext<DCRTPoly>> sorted(numct);

    if(numk * nextsize * mergecrit < batchSize){  //merge to one
        int32_t next_arrsize = min(nextsize, batchSize / nextsize);
        int32_t arrsizesquare = arrsize * arrsize ;
        int32_t next_arrsizesquare = next_arrsize * next_arrsize ;

        vector<Plaintext> masks = maskPrecompute(arrsize, batchSize, 1, cc, false, numk); // 0: colmask 1: rowmask 2: topkmask
        Plaintext targetranks = GenIndicatorCheckerForSort(arrsize, cc, numk);

        //Variables for merge
        int32_t var0 = max(numk/max(next_arrsizesquare / arrsize,1),1);
        int32_t var1 = next_arrsize * var0 / numk; 
        int32_t var2 = var0 * max(next_arrsizesquare / arrsize,1) * arrsize;
        int32_t var3 = max(var0, next_arrsize/ arrsize);

        vector<Plaintext> masks_next = maskPrecompute(next_arrsize, batchSize, 1, cc, false, numk); // 0: colmask 1: rowmask
        
        cout << "Merge Phase, num of ciphertext: " << numct << ", current sub-array size: " << arrsize << ", next sub-array size: " << next_arrsize << endl;

        sorted = sort_SIMD(ciphertexts, arrsize, bound, scaleModSize, keys, numct, boot1, targetranks, 1, multithread, true);

        // cout << "sorted" << var0 << ", " << var1 <<", " << var2 <<", " << var3 << endl;

        #pragma omp parallel for if(multithread)
        for(int32_t i=0;i<numct;i++){
            // bool verbose = false;
            // if(i==0)verbose = true;
            // sorted[i]=sort(ciphertexts[i], arrsize, bound, scaleModSize, keys,  boot1, boot2, targetranks, 1, verbose);
            sorted[i]=cc->EvalMult(sorted[i], masks[2]);

            // if(verbose)cout << next_arrsize << ", " << arrsize << endl;
            if(i!=0){
                int32_t idx0 = i % var1;
                int32_t idx1 = i / var1;
                sorted[i]=cc->EvalRotate(sorted[i], -(idx0 + var2 * idx1));
            }
        }
    


        for(int32_t i=0;i< log2(numct);i++){
            int32_t interval = 1 << i;
            #pragma omp parallel for if(multithread)
            for(int32_t j=0;j < numct / (2*interval) ; j++){
                sorted[j*interval*2]=cc->EvalAdd(sorted[j*interval*2], sorted[j*interval*2 + interval]);
            }
        }

        vector<Ciphertext<DCRTPoly>> result(1);
 
        //Merge ciphertexts to one
        // #pragma omp parallel for if(multithread && (numct/multiplier > 1))
        result[0] = sorted[0]->Clone();
        cc->ModReduceInPlace(result[0]);

        int32_t rotvar0 = next_arrsize / (numk * numct);
        int32_t rotvar1 = numk / var3;
        if(rotvar0 > 1)result[0]=RotAndSum(result[0], (arrsizesquare - numct) * rotvar0, (arrsizesquare - numct) ); // gathering
        if(rotvar1 > 1)result[0]=RotAndSum(result[0], (max(arrsize, next_arrsize) - var1) * rotvar1 , max(arrsize, next_arrsize) - var1); // col transpose
        result[0]=cc->EvalMult(result[0], masks_next[1]);
        cc->ModReduceInPlace(result[0]);
        
        return result;
    }else{
        cout << "num of ciphertext: " << numct << ", current sub-array size: " << arrsize << endl;
        Plaintext targetranks = GenIndicatorCheckerForSort(arrsize, cc, numk);

        sorted = sort_SIMD(ciphertexts, arrsize, bound, scaleModSize, keys, numct, boot1, targetranks, 1, multithread, true);
        vector<Ciphertext<DCRTPoly>> result = rearrange_SIMD(sorted, arrsize, bound, scaleModSize, keys, numk, numct, boot1, multithread, true);

        return result;
    }


    // return 0;
}


vector<int32_t> schedule_topk(const int32_t size, const int32_t numk, int32_t init_arrsize_bound, int32_t init_numct_bound, bool binarymerge,  const int32_t boot_threshold, const int32_t mergecrit, const bool mute){
    int32_t batchSize = 1 << 16;
    if(mute == false)std::cout << "\nTest on Size: " << size << ", k: " << numk << ", bootthreads: " << boot_threshold << std::endl;

    //Round Estimation
    int32_t arrsize=2*numk;
    int32_t currentsize=size;
    int32_t init_arrsize = 2*numk;
    int32_t numct = 1; 
    int32_t roundM = 0;
    int32_t roundS = 0;

    if(binarymerge){
        init_arrsize = 2 * numk;
        arrsize = 2 * numk;
        while(currentsize >= batchSize / numk){
            if(mute == false)cout << " Binary merge: Size " <<  currentsize << ", length of sub-array: " << arrsize << endl;
            currentsize /=2;
        }
    }else{

        numct = init_arrsize_bound * currentsize / batchSize;
        init_arrsize = init_arrsize_bound;
        if(init_numct_bound < numct){
            numct = init_numct_bound;
            init_arrsize = numct * batchSize / currentsize;
        }
        if(numct < 1){
            numct = 1;
            init_arrsize = batchSize / currentsize;
        }

        arrsize = init_arrsize;
        if(mute == false)cout << "initial number of ciphertexts: " << numct << ", number of candidates: "<< currentsize << ", length of sub-array: " << arrsize << endl;

        while(numct > 1){
            int32_t multiplier = 1;
            if(mute == false)cout << "BlockTopkM: number of ciphertexts: " << numct << ", number of candidates: "<< currentsize << ", length of sub-array: " << arrsize << endl;
            currentsize = currentsize * numk / arrsize;
            
            //condition check
            if(numk * currentsize * mergecrit < batchSize){ //merge to one
                multiplier = numct;
            }else{
                int32_t target_output_bound = std::min((int32_t)(256*currentsize / batchSize), boot_threshold); // number of output: (condition 2, empty sub-matrix occurence) (condition 3, initial merge for threads_actual)
                if( target_output_bound < numct ){
                    multiplier = numct / target_output_bound; 
                }
            }
            numct /=multiplier;
            arrsize = numct * batchSize / currentsize;
            if(arrsize <= numk)return {256,256};
            roundM +=1;
        }
    }
    
    while(currentsize*currentsize > batchSize){
        if(mute == false)cout << "BlockTopkS: single ciphertext: number of candidates: "<< currentsize << ", length of sub-array: " << arrsize << endl;
        currentsize = currentsize*numk/arrsize;
        arrsize = batchSize / currentsize;
        if(arrsize <= numk)return {256,256};
        roundS +=1;
    }

    arrsize = currentsize;
    if(mute == false){
        cout <<  "Final Rankselect with number of candidates:" <<  currentsize << endl;
        return {init_arrsize, arrsize};
    }else{
        return {roundM, roundS};
    }

}

// Round Estimator by summation
// int32_t roundestimator(const int32_t size, const int32_t numk, int32_t init_arrsize_bound, int32_t init_numct_bound, int32_t rounds_budget){

//     vector<int32_t> best_rounds = {100, 100};
//     for(int32_t i =0 ; i < log2(init_numct_bound); i++){
//         int32_t init_numct = init_numct_bound >> i;
//         vector<int32_t> rounds = schedule_topk(size, numk, init_arrsize_bound, init_numct, false,  128, true);
//         if(best_rounds[0] > rounds[0]){
//             best_rounds[0] = rounds[0];
//             best_rounds[1] = rounds[1];
//         }else{
//             if((best_rounds[0] < rounds[0]) || (best_rounds[1] + rounds_budget  < rounds[1])){
//                 cout << "Best number of rounds :" << best_rounds[0] << ", " << best_rounds[1] << " on inital number of ciphertext: " <<  (init_numct << 1) << endl;
//                 return init_numct << 1;
//             }
//         }
//     }
    
//     return init_numct_bound;
// }


int32_t roundestimator(const int32_t size, const int32_t numk, int32_t init_arrsize_bound, int32_t init_numct_bound, const int32_t mergecrit){

    int32_t result_init_numct = 256;

    vector<int32_t> best_rounds = {255, 255};
    for(int32_t i =0 ; i < log2(init_numct_bound); i++){
        int32_t init_numct = init_numct_bound >> i;
        vector<int32_t> rounds = schedule_topk(size, numk, init_arrsize_bound, init_numct, false,  128, mergecrit, true);
        
        if((best_rounds[0] == rounds[0]) && (best_rounds[1] == rounds[1])){
            result_init_numct = init_numct;
        }

        if(best_rounds[0] + best_rounds[1] > rounds[0] + rounds[1]){
            best_rounds[0] = rounds[0];
            best_rounds[1] = rounds[1];
            result_init_numct = init_numct;

        }else{
            if(best_rounds[0] + best_rounds[1] < rounds[0] + rounds[1]){
                cout << "Best number of rounds :" << best_rounds[0] << ", " << best_rounds[1] << " on inital number of ciphertext: " <<  result_init_numct << endl;
                return result_init_numct;
            }
        }
    }
    
    return init_numct_bound;
}


Ciphertext<DCRTPoly> topk_multithread_opt(const vector<Ciphertext<DCRTPoly>> ciphertexts, const int32_t size, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys, const int32_t numk, const bool boot1, const bool boot2, const bool multithread, const int32_t mergecrit){
	const CryptoContext<DCRTPoly> cc = ciphertexts[0]->GetCryptoContext();
    int32_t batchSize = cc->GetEncodingParams()->GetBatchSize(); 

    
    int32_t arrsize = 2*numk;
    int32_t currentsize = size;
    int32_t round = 1;

    int32_t numofct = ciphertexts.size();

    vector<Ciphertext<DCRTPoly>> sorted(numofct);
    bool init=true;

    // cout << "numcofct: " << numofct << endl;


    while(numofct > 1){

        arrsize = batchSize * numofct / currentsize;
        if(arrsize * arrsize > batchSize) arrsize = batchSize / arrsize;
        cout << "!!!!!!!! ROUND " << round << ": number of ciphertexts:" << numofct << ", arrsize: " << arrsize << ", number of elements: " << currentsize << endl;
        round +=1;
        if(init){
            sorted = blocktopkM(ciphertexts, currentsize, arrsize, bound, scaleModSize, keys, numk, boot1, multithread, mergecrit);
            init = false;
        }else{
            sorted = blocktopkM(sorted, currentsize, arrsize, bound, scaleModSize, keys, numk, boot1, multithread, mergecrit);
        }
        currentsize = currentsize * numk / arrsize;
        numofct = sorted.size(); 

    }

    while(currentsize*currentsize > batchSize){
        arrsize = batchSize/currentsize;
        cout << "!!!!!!!! ROUND " << round <<  ": arrsize:" << arrsize << ", number of elements: " << currentsize << endl;
        round +=1;
        if(init){
            sorted[0] = blocktopkS(ciphertexts[0], arrsize, bound, scaleModSize, keys, numk, boot1, boot2);
            init = false;
        }else{
            sorted[0] = blocktopkS(sorted[0], arrsize, bound, scaleModSize, keys, numk, boot1, boot2);
        }

        // Plaintext pt;
        // cc->Decrypt(keys.secretKey, sorted[0], &pt);
        // std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
        // pt->SetLength(256);
        // cout << pt << endl;

        currentsize = currentsize*numk/arrsize;
    }

    cout << "!!!!!!!! ROUND " << round << ": Final arrsize:" << currentsize << ", number of elements: " << currentsize << endl;

    Plaintext targetranks = GenIndicatorCheckerForSort(currentsize, cc, numk);
    sorted[0] = sort(sorted[0], currentsize, bound, scaleModSize, keys, boot1, boot2, targetranks, 1);
 
    return sorted[0];
}



Ciphertext<DCRTPoly> nexus(const Ciphertext<DCRTPoly> ciphertext, const int32_t size, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys){
	const CryptoContext<DCRTPoly> cc = ciphertext->GetCryptoContext();
    // int32_t batchSize = cc->GetEncodingParams()->GetBatchSize();
    // uint32_t currentsize = size;

    cout << "Start on array size: " << size << endl;
    Ciphertext<DCRTPoly> tmp = cc->EvalRotate(ciphertext, -size);
    Ciphertext<DCRTPoly> result = cc->EvalAdd(ciphertext, tmp);
    Ciphertext<DCRTPoly> cmp;

    for(usint i=0;i< log2(size); i++){
        cout << "Step " << i << " per " << log2(size) << endl;

        tmp = cc->EvalRotate(result, 1 << i);
        result = cc->EvalSub(result, tmp);
        cmp = comp(result, bound, true, 45);

        result = cc->EvalMult(result, cmp);
        cc->ModReduceInPlace(result);
        result = cc->EvalAdd(result, tmp);

        

        // Plaintext pt;
        // cc->Decrypt(keys.secretKey, result, &pt);
        // std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
        // pt->SetLength(256);
        // cout << pt << endl;
    
        if((result->GetLevel()) > 30){
            cout << "Boot: " <<  result->GetLevel() << endl;
            // cc->Decrypt(keys.secretKey, copy, &pt);
            // std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
            // pt->SetLength(256);
            // cout << pt << endl;
            result = cc->EvalBootstrap(result,2,10);
            // cc->Decrypt(keys.secretKey, copy, &pt);
            // std::cout << "Estimated level: " << pt->GetLevel() << std::endl;
            // pt->SetLength(256);
            // cout << pt << endl;

            // result = fakeboot(result, cc, keys);

        }
    }

    cout << "Final " << endl;

    result = cc->EvalSub(ciphertext, result);
    result = comp(result, bound, true, 45);
    result = cc->EvalAdd(result,result);

    return result;
}


}