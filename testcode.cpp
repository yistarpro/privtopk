#include "openfhe.h"
#include "testcode.h"
#include "algorithms.h"
#include "utils.h"

#include <omp.h>

using namespace lbcrypto;
using namespace std;

namespace ckkseif {

    string statTime(const vector<double> times, const usint iteration){
        double avg=0.0;
        double std=0.0;

        if(iteration!=1){
            for(long j=0;j<iteration;j++)avg+=times[j];
            avg/=iteration;
            for(long j=0;j<iteration;j++)std+=(times[j]-avg)*(times[j]-avg);
            std/=iteration;
            std=sqrt(std);
            cout << "Average time = " << avg << ", Std =" << std << endl;
            return "Average time = "+to_string(avg) + ", Std =" +to_string(std);
        }else{
            cout << "Average time = " << times[0] << endl;
            return "Average time = "+to_string(times[0]);
        }
    }



    void SortTest(const uint32_t scaleModSize, const usint size, const int32_t arraybound,  const usint levelBudgetElmt, const usint logbatchsize, usint iteration) {
        TimeVar t;
        vector<double> timeEval(iteration);

        // const size_t numThreads = omp_get_num_procs();;
        // std::cout << "Number of threads     : " << numThreads << std::endl;
        // omp_set_num_threads(numThreads);
        // omp_set_max_active_levels(10);

        // if (numCiphertext <= numThreads / 16)omp_set_max_active_levels(10);

        // uint32_t multDepth = 2320/scaleModSize;
        uint32_t multDepth = 55;

        cout << "multDepth = " << multDepth << endl;
        // uint32_t multDepth = 58;
        uint32_t batchSize = 1 << logbatchsize;


        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetMultiplicativeDepth(multDepth);
        parameters.SetScalingModSize(scaleModSize);
        parameters.SetRingDim(batchSize << 1);
        parameters.SetBatchSize(batchSize);
        
        parameters.SetSecretKeyDist(SPARSE_TERNARY);

        // parameters.SetSecurityLevel(HEStd_NotSet);
        // /// Bootstrap block 1 ////
        // parameters.SetFirstModSize(scaleModSize+1);
        // parameters.SetScalingTechnique(FIXEDAUTO);
        

        //cout << "CKKS standard deviation " << parameters.GetStandardDeviation() << endl;
        //cout << "CKKS security level " <<  parameters.GetSecurityLevel() << endl;

        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

        // Enable the features that you wish to use
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);

        if(levelBudgetElmt > 0) cc->Enable(FHE);



        //std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() << std::endl << std::endl;

        paramcheck(cc);

        // B. Step 2: Key Generation
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
        AddRotKeyForSort(keys.secretKey, cc, size);
	    //bootSet2(cc, keys.secretKey, batchSize);
        if(levelBudgetElmt > 0) {
            std::vector<uint32_t> levelBudget = {levelBudgetElmt, levelBudgetElmt};
            cc->EvalBootstrapSetup(levelBudget);
            cout << "Boot Setup Done" << endl;
            cc->EvalBootstrapKeyGen(keys.secretKey, batchSize);
        }

        // Step 3: Encoding and encryption of inputs
        std::cout << "!!!!!!!!!!!!!!! Sort Test !!!!!!!!!!!!!!!" << std::endl;
        std::cout << "\nTest on bound: " << arraybound << ", Size: " << size <<  std::endl;

        // Inputs
        std::vector<double> x1(batchSize);
        x1 = fixedDiscreteArrayHalf(batchSize, arraybound);
        for(usint i=size;i<batchSize;i++)x1[i]=0.0;
        Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x1);
        ptxt1->SetLength(8);
        std::cout << "\n Input x1: " << ptxt1 << std::endl;
        ptxt1 = cc->MakeCKKSPackedPlaintext(x1);

        // Encrypt the encoded vectors
        auto c1 = cc->Encrypt(keys.publicKey, ptxt1);
        Plaintext result;

        // Plaintext numstocheck = GenIndicatorCheckerForSort(size, cc, 0);
        // numstocheck->SetLength(256);
        // cout << numstocheck << endl;
        Plaintext targetranks = GenIndicatorCheckerForSort(size, cc, 0);

        for(usint i=0;i<iteration;i++){
            double diff = 0.0;
            TIC(t);
            auto c2 = sort(c1, size, arraybound, scaleModSize, keys, false, false, targetranks);
            timeEval[i] = TOC(t);
            cc->Decrypt(keys.secretKey, c2, &result);
            vector<double> vals1 = result->GetRealPackedValue();
            vector<double> vals2(size);
            bool sorted = true;
            for(usint j=0;j<size;j++){
                vals2[j]=vals1[j*size];
            }
            for(usint j=1;j<size;j++){
                if(vals2[j] - vals2[j-1] > diff)diff=vals2[j] - vals2[j-1];
            }

            if(diff > 1.0/((double)arraybound))sorted = false;
            //precision(result, batchSize);
            result->SetLength(size);
            cout << vals2 << " :: Sorted: " << sorted << ", " << diff <<  " :: " << "\nEstimated level: " << result->GetLevel() << std::endl;
            cout << "Time: " << timeEval[i] << endl;
        }
        statTime(timeEval, iteration);

    } 

    void TopkTest_MT_opt(const uint32_t scaleModSize, const usint size, const int32_t arraybound, const usint numk, const usint levelBudgetElmt, usint iteration, const int32_t mtenvironment, const int32_t initial_subsize_bound, const int32_t init_numct_bound, const bool optimized, const int32_t custom_mergecrit ) {
        TimeVar t;
        vector<double> timeEval(iteration);

        // uint32_t multDepth = 2320/scaleModSize;
        uint32_t multDepth = 2500/scaleModSize;

        bool multithread = true;

        if(mtenvironment!=0){
            if(mtenvironment==1){
                omp_set_num_threads(1);
                std::cout << "Single Thread Version " << std::endl;
                multithread = false;
            }else{
                omp_set_num_threads(mtenvironment);
            }

        }
        // omp_set_max_active_levels(10);
        const size_t numThreads_max = omp_get_num_procs();
        const size_t numThreads = omp_get_max_threads();

        std::cout << "Number of threads     : " << numThreads <<  " out of " << numThreads_max << std::endl;


        // uint32_t multDepth = 55;

        cout << "multDepth = " << multDepth << endl;
        uint32_t batchSize = 1 << 16;
        bool boot1=false;


        std::cout << "\nTest on bound: " << arraybound << ", Size: " << size << ", k: " << numk << ", threads: " << numThreads <<", initial sub-array size boud: " << initial_subsize_bound << ", initial number of ciphertexts bound: " << init_numct_bound << std::endl;

        //Round Estimation
        int32_t init_numct_bound_opt = init_numct_bound;
        int32_t mergecrit = custom_mergecrit;
        if(optimized){
            if(numk % 2 == 0){
                mergecrit = 16 / sqrt(2*numk);
            }else{
                mergecrit = 16 / sqrt(4*numk);
            }
            std::cout << "Merge Cirterion updated to : "<<  batchSize / (mergecrit * numk) << endl; 
        
            init_numct_bound_opt = roundestimator(size, numk, initial_subsize_bound, init_numct_bound, mergecrit);
        }
       	vector<int32_t> subsizes = schedule_topk(size, numk, initial_subsize_bound, init_numct_bound_opt, false, 128, mergecrit, false); 


        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecretKeyDist(SPARSE_TERNARY);
        parameters.SetMultiplicativeDepth(multDepth);
        parameters.SetScalingModSize(scaleModSize);
        parameters.SetRingDim(batchSize << 1);
        parameters.SetBatchSize(batchSize);
        
        parameters.SetFirstModSize(scaleModSize+1);
    
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

        // Enable the features that you wish to use
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        if(levelBudgetElmt > 0) cc->Enable(FHE);


        //std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() << std::endl << std::endl;

        paramcheck(cc);
        std::vector<uint32_t> levelBudget = {levelBudgetElmt, levelBudgetElmt};
        if(levelBudgetElmt > 0) {
            usint depth = FHECKKSRNS::GetBootstrapDepth(levelBudget, parameters.GetSecretKeyDist());
            cout << "scaleModSize: " << scaleModSize << endl;
            cout << "bootdepth: " << depth << ", levelBudget: " << levelBudgetElmt << endl;
            cout << "budgetdepth: " << multDepth-depth << endl;
            cout << "SlotDim: " << 16 << endl;
            boot1=true;
        }
        // B. Step 2: Key Generation
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
	    //bootSet2(cc, keys.secretKey, batchSize);
        if(levelBudgetElmt > 0) {
            cc->EvalBootstrapSetup(levelBudget);
            cc->EvalBootstrapKeyGen(keys.secretKey, batchSize);
            cout << "Boot Setup Done" << endl;
        }
        AddRotKeyForTopk_multithread_opt(keys.secretKey, cc, size, numk, subsizes[0], mergecrit);


        // Step 3: Encoding and encryption of inputs
        std::cout << "!!!!!!!!!!!!!!! Topk Multi-threading Test !!!!!!!!!!!!!!!" << std::endl;

        // Inputs
        std::vector<double> x1(batchSize);
        // x1 = fixedUBDiscreteArrayHalf(size, arraybound, batchSize, true);
        x1 = fixedUBDiscreteArrayHalf(size, arraybound, batchSize);

        // x1 = equalvalueArray(size, arraybound, batchSize);

 

        vector<Ciphertext<DCRTPoly>> c1 = encryptForTopk(x1, size, numk, keys.publicKey, cc, subsizes[0]);
        Plaintext result;

        // cc->Decrypt(keys.secretKey, c1[1], &result);
        // result->SetLength(256);
        // cout << result << endl;


        for(usint i=0;i<iteration;i++){
            TIC(t);
            auto c2 = topk_multithread_opt(c1, size, arraybound, scaleModSize, keys, numk, boot1, false, multithread, mergecrit);
            timeEval[i] = TOC(t);
            cout << "Time: " << timeEval[i] << endl;

            cc->Decrypt(keys.secretKey, c2, &result);
            //precision(result, batchSize);
            printpt(result, numk, 0, subsizes[1], false, true);
            printpt(result, numk, 0, subsizes[1], false, false, 256);
            vector<double> vals1 = result->GetRealPackedValue();
            double err = 0;
            int32_t counter = 0;
            int32_t c = arraybound-1;
            for(usint j = 0; j < numk; j++){
                double val = ((double)c)/((double)arraybound*2.0) - vals1[j*subsizes[1]];
                if(val < 0)val = -val;
                if(val > err)err=val;                

                counter+=1;
                if(counter+c==arraybound){
                    counter=0;
                    if(c>0)c-=1;
                }
            }

            cout << "\nEstimated level: " << result->GetLevel() << ", precision: " << log2(err) << std::endl;
            // cout << "Time: " << timeEval[i] << endl;
        }

        statTime(timeEval, iteration);

    } 

	void nexusTest(const uint32_t scaleModSize, const usint size, const int32_t arraybound, const usint levelBudgetElmt, usint iteration){
        TimeVar t;
        vector<double> timeEval(iteration);

        // uint32_t multDepth = 2320/scaleModSize;
        uint32_t multDepth = 2500/scaleModSize;

        // uint32_t multDepth = 55;

        cout << "multDepth = " << multDepth << endl;
        uint32_t batchSize = 1 << 16;


        std::cout << "\nTest NEXUS on bound: " << arraybound << ", Size: " << size << std::endl;


        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecretKeyDist(SPARSE_TERNARY);

        parameters.SetMultiplicativeDepth(multDepth);
        parameters.SetScalingModSize(scaleModSize);
        parameters.SetRingDim(batchSize << 1);
        parameters.SetBatchSize(batchSize);
        

        // parameters.SetSecurityLevel(HEStd_NotSet);
        // /// Bootstrap block 1 ////
        parameters.SetFirstModSize(scaleModSize+1);
        // parameters.SetScalingTechnique(FIXEDAUTO);
        

        //cout << "CKKS standard deviation " << parameters.GetStandardDeviation() << endl;
        //cout << "CKKS security level " <<  parameters.GetSecurityLevel() << endl;

        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

        // Enable the features that you wish to use
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        if(levelBudgetElmt > 0) cc->Enable(FHE);


        //std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() << std::endl << std::endl;

        paramcheck(cc);
        std::vector<uint32_t> levelBudget = {levelBudgetElmt, levelBudgetElmt};
        if(levelBudgetElmt > 0) {
            usint depth = FHECKKSRNS::GetBootstrapDepth(levelBudget, parameters.GetSecretKeyDist());
            cout << "scaleModSize: " << scaleModSize << endl;
            cout << "bootdepth: " << depth << ", levelBudget: " << levelBudgetElmt << endl;
            cout << "budgetdepth: " << multDepth-depth << endl;
            cout << "SlotDim: " << 16 << endl;
        }
        // B. Step 2: Key Generation
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
	    //bootSet2(cc, keys.secretKey, batchSize);
        if(levelBudgetElmt > 0) {
            cc->EvalBootstrapSetup(levelBudget);
            cc->EvalBootstrapKeyGen(keys.secretKey, batchSize);
            cout << "Boot Setup Done" << endl;
        }
        AddRotKeyForPo2(keys.secretKey, cc, batchSize);


        // Step 3: Encoding and encryption of inputs
        std::cout << "!!!!!!!!!!!!!!! NEXUS Test !!!!!!!!!!!!!!!" << std::endl;

        // Inputs
        std::vector<double> x1(batchSize);
        x1 = fixedUBDiscreteArray(size, arraybound, batchSize);
        // x1 = fixedDiscreteArrayHalf(batchSize, arraybound);
        Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x1);
        // ptxt1->SetLength(256);
        // std::cout << "\n Input x1: " << ptxt1 << std::endl;

        // Encrypt the encoded vectors
        
        // cc->Decrypt(keys.secretKey, c1[0], &result);
        // result->SetLength(512);
        // cout << result << endl;

       
 

        auto c1 = cc->Encrypt(keys.publicKey, ptxt1);
        Plaintext result;

        // cc->Decrypt(keys.secretKey, c1[0], &result);
        // result->SetLength(256);
        // cout << result << endl;


        for(usint i=0;i<iteration;i++){
            TIC(t);
            auto c2 = nexus(c1, size, arraybound, scaleModSize, keys);
            timeEval[i] = TOC(t);
            cout << "Time: " << timeEval[i] << endl;

            cc->Decrypt(keys.secretKey, c2, &result);
            //precision(result, batchSize);
            // printpt(result, numk, 0, finalsub, false, false);
            vector<double> vals1 = result->GetRealPackedValue();
            binaryprecision(result, size);

            for(usint j = 0; j < size; j++){
                if(vals1[j] > 0.6){
                    cout << "Location: " << j << ", mask: " << vals1[j] << endl;
                    break;
                }
            }

            // cout << "\nEstimated level: " << result->GetLevel() << ", precision: " << log2(err) << std::endl;
            // cout << "Time: " << timeEval[i] << endl;
        }

        statTime(timeEval, iteration);

    } 

}
