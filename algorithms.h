#ifndef EIF_ALGORITHMS_H
#define EIF_ALGORITHMS_H

#include "openfhe.h"
#include <iostream>
#include <map>


using namespace lbcrypto;
using namespace std;

namespace ckkseif {

	//----------------------------------------------------------------------------------
	//   Add Rotation Keys for Each Algorithm
	//----------------------------------------------------------------------------------

    std::vector<int32_t> GenIdxForRotsum(const int32_t from, const int32_t to, const bool bidirectional = false);

	std::vector<int32_t> GenIdxForMultiples(const int32_t base, const int32_t number, const bool bidirectional = false);

    void AddRotKeyForPo2(const PrivateKey<DCRTPoly> privateKey, CryptoContext<DCRTPoly> cc, const int32_t batchSize);

    void AddRotKeyForSort(const PrivateKey<DCRTPoly> privateKey, CryptoContext<DCRTPoly> cc, const int32_t arrsize);

	void AddRotKeyForTopk_multithread_opt(const PrivateKey<DCRTPoly> privateKey, CryptoContext<DCRTPoly> cc, const int32_t size, const int32_t numk, const int32_t arrsize_init, const int32_t mergecrit = 1);

    void AddRotKeyForBoot(const PrivateKey<DCRTPoly> privateKey, CryptoContext<DCRTPoly> cc, const int32_t size);

	//----------------------------------------------------------------------------------
	//   ADVANCED HOMOMORPHIC OPERATIONS
	//----------------------------------------------------------------------------------

	Ciphertext<DCRTPoly> BootAuto(const Ciphertext<DCRTPoly> ciphertext);

	Ciphertext<DCRTPoly> BootWithPrec(const Ciphertext<DCRTPoly> ciphertext, const usint verbose, const KeyPair<DCRTPoly> keys);

	Ciphertext<DCRTPoly> Product(const vector<Ciphertext<DCRTPoly>> ciphertext);

    Ciphertext<DCRTPoly> EvalLog(const Ciphertext<DCRTPoly> ciphertext, const double bound, const double base, const usint degree);

    Ciphertext<DCRTPoly> EvalLogLike(const Ciphertext<DCRTPoly> ciphertext, const double bound);

    Ciphertext<DCRTPoly> EvalInverse(const Ciphertext<DCRTPoly> ciphertext, const double bound, const usint degree);

    Ciphertext<DCRTPoly> RotAndSum(const Ciphertext<DCRTPoly> ciphertext, const int32_t from, const int32_t to, const bool verbose = false);
    //16 to 1 : sum up 16 slots
	//Col Sum: RotAndSum(ciphertext, size, 1);
	//Row Sum: RotAndSum(ciphertext, sizesquare, size);
	//Transpose Row to Col: RotAndSum(ciphertext, -(sizesquare-size), -(size-1));
	//Transpose Col to Row: RotAndSum(ciphertext, sizesquare-size, size-1);
	//Copy Col: RotAndSum(ciphertext, -size, -1);
	//Copy Row: RotAndSum(ciphertext, -sizesquare, -size);
	//fullCopy: RotAndSum(ciphertext, -batchSize, -size);


	//----------------------------------------------------------------------------------
	//   Encrypted Indicator Function
	//----------------------------------------------------------------------------------

	vector<double> GetCoeff(const usint bound);

    Plaintext GenIndicatorChecker(const usint bound, const CryptoContext<DCRTPoly> cc);

    Plaintext GenIndicatorCheckerInterval(const usint from, const usint size, const CryptoContext<DCRTPoly> cc);

	//gen pt filled with: from ~ to-1
    Plaintext GenIndicatorCheckerIntervalRecursive(const usint from, const usint to, const usint size, const CryptoContext<DCRTPoly> cc);

	//gen pt filled with: [slicelength * iter] * size, [slicelength * iter + 1] * size .....
	Plaintext GenIndicatorCheckerForSort(const int32_t arrsize, const CryptoContext<DCRTPoly> cc, const int32_t numk=0);

	//gen pt filled with: from ~ to-1
    vector<Plaintext> GenIndicatorCheckerForSIMDCOUNT(const usint base, const usint size, const usint paral, const CryptoContext<DCRTPoly> cc);

    Plaintext GenIndicatorCheckerPartialArray(const usint from, const usint size, const CryptoContext<DCRTPoly> cc, const vector<usint> list);

	//Set indicator rounds for plaintext modulus 35,40,45,50,59
    vector<usint> GenIndicatorRounds(const usint bound, const usint scaleModSize);

	//Set indicator rounds for plaintext modulus 35,40,45,50,59 - not optimized 
    vector<usint> GenZeroTestRounds(const usint bound, const usint scaleModSize);


	Ciphertext<DCRTPoly> Cleanse(const Ciphertext<DCRTPoly> ciphertext, const usint round=1);

	Ciphertext<DCRTPoly> Indicator(const Ciphertext<DCRTPoly> ciphertext, const usint bound, const vector<usint> rounds, const double numtocheck, const usint levlimit=100);

	//Indcator for sort, checking entry of M_comp equal 0. For domain -1, 0, 1, only 0 maps to 1.
	Ciphertext<DCRTPoly> IndicatorBinary(const Ciphertext<DCRTPoly> ciphertext, const vector<usint> rounds);

	Ciphertext<DCRTPoly> IndicatorSIMD(const Ciphertext<DCRTPoly> ciphertext, const usint bound, const vector<usint> rounds, const Plaintext numtocheck, const usint levlimit=100);

	Ciphertext<DCRTPoly> ZeroTest(const Ciphertext<DCRTPoly> ciphertext, const usint bound, const vector<usint> rounds);


	//----------------------------------------------------------------------------------
	//   Comparison & Another Indicator & Parity
	//----------------------------------------------------------------------------------
	Ciphertext<DCRTPoly> normalize(const Ciphertext<DCRTPoly> ciphertext, const double bound);

	//customizable comparison
	Ciphertext<DCRTPoly> comparison(const Ciphertext<DCRTPoly> ciphertext, const usint degf, const usint degg, const double bound, const usint ver);

	//optimized comparison
	Ciphertext<DCRTPoly> comp(const Ciphertext<DCRTPoly> ciphertext, const double bound=1.0, const bool lastmod=false, const usint levlimit=100);
	
	Ciphertext<DCRTPoly> compandUp(const Ciphertext<DCRTPoly> ciphertext, const double bound, const bool boot, const usint up);

	Ciphertext<DCRTPoly> compDecrete(const Ciphertext<DCRTPoly> ciphertext, const int32_t bound);

	//Recrypt by decrypt-encrypt. Used for tests
	Ciphertext<DCRTPoly> fakeboot(const Ciphertext<DCRTPoly> ciphertext, const CryptoContext<DCRTPoly> cc, const KeyPair<DCRTPoly> keys, bool verbose = true);

	Ciphertext<DCRTPoly> boot(const Ciphertext<DCRTPoly> ciphertext, const CryptoContext<DCRTPoly> cc, const KeyPair<DCRTPoly> keys);


	// //----------------------------------------------------------------------------------
	// //   Sorting
	// //----------------------------------------------------------------------------------

	vector<Plaintext> maskPrecompute(const int32_t arrsize, const int32_t batchSize, const int32_t skewnormalizer, const CryptoContext<DCRTPoly> cc, const bool skew = true, const int32_t numk=0);

	vector<Plaintext> maskPrecompute_full(const int32_t size, const int32_t batchSize, const int32_t skewnormalizer, const CryptoContext<DCRTPoly> cc);

	// vector<Plaintext> maskPrecomputeForTopk(const usint size, const int32_t batchSize, const int32_t numk, const CryptoContext<DCRTPoly> cc);

	//boot1: eval boot at appropriate time, boot2: eval fakeboot at appropriate time
	Ciphertext<DCRTPoly> sort(const Ciphertext<DCRTPoly> ciphertext, const int32_t arrsize,  const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys, const bool boot1, const bool boot2, const Plaintext targetranks, const int32_t padding=0, const bool verbose = true);

	vector<Ciphertext<DCRTPoly>> sort_SIMD(const vector<Ciphertext<DCRTPoly>> ciphertext, const int32_t arrsize,  const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys, const int32_t numct, const bool boot1, const Plaintext targetranks, const int32_t padding=0, const bool multithread = true, const bool verbose = true);

	vector<Ciphertext<DCRTPoly>> rearrange_SIMD(const vector<Ciphertext<DCRTPoly>> ciphertext, const int32_t arrsize,  const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys,  const int32_t numk, const int32_t numct, const bool boot1, const bool multithread, const bool verbose = true);

	vector<Ciphertext<DCRTPoly>> bootPacked(const vector<Ciphertext<DCRTPoly>> ciphertext, const int32_t arrsize,  const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys,  const int32_t numk, const int32_t numct, const bool bootstrapping_int, const bool boot1, const bool multithread, const bool verbose = true);


	vector<Ciphertext<DCRTPoly>> encryptForTopk(const vector<double> vals, const int32_t size, const int32_t numk, const PublicKey<DCRTPoly> publicKey, CryptoContext<DCRTPoly> cc, const int32_t arrsize_init=0);

	Ciphertext<DCRTPoly> blocktopkS(const Ciphertext<DCRTPoly> ciphertext, const int32_t arrsize,  const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys, const int32_t numk, const bool boot1, const bool boot2, const bool verbose=true);

	vector<Ciphertext<DCRTPoly>> blocktopkM(const vector<Ciphertext<DCRTPoly>> ciphertexts, const int32_t size, const int32_t arrsize, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys, const int32_t numk, const bool boot1, const bool multithread = true, const int32_t mergecrit = 1);

	vector<int32_t> schedule_topk(const int32_t size, const int32_t numk, int32_t init_arrsize_bound = 64, int32_t init_numct_bound = 128, bool binarymerge = false, const int32_t boot_threshold = 32, const int32_t mergecrit = 1, const bool mute = false); 

	int32_t roundestimator(const int32_t size, const int32_t numk, int32_t init_arrsize_bound = 64, int32_t init_numct_bound = 128, const int32_t mergecrit = 1); 


	Ciphertext<DCRTPoly> topk_multithread_opt(const vector<Ciphertext<DCRTPoly>> ciphertexts, const int32_t size, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys, const int32_t numk, const bool boot1, const bool boot2, const bool multithread = true, const int32_t mergecrit = 1); 


	Ciphertext<DCRTPoly> nexus(const Ciphertext<DCRTPoly> ciphertext, const int32_t size, const int32_t bound, const usint scaleModSize, const KeyPair<DCRTPoly> keys);


	
}
#endif
