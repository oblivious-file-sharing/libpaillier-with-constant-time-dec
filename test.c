#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include "paillier.h"

int main(){
	int modulus_bits = 2048;
	paillier_pubkey_t *ppk;
	paillier_prvkey_t *psk;

	paillier_keygen(modulus_bits, &ppk, &psk, paillier_get_rand_devurandom);

	paillier_plaintext_t *pt;
	char test[100] = "The Paillier encryption and decryption is correct.";
	pt = paillier_plaintext_from_bytes(test, 100);

	paillier_ciphertext_t *ct;
	ct = paillier_enc(NULL, ppk, pt, paillier_get_rand_devurandom);

	printf("Start the rerandomization.\n");
	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);

	for(int i = 0; i < 1000; i++){
		paillier_rerand(ct, ppk, paillier_get_rand_devurandom);
	}
	clock_gettime(CLOCK_MONOTONIC, &end);
	printf("End the rerandomization.\n");
	printf("Paillier rerandomization per 2048 bits: %lf\n", (end.tv_sec - start.tv_sec + (end.tv_nsec - start.tv_nsec) / 1000000000.) / 1000);

	paillier_plaintext_t *pt_dec;
	pt_dec = paillier_dec(NULL, ppk, psk, ct);

	char *test_dec = paillier_plaintext_to_bytes(100, pt_dec);

	printf("%s\n", test_dec);

	paillier_freeplaintext(pt_dec);
	paillier_freeplaintext(pt);

	paillier_freeciphertext(ct);

	paillier_freepubkey(ppk);
	paillier_freeprvkey(psk);

	return 0;
}
