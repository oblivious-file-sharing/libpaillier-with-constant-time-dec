/*
	libpaillier - A library implementing the Paillier cryptosystem.

	Copyright (C) 2006 SRI International.

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful, but
	WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
	General Public License for more details.
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "paillier.h"

void paillier_inline_mpz_recalloc (mpz_t rp, mp_size_t target) {
	mp_size_t original = mpz_size(rp);
	
	if (original < target) {
		mp_limb_t *tmp = mpz_limbs_modify(rp, target);
		mpn_zero(&tmp[original], target - original);
		mpz_limbs_finish(rp, target);
	}
}

void init_rand (gmp_randstate_t rand, paillier_get_rand_t get_rand, int bytes) {
	void *buf;
	mpz_t s;
	
	buf = malloc(bytes);
	get_rand(buf, bytes);
	
	gmp_randinit_default(rand);
	mpz_init(s);
	mpz_import(s, bytes, 1, 1, 0, 0, buf);
	gmp_randseed(rand, s);
	mpz_clear(s);
	
	free(buf);
}

void complete_pubkey (paillier_pubkey_t *pub) {
	mpz_mul(pub->n_squared, pub->n, pub->n);
	mpz_add_ui(pub->n_plusone, pub->n, 1);
}

void complete_prvkey (paillier_prvkey_t *prv, paillier_pubkey_t *pub) {
	mpz_powm(prv->x, pub->n_plusone, prv->lambda, pub->n_squared);
	mpz_sub_ui(prv->x, prv->x, 1);
	mpz_div(prv->x, prv->x, pub->n);
	mpz_invert(prv->x, prv->x, pub->n);
}

void paillier_keygen (int modulusbits,
                      paillier_pubkey_t **pub,
                      paillier_prvkey_t **prv,
                      paillier_get_rand_t get_rand) {
	mpz_t p;
	mpz_t q;
	gmp_randstate_t rand;
	
	/* allocate the new key structures */
	
	*pub = (paillier_pubkey_t *) malloc(sizeof(paillier_pubkey_t));
	*prv = (paillier_prvkey_t *) malloc(sizeof(paillier_prvkey_t));
	
	/* initialize our integers */
	
	mpz_init((*pub)->n);
	mpz_init((*pub)->n_squared);
	mpz_init((*pub)->n_plusone);
	mpz_init((*prv)->lambda);
	mpz_init((*prv)->x);
	mpz_init(p);
	mpz_init(q);
	
	/* pick random (modulusbits/2)-bit primes p and q */
	init_rand(rand, get_rand, modulusbits / 8 + 1);
	do {
		do
			mpz_urandomb(p, rand, modulusbits / 2);
		while (!mpz_probab_prime_p(p, 10));
		
		do
			mpz_urandomb(q, rand, modulusbits / 2);
		while (!mpz_probab_prime_p(q, 10));
		
		/* compute the public modulus n = p q */
		
		mpz_mul((*pub)->n, p, q);
	} while (!mpz_tstbit((*pub)->n, modulusbits - 1));
	complete_pubkey(*pub);
	(*pub)->bits = modulusbits;
	
	/* compute the private key lambda = lcm(p-1,q-1) */
	
	mpz_sub_ui(p, p, 1);
	mpz_sub_ui(q, q, 1);
	mpz_lcm((*prv)->lambda, p, q);
	complete_prvkey(*prv, *pub);
	
	/* clear temporary integers and randstate */
	
	mpz_clear(p);
	mpz_clear(q);
	gmp_randclear(rand);
}

paillier_ciphertext_t *paillier_enc (paillier_ciphertext_t *res,
                                     paillier_pubkey_t *pub,
                                     paillier_plaintext_t *pt,
                                     paillier_get_rand_t get_rand) {
	mpz_t r;
	gmp_randstate_t rand;
	mpz_t x;
	
	/* pick random blinding factor */
	
	mpz_init(r);
	init_rand(rand, get_rand, pub->bits / 8 + 1);
	do
		mpz_urandomb(r, rand, pub->bits);
	while (mpz_cmp(r, pub->n) >= 0);
	
	/* compute ciphertext */
	
	if (!res) {
		res = (paillier_ciphertext_t *) malloc(sizeof(paillier_ciphertext_t));
		mpz_init(res->c);
	}
	
	mpz_init(x);
	mpz_powm(res->c, pub->n_plusone, pt->m, pub->n_squared);
	mpz_powm(x, r, pub->n, pub->n_squared);
	
	mpz_mul(res->c, res->c, x);
	mpz_mod(res->c, res->c, pub->n_squared);
	
	mpz_clear(x);
	mpz_clear(r);
	gmp_randclear(rand);
	
	return res;
}

/* add: rerandomization */
paillier_ciphertext_t *paillier_rerand (paillier_ciphertext_t *ct,
                                        paillier_pubkey_t *pub,
                                        paillier_get_rand_t get_rand) {
	mpz_t r;
	gmp_randstate_t rand;
	
	/* pick random blinding factor */
	mpz_init(r);
	init_rand(rand, get_rand, pub->bits / 8 + 1);
	do
		mpz_urandomb(r, rand, pub->bits);
	while (mpz_cmp(r, pub->n) >= 0);
	
	mpz_powm(r, r, pub->n, pub->n_squared);
	
	mpz_mul(ct->c, ct->c, r);
	mpz_mod(ct->c, ct->c, pub->n_squared);
	
	mpz_clear(r);
	gmp_randclear(rand);
	
	return ct;
}

paillier_plaintext_t *paillier_dec (paillier_plaintext_t *res,
                                    paillier_pubkey_t *pub,
                                    paillier_prvkey_t *prv,
                                    paillier_ciphertext_t *ct) {
	if (!res) {
		res = (paillier_plaintext_t *) malloc(sizeof(paillier_plaintext_t));
		mpz_init(res->m);
	}
	
	int nbits = pub->bits;
	int nlimbs = (nbits + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
	
	int nbits_2 = nbits * 2;
	int nlimbs_2 = (nbits_2 + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
	
	mpz_t tmp_1;
	mpz_init_set_ui(tmp_1, 1);
	
	mpz_t tmp;
	mpz_init2(tmp, nbits_2);
	
	paillier_inline_mpz_recalloc(tmp_1, nlimbs);
	
	const mp_limb_t (*c_ro);
	const mp_limb_t (*n_squared_ro);
	const mp_limb_t (*n_ro);
	const mp_limb_t (*lambda_ro);
	mp_limb_t *m_w;
	const mp_limb_t (*x_ro);
	const mp_limb_t (*m_ro);
	const mp_limb_t (*tmp_1_ro);
	const mp_limb_t (*tmp_ro);
	mp_limb_t *tmp_w;
	
	{
		// res->m = (ct->c)^priv mod n^2
		
		paillier_inline_mpz_recalloc(res->m, nlimbs_2);
		paillier_inline_mpz_recalloc(ct->c, nlimbs_2);
		paillier_inline_mpz_recalloc(pub->n_squared, nlimbs_2);
		paillier_inline_mpz_recalloc(prv->lambda, (nbits + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS);
		
		mp_limb_t scratch_powm[mpn_sec_powm_itch(nlimbs_2, nbits, nlimbs_2)];
		
		c_ro = mpz_limbs_read(ct->c);
		n_squared_ro = mpz_limbs_read(pub->n_squared);
		lambda_ro = mpz_limbs_read(prv->lambda);
		m_w = mpz_limbs_write(res->m, nlimbs_2);
		
		mpn_sec_powm(m_w, c_ro, nlimbs_2, lambda_ro, nbits, n_squared_ro, nlimbs_2, scratch_powm);
		
		mpz_limbs_finish(res->m, nlimbs_2);
	}
	
	{
		// tmp = res->m - 1 mod n^2
		
		mp_limb_t scratch_sub[mpn_sec_sub_1_itch(nlimbs_2)];
		
		tmp_1_ro = mpz_limbs_read(tmp_1);
		m_ro = mpz_limbs_read(res->m);
		tmp_w = mpz_limbs_write(tmp, nlimbs_2);
		
		mpn_sec_sub_1(tmp_w, m_ro, nlimbs_2, tmp_1_ro[0], scratch_sub);
		mpz_limbs_finish(tmp, nlimbs_2);
	}
	
	{
		// res->m = tmp / n  mod n
		paillier_inline_mpz_recalloc(pub->n, nlimbs);
		
		mp_limb_t scratch_div_qr[mpn_sec_div_qr_itch(nlimbs_2, nlimbs)];
		
		n_ro = mpz_limbs_read(pub->n);
		m_w = mpz_limbs_write(res->m, nlimbs_2 - nlimbs);
		tmp_w = mpz_limbs_modify(tmp, nlimbs_2);
		
		mpn_sec_div_qr(m_w, tmp_w, nlimbs_2, n_ro, nlimbs, scratch_div_qr);
		mpz_limbs_finish(res->m, nlimbs_2 - nlimbs);
		mpz_limbs_finish(tmp, nlimbs_2);
	}
	
	{
		//tmp = res->m * prv->x
		
		paillier_inline_mpz_recalloc(prv->x, nlimbs);
		
		tmp_w = mpz_limbs_write(tmp, nlimbs_2);
		x_ro = mpz_limbs_read(prv->x);
		m_ro = mpz_limbs_read(res->m);
		
		mp_limb_t scratch_mul[mpn_sec_mul_itch(nlimbs, nlimbs)];
		
		mpn_sec_mul(tmp_w, m_ro, nlimbs, x_ro, nlimbs, scratch_mul);
		
		mpz_limbs_finish(tmp, nlimbs_2);
	}
	
	{
		// tmp = tmp mod n
		
		tmp_w = mpz_limbs_modify(tmp, nlimbs_2);
		
		mp_limb_t scratch_div_r[mpn_sec_div_r_itch(nlimbs_2, nlimbs)];
		mpn_sec_div_r(tmp_w, nlimbs_2, n_ro, nlimbs, scratch_div_r);
		
		mpz_limbs_finish(tmp, nlimbs);
	}
	
	{
		// res->m = tmp
		
		m_w = mpz_limbs_write(res->m, nlimbs);
		tmp_ro = mpz_limbs_read(tmp);
		
		mpn_copyi(m_w, tmp_ro, nlimbs);
		
		mpz_limbs_finish(res->m, nlimbs);
	}
	
	mpz_clears(tmp, tmp_1, NULL);
	
	return res;
}

void paillier_mul (paillier_pubkey_t *pub,
                   paillier_ciphertext_t *res,
                   paillier_ciphertext_t *ct0,
                   paillier_ciphertext_t *ct1) {
	mpz_mul(res->c, ct0->c, ct1->c);
	mpz_mod(res->c, res->c, pub->n_squared);
}

void paillier_exp (paillier_pubkey_t *pub,
                   paillier_ciphertext_t *res,
                   paillier_ciphertext_t *ct,
                   paillier_plaintext_t *pt) {
	mpz_powm(res->c, ct->c, pt->m, pub->n_squared);
}

paillier_plaintext_t *paillier_plaintext_from_ui (unsigned long int x) {
	paillier_plaintext_t *pt;
	
	pt = (paillier_plaintext_t *) malloc(sizeof(paillier_plaintext_t));
	mpz_init_set_ui(pt->m, x);
	
	return pt;
}

paillier_plaintext_t *paillier_plaintext_from_bytes (void *m, int len) {
	paillier_plaintext_t *pt;
	
	pt = (paillier_plaintext_t *) malloc(sizeof(paillier_plaintext_t));
	mpz_init(pt->m);
	mpz_import(pt->m, len, 1, 1, 0, 0, m);
	
	return pt;
}

void *paillier_plaintext_to_bytes (int len,
                                   paillier_plaintext_t *pt) {
	void *buf0;
	void *buf1;
	size_t written;
	
	buf0 = mpz_export(0, &written, 1, 1, 0, 0, pt->m);
	
	if (written == len)
		return buf0;
	
	buf1 = malloc(len);
	memset(buf1, 0, len);
	
	if (written == 0)
		/* no need to copy anything, pt->m = 0 and buf0 was not allocated */
		return buf1;
	else if (written < len)
		/* pad with leading zeros */
		memcpy(buf1 + (len - written), buf0, written);
	else
		/* truncate leading garbage */
		memcpy(buf1, buf0 + (written - len), len);
	
	free(buf0);
	
	return buf1;
}

paillier_plaintext_t *paillier_plaintext_from_str (char *str) {
	return paillier_plaintext_from_bytes(str, strlen(str));
}

char *paillier_plaintext_to_str (paillier_plaintext_t *pt) {
	char *buf;
	size_t len;
	
	buf = (char *) mpz_export(0, &len, 1, 1, 0, 0, pt->m);
	buf = (char *) realloc(buf, len + 1);
	buf[len] = 0;
	
	return buf;
}

paillier_ciphertext_t *paillier_ciphertext_from_bytes (void *c, int len) {
	paillier_ciphertext_t *ct;
	
	ct = (paillier_ciphertext_t *) malloc(sizeof(paillier_ciphertext_t));
	mpz_init(ct->c);
	mpz_import(ct->c, len, 1, 1, 0, 0, c);
	
	return ct;
}

void *paillier_ciphertext_to_bytes (int len,
                                    paillier_ciphertext_t *ct) {
	void *buf;
	int cur_len;
	
	cur_len = mpz_sizeinbase(ct->c, 2);
	cur_len = PAILLIER_BITS_TO_BYTES(cur_len);
	buf = malloc(len);
	memset(buf, 0, len);
	mpz_export(buf + (len - cur_len), 0, 1, 1, 0, 0, ct->c);
	
	return buf;
}

char *paillier_pubkey_to_hex (paillier_pubkey_t *pub) {
	return mpz_get_str(0, 16, pub->n);
}

char *paillier_prvkey_to_hex (paillier_prvkey_t *prv) {
	return mpz_get_str(0, 16, prv->lambda);
}

paillier_pubkey_t *paillier_pubkey_from_hex (char *str) {
	paillier_pubkey_t *pub;
	
	pub = (paillier_pubkey_t *) malloc(sizeof(paillier_pubkey_t));
	mpz_init_set_str(pub->n, str, 16);
	pub->bits = mpz_sizeinbase(pub->n, 2);
	mpz_init(pub->n_squared);
	mpz_init(pub->n_plusone);
	complete_pubkey(pub);
	
	return pub;
}

paillier_prvkey_t *paillier_prvkey_from_hex (char *str, paillier_pubkey_t *pub) {
	paillier_prvkey_t *prv;
	
	prv = (paillier_prvkey_t *) malloc(sizeof(paillier_prvkey_t));
	mpz_init_set_str(prv->lambda, str, 16);
	mpz_init(prv->x);
	complete_prvkey(prv, pub);
	
	return prv;
}

void paillier_freepubkey (paillier_pubkey_t *pub) {
	mpz_clear(pub->n);
	mpz_clear(pub->n_squared);
	mpz_clear(pub->n_plusone);
	free(pub);
}

void paillier_freeprvkey (paillier_prvkey_t *prv) {
	mpz_clear(prv->lambda);
	mpz_clear(prv->x);
	free(prv);
}

void paillier_freeplaintext (paillier_plaintext_t *pt) {
	mpz_clear(pt->m);
	free(pt);
}

void paillier_freeciphertext (paillier_ciphertext_t *ct) {
	mpz_clear(ct->c);
	free(ct);
}

void paillier_get_rand_file (void *buf, int len, char *file) {
	FILE *fp;
	void *p;
	
	fp = fopen(file, "r");
	
	p = buf;
	while (len) {
		size_t s;
		s = fread(p, 1, len, fp);
		p += s;
		len -= s;
	}
	
	fclose(fp);
}

void paillier_get_rand_devrandom (void *buf, int len) {
	paillier_get_rand_file(buf, len, "/dev/random");
}

void paillier_get_rand_devurandom (void *buf, int len) {
	paillier_get_rand_file(buf, len, "/dev/urandom");
}

paillier_ciphertext_t *paillier_create_enc_zero () {
	paillier_ciphertext_t *ct;
	
	/* make a NON-RERANDOMIZED encryption of zero for the purposes of
		 homomorphic computation */
	
	/* note that this is just the number 1 */
	
	ct = (paillier_ciphertext_t *) malloc(sizeof(paillier_ciphertext_t));
	mpz_init_set_ui(ct->c, 1);
	
	return ct;
}
