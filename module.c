#include <linux/module.h>
#include <linux/types.h>
#include <crypto/akcipher.h>
#include <crypto/internal/akcipher.h>

/*
 * Wegen scatterlist: scatterwalk.h in include/?
 * In der Beschreibung von der Userspace-Lib steht, dass man für asymmetrische Verfahren einen Kernel Patch einspielen muss. Wenn das stimmt, muss ich darauf definitv achten und darauf in der Arbeit eingehen!
 */

struct walnut_ctx {
	union walnut_key {
		struct walnut_pub_key *pub_key;
		int8_t *priv_key;
	} key;
};

struct walnut_pub_key {
	uint8_t *matrix;
	uint8_t *permutation;
};

// encrypt und decrypt funktioniert gar nicht? Im Paper ist nur signatur und verifizierung spezifiziert, aber für die andern zwei müssten die operationen ja mit den umgekehrten schlüsseln verwendet werden, was aber nicht so einfach geht weil die ja komplett anders sind... also wahrscheinlich ist das gar nicht möglich!
static int walnut_enc(struct akcipher_request *req)
{
	struct crypto_akcipher *tfm = crypto_akcipher_reqtfm(req);
	const struct walnut_pub_key *pub_key = akcipher_tfm_ctx(tfm);
	// so was wie
	// return -1
	// und
	// print Operation not supported...

	return -1;
}

static int walnut_dec(struct akcipher_request *req)
{
	struct crypto_akcipher *tfm = crypto_akcipher_reqtfm(req);
	const int8_t *priv_key = akcipher_tfm_ctx(tfm);

	return -1;
}

// Signatur braucht eigentlich Permutation des public keys, aber die kann man auch aus dem private key gewinnen. Vielleicht ist es so einfacher, dann brauch ich hier nur den private key mit zu übergeben.
static int walnut_sign(struct akcipher_request *req)
{
	struct crypto_akcipher *tfm = crypto_akcipher_reqtfm(req);
	const int8_t *priv_key = akcipher_tfm_ctx(tfm);

	/* signatur anschmeißen
	   Daten sind in req->src mit req->src_len
	   resultat in req->dst mit req->dst_len schreiben
	   src und dst sind scatterlists... Wie kann ich die befüllen?
	   */

	return 0;
}

static int walnut_verify(struct akcipher_request *req)
{
	struct crypto_akcipher *tfm = crypto_akcipher_reqtfm(req);
	const struct walnut_pub_key *pub_key = akcipher_tfm_ctx(tfm);

	return 0;
}

static int walnut_set_pub_key(struct crypto_akcipher *tfm, const void *key,
	unsigned int keylen)
{
	struct walnut_pub_key *pub_key = akcipher_tfm_ctx(tfm);

	// vllt. key verifizieren? dafür brauch ich aber auch N, wie wird das
	// überhaupt gemacht? Die Konfiguration?

	// check for maximum length overflow of key...

	//memcpy key to pub_key
	
	return 0;
}

static int walnut_set_priv_key(struct crypto_akcipher *tfm, const void *key,
	unsigned int keylen)
{
	int8_t *priv_key = akcipher_tfm_ctx(tfm);
	// s. pub_key

	return 0;
}

static unsigned int walnut_max_size(struct crypto_akcipher *tfm)
{
	// how can we figure out the max size? maybe take the worst case length
	// of the BKL normal form, as if when there was no dehornoy reduction
	// afterwards. because BKL increases the size of the signature afaik.
	return 0;
}

static void walnut_exit_tfm(struct crypto_akcipher *tfm)
{
	// free whatever needs to be freed here...
	// the pub_key struct, for example.
	
	struct walnut_pub_key *pub_key = akcipher_tfm_ctx(tfm);
	kfree(pub_key);
}

static struct akcipher_alg walnut_dsa = {
	.encrypt = walnut_enc,
	.decrypt = walnut_dec,
	.sign = walnut_sign,
	.verify = walnut_verify,
	.set_pub_key = walnut_set_pub_key,
	.set_priv_key = walnut_set_priv_key,
	.max_size = walnut_max_size,
	.exit = walnut_exit_tfm,
	.base = {
		.cra_name = "walnut_dsa",
		.cra_driver_name = "walnut_dsa",
		.cra_priority = 100,
		.cra_module = THIS_MODULE,
		.cra_ctxsize = sizeof(struct walnut_ctx)
	}
};

static int walnut_init(void)
{
	printk(KERN_ALERT "Hello, world\n");
	int err = crypto_register_akcipher(&walnut_dsa);
	if (err)
		return err;
	return 0;
}

static void walnut_exit(void)
{
	crypto_unregister_akcipher(&walnut_dsa);
}

module_init(walnut_init);
module_exit(walnut_exit);
MODULE_ALIAS_CRYPTO("walnut_dsa");
MODULE_LICENSE("GPL");
MODULE_AUTHOR("Stefan Seil <stefan.seil@campus.lmu.de>");
MODULE_DESCRIPTION("WalnutDSA C implementation");

