#ifndef WAVEMIN_H_
#define WAVEMIN_H_

#include "waveaux.h"

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_MSC_VER)
#pragma warning(disable : 4200)
#pragma warning(disable : 4996)
#endif

#ifndef cplx_type
#define cplx_type float
#endif

typedef struct cplx_t {
	cplx_type re;
	cplx_type im;
} cplx_data;

typedef struct wave_set* wave_object;

wave_object wave_init(char* wname);

struct wave_set{
	char wname[50];
	int filtlength;// When all filters are of the same length. [Matlab uses zero-padding to make all filters of the same length]
	int lpd_len;// Default filtlength = lpd_len = lpr_len = hpd_len = hpr_len
	int hpd_len;
	int lpr_len;
	int hpr_len;
	float *lpd;
	float *hpd;
	float *lpr;
	float *hpr;
	float params[0];
};

typedef struct wt_set* wt_object;

wt_object wt_init(wave_object wave, int siglength, int J);

struct wt_set{
	wave_object wave;
	char method[10];
	int siglength;// Length of the original signal.
	int outlength;// Length of the output DWT vector
	int lenlength;// Length of the Output Dimension Vector "length"
	int J; // Number of decomposition Levels
	int MaxIter;// Maximum Iterations J <= MaxIter
	int even;// even = 1 if signal is of even length. even = 0 otherwise
	char ext[10];// Type of Extension used - "per" or "sym"
	char cmethod[10]; // Convolution Method - "direct" or "FFT"

	int N; //
	int cfftset;
	int zpad;
	int length[102];
	float *output;
	float params[0];
};

	void dwt(wt_object wt, float *inp);

	void idwt_sym_direct(wt_object wt, float *dwtop);

	void wave_free(wave_object object);

	void wt_free(wt_object object);

#ifdef __cplusplus
}
#endif

#endif /* WAVEMIN_H_ */
