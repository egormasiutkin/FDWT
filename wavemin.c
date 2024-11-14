#include "wavemin.h"

wave_object wave_init(char* wname) {
	wave_object obj = NULL;
	int retval;
	retval = 0;

	if (wname != NULL) {
		retval = filtlength(wname);
	}

	obj = (wave_object)malloc(sizeof(struct wave_set) + sizeof(float) * 4 * retval);

	obj->filtlength = retval;
	obj->lpd_len = obj->hpd_len = obj->lpr_len = obj->hpr_len = obj->filtlength;
	strcpy(obj->wname, wname);
	if (wname != NULL) {
		filtcoef(wname, obj->params, obj->params + retval, obj->params + 2 * retval, obj->params + 3 * retval);
	}
	obj->lpd = &obj->params[0];
	obj->hpd = &obj->params[retval];
	obj->lpr = &obj->params[2 * retval];
	obj->hpr = &obj->params[3 * retval];
	return obj;
}

wt_object wt_init(wave_object wave, int siglength, int J) {
	int size, i, MaxIter;
	wt_object obj = NULL;

	size = wave->filtlength;

	if (J > 100) {
		J = 100;
	}

	MaxIter = wmaxiter(siglength, size);

	if (J > MaxIter) {
		J = MaxIter;
	}

	obj = (wt_object)malloc(sizeof(struct wt_set) + sizeof(float)* (siglength + 2 * J * (size + 1)));
	obj->outlength = siglength + 2 * J * (size + 1); // Default
	strcpy(obj->ext, "sym"); // Default

	obj->wave = wave;
	obj->siglength = siglength;
	obj->J = J;
	obj->MaxIter = MaxIter;
	strcpy(obj->method, "dwt");

	if (siglength % 2 == 0) {
		obj->even = 1;
	}
	else {
		obj->even = 0;
	}

	strcpy(obj->cmethod, "direct"); // Default
	obj->cfftset = 0;
	obj->lenlength = J + 2;
	obj->output = &obj->params[0];
	for (i = 0; i < siglength + 2 * J * (size + 1); ++i) {
		obj->params[i] = 0.0;
	}

	return obj;
}

static void wconv(wt_object wt, float *sig, int N, float *filt, int L, float *oup) {
	if (!strcmp(wt->cmethod, "direct")) {
		conv_direct(sig, N, filt, L, oup);
	}
	else {
		printf("Convolution Only accepts direct convolution");
		exit(-1);
	}
}

static void dwt_sym(wt_object wt, float *inp, int N, float *cA, int len_cA, float *cD, int len_cD) {
	int i, l, t, len_avg;

	len_avg = wt->wave->lpd_len;

	for (i = 0; i < len_cA; ++i) {
		t = 2 * i + 1;
		cA[i] = 0.0;
		cD[i] = 0.0;
		for (l = 0; l < len_avg; ++l) {
			if ((t - l) >= 0 && (t - l) < N) {
				cA[i] += wt->wave->lpd[l] * inp[t - l];
				cD[i] += wt->wave->hpd[l] * inp[t - l];
			}
			else if ((t - l) < 0) {
				cA[i] += wt->wave->lpd[l] * inp[-t + l - 1];
				cD[i] += wt->wave->hpd[l] * inp[-t + l - 1];
			}
			else if ((t - l) >= N) {
				cA[i] += wt->wave->lpd[l] * inp[2 * N - t + l - 1];
				cD[i] += wt->wave->hpd[l] * inp[2 * N - t + l - 1];
			}
		}
	}
}

void dwt(wt_object wt, float *inp) {
	int i, J, temp_len, iter, N, lp;
	int len_cA;
	float *orig, *orig2;

	temp_len = wt->siglength;
	J = wt->J;
	wt->length[J + 1] = temp_len;
	wt->outlength = 0;
	wt->zpad = 0;
	orig = (float*)malloc(sizeof(float)* temp_len);
	orig2 = (float*)malloc(sizeof(float)* temp_len);

	for (i = 0; i < wt->siglength; ++i) {
		orig[i] = inp[i];
	}

	if (wt->zpad == 1) {
		orig[temp_len - 1] = orig[temp_len - 2];
	}

	N = temp_len;
	lp = wt->wave->lpd_len;

	if (!strcmp(wt->ext, "sym")) {

		i = J;
		while (i > 0) {
			N = N + lp - 2;
			N = (int)ceil((float)N / 2.0);
			wt->length[i] = N;
			wt->outlength += wt->length[i];
			i--;
		}
		wt->length[0] = wt->length[1];
		wt->outlength += wt->length[0];
		N = wt->outlength;

		for (iter = 0; iter < J; ++iter) {
			len_cA = wt->length[J - iter];
			N -= len_cA;
			dwt_sym(wt, orig, temp_len, orig2, len_cA, wt->params + N, len_cA);
			temp_len = wt->length[J - iter];

			if (iter == J - 1) {
				for (i = 0; i < len_cA; ++i) {
					wt->params[i] = orig2[i];
				}
			}
			else {
				for (i = 0; i < len_cA; ++i) {
					orig[i] = orig2[i];
				}
			}
		}
	}
	else {
		printf("Signal extension error");
		exit(-1);
	}

	free(orig);
	free(orig2);
}

static void idwt_sym(wt_object wt, float *cA, int len_cA, float *cD, int len_cD, float *X) {
	int len_avg, i, l, m, n, t, v;

	len_avg = (wt->wave->lpr_len + wt->wave->hpr_len) / 2;
	m = -2;
	n = -1;

	for (v = 0; v < len_cA; ++v) {
		i = v;
		m += 2;
		n += 2;
		X[m] = 0.0;
		X[n] = 0.0;
		for (l = 0; l < len_avg / 2; ++l) {
			if ((cA[i - l] != 0) || (cD[i - l] != 0)) 
			{
				if (((i - l) >= 0) && ((i - l) < len_cA)) {
					t = 2 * l;
					X[m] += wt->wave->lpr[t] * cA[i - l] + wt->wave->hpr[t] * cD[i - l];
					X[n] += wt->wave->lpr[t + 1] * cA[i - l] + wt->wave->hpr[t + 1] * cD[i - l];
				}
			}			
		}
	}
}

void idwt_sym_direct(wt_object wt, float *dwtop) {
	int J, U, i, lf, N, iter, k;
	int app_len, det_len;
	float *X_lp, *out;

	J = wt->J;
	U = 2;
	app_len = wt->length[0];
	out = (float*)malloc(sizeof(float)* (wt->siglength + 1));

	app_len = wt->length[0];
	det_len = wt->length[1];
	N = 2 * wt->length[J] - 1;
	lf = (wt->wave->lpr_len + wt->wave->hpr_len) / 2;

	X_lp = (float*)malloc(sizeof(float)* (N + 2 * lf - 1));
	iter = app_len;

	for (i = 0; i < app_len; ++i) {
		out[i] = wt->output[i];
	}

	for (i = 0; i < J; ++i) {

		idwt_sym(wt, out, det_len, wt->output + iter, det_len, X_lp);
		for (k = lf - 2; k < 2 * det_len; ++k) {
			out[k - lf + 2] = X_lp[k];
		}

		iter += det_len;
		det_len = wt->length[i + 2];
	}

	free(X_lp);

	for (i = 0; i < wt->siglength; ++i) {
		dwtop[i] = out[i];
	}
	free(out);
}

void wave_free(wave_object object) {
	free(object);
}

void wt_free(wt_object object) {
	free(object);
}
