package fft

import "math"

const (
	w71 = 0.62348980185873353052500488400423981063227473089640210536 - 0.78183148246802980870844452667405775023233451870868752898i // exp(-2πi*1/7)
	w72 = -0.2225209339563144042889025644967947594663555687645449553 - 0.9749279121818236070181316829939312172327858006199974376i  // exp(-2πi*2/7)
	w73 = -0.9009688679024191262361023195074450511659191621318571500 - 0.4338837391175581204757683328483587546099907277874598764i  // exp(-2πi*3/7)
	w74 = -0.9009688679024191262361023195074450511659191621318571500 + 0.4338837391175581204757683328483587546099907277874598764i  // exp(-2πi*4/7)
	w75 = -0.2225209339563144042889025644967947594663555687645449553 + 0.9749279121818236070181316829939312172327858006199974376i  // exp(-2πi*5/7)
	w76 = 0.62348980185873353052500488400423981063227473089640210536 + 0.78183148246802980870844452667405775023233451870868752898i // exp(-2πi*6/7)
)

func radix7(x []complex128, s int) []complex128 {
	const r = 7
	n := len(x)

	// Copy input into new storage.
	y := make([]complex128, n)
	copy(y, x)

	// Reorder input using base-r digit reversal permutation.
	for i, j := 0, 0; i < n-1; i++ {
		if i < j {
			y[i], y[j] = y[j], y[i]
		}
		k := (r - 1) * n / r
		for k <= j {
			j -= k
			k /= r
		}
		j += k / (r - 1)
	}

	for m := r; m <= n; m *= r {
		// Calculate twiddle factor.
		wim, wre := math.Sincos(-2 * math.Pi * float64(s) / float64(m))
		w := complex(wre, wim)

		mr := m / r
		for i, wi := 0, complex(1, 0); i < mr; i++ {
			for j := 0; j < n; j += m {
				// Retrieve subset of points.
				t0, t1, t2, t3, t4, t5, t6 := y[i+j], y[i+j+mr], y[i+j+2*mr], y[i+j+3*mr], y[i+j+4*mr], y[i+j+5*mr], y[i+j+6*mr]

				// Apply twiddle factors w**(i+k) for 1 ≤ k < r.
				t1 *= wi
				t2 *= wi * wi
				t3 *= wi * wi * wi
				t4 *= wi * wi * wi * wi
				t5 *= wi * wi * wi * wi * wi
				t6 *= wi * wi * wi * wi * wi * wi

				// Transform points using r-point DFT.
				y[i+j] += t1 + t2 + t3 + t4 + t5 + t6
				if s > 0 {
					y[i+j+mr], y[i+j+2*mr], y[i+j+3*mr], y[i+j+4*mr], y[i+j+5*mr], y[i+j+6*mr] =
						t0+t1*w71+t2*w72+t3*w73+t4*w74+t5*w75+t6*w76,
						t0+t1*w72+t2*w74+t3*w76+t4*w71+t5*w73+t6*w75,
						t0+t1*w73+t2*w76+t3*w72+t4*w75+t5*w71+t6*w74,
						t0+t1*w74+t2*w71+t3*w75+t4*w72+t5*w76+t6*w73,
						t0+t1*w75+t2*w73+t3*w71+t4*w76+t5*w74+t6*w72,
						t0+t1*w76+t2*w75+t3*w74+t4*w73+t5*w72+t6*w71
				} else {
					// 1/w71 = w76, etc.
					y[i+j+mr], y[i+j+2*mr], y[i+j+3*mr], y[i+j+4*mr], y[i+j+5*mr], y[i+j+6*mr] =
						t0+t1*w76+t2*w75+t3*w74+t4*w73+t5*w72+t6*w71,
						t0+t1*w75+t2*w73+t3*w71+t4*w76+t5*w74+t6*w72,
						t0+t1*w74+t2*w71+t3*w75+t4*w72+t5*w76+t6*w73,
						t0+t1*w73+t2*w76+t3*w72+t4*w75+t5*w71+t6*w74,
						t0+t1*w72+t2*w74+t3*w76+t4*w71+t5*w73+t6*w75,
						t0+t1*w71+t2*w72+t3*w73+t4*w74+t5*w75+t6*w76
				}
			}
			wi *= w
		}
	}
	return y
}
