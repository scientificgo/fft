package fft

import "math"

const (
	sqrt3 = 1.7320508075688772935274463415058723669428052538103806280558 // https://oeis.org/A002194
	w31   = (-1 - sqrt3*1i) / 2                                          // exp(-2πi*1/3)
	w32   = (-1 + sqrt3*1i) / 2                                          // exp(-2πi*2/3)
)

func radix3(x []complex128, s int) []complex128 {
	const r = 3
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
				t0, t1, t2 := y[i+j], y[i+j+mr], y[i+j+2*mr]

				// Apply twiddle factors w**(i+k) for 1 ≤ k < r.
				t1 *= wi
				t2 *= wi * wi

				// Transform points using r-point DFT.
				y[i+j] += t1 + t2
				if s > 0 {
					y[i+j+mr], y[i+j+2*mr] =
						t0+t1*w31+t2*w32,
						t0+t1*w32+t2*w31
				} else {
					// 1/w31 = w32, etc.
					y[i+j+mr], y[i+j+2*mr] =
						t0+t1*w32+t2*w31,
						t0+t1*w31+t2*w32
				}
			}
			wi *= w
		}
	}
	return y
}
