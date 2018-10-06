package fft

import "math"

// stockham returns the discrete Fourier transform or its unscaled inverse of x for
// s = 1 or -1 respectively. It is very efficient but limited to len(x) = 2**n.
//
// Stockham autosort algorithm based on
// "Computational Frameworks for the Fast Fourier Transform", Charles Van Loan, SIAM, 1992.
func stockham(x []complex128, s int) []complex128 {
	n := len(x)
	n2 := n >> 1

	tmp := make([]complex128, n)
	y := make([]complex128, n)
	copy(y, x)

	// Iterate log2(n/2) times.
	for r, l := n2, 1; r >= 1; r >>= 1 {
		y, tmp = tmp, y

		// Calculate twiddle factor w = exp(-s*iπ/l) components (wr, wi).
		wim, wre := math.Sincos(-float64(s) * math.Pi / float64(l))

		for j, wj := 0, complex(1, 0); j < l; j++ {
			jrs := j * (r << 1)
			for k, m := jrs, jrs>>1; k < jrs+r; k++ {
				t := wj * tmp[k+r]
				y[m] = tmp[k] + t
				y[m+n2] = tmp[k] - t
				m++
			}
			// Increment twiddle factor using w**j+1 = w**j * w.
			wj *= complex(wre, wim)
		}
		l <<= 1
	}
	return y
}
