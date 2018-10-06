package fft

import "math"

// Bluestein algorithm for prime-n sized FFTs, i.e. pad input to a size m = 2**p ≥ 2*n-1
// for some integer p and use the circular convolution theorem to turn the size n DFT
// into two size m DFTs and a size m inverse DFT.

func bluesteinFwd(x []complex128) []complex128 {
	n := len(x)

	m := 1 << uint(math.Ilogb(float64(2*n-1)))
	if m < 2*n-1 {
		m <<= 1
	}

	w := make([]complex128, m)
	y := make([]complex128, m)
	copy(y, x)

	a0 := math.Pi / float64(n)
	w[0] = 1
	for i := 1; i < n; i++ {
		s, c := math.Sincos(a0 * float64(i*i))
		w[i] = complex(c, s)
		w[m-i] = complex(c, s)
		y[i] *= complex(c, -s)
	}

	y = stockham(y, 1)
	for i, ww := range stockham(w, 1) {
		y[i] *= ww
	}
	y = stockham(y, -1)

	for i := 0; i < n; i++ {
		y[i] *= complex(real(w[i])/float64(m), -imag(w[i])/float64(m))
	}

	return y[:n]
}

func bluesteinBwd(x []complex128) []complex128 {
	n := len(x)
	y := make([]complex128, n)
	for i, xi := range x {
		y[i] = complex(real(xi), -imag(xi))
	}
	y = bluesteinFwd(y)
	for i, yi := range y {
		y[i] = complex(real(yi)/float64(n), -imag(yi)/float64(n))
	}
	return y
}
