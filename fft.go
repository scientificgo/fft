// Copyright (c) 2018, Jack Parkinson. All rights reserved.
// Use of this source code is governed by the BSD 3-Clause
// license that can be found in the LICENSE file.

package fft

import "math"

// Fft returns a discrete Fourier transform of x.
// It does not check x for NaN or Inf values; these checks should be done separately.
func Fft(x []complex128) []complex128 {
	n := len(x)
	if n < 2 {
		return x
	}
	var res []complex128
	switch p := base(n); p {
	case 2:
		res = fftStockham(x, 1)
	case 3:
		res = fftRad3(x, 1)
	case 5:
		res = fftRad5(x, 1)
	case 6:
		res = fftRad6(x, 1)
	case 7:
		res = fftRad7(x, 1)
	default:
		res = fftBluestein(x)
	}

	return res
}

// Ifft returns the inverse discrete Fourier transform of x.
// It does not check x for NaN or Inf values; these checks should be done separately.
func Ifft(x []complex128) []complex128 {
	n := len(x)

	if n < 2 {
		return x
	}

	var res []complex128
	var scaled bool
	switch p := base(n); p {
	case 2:
		res = fftStockham(x, -1)
	case 3:
		res = fftRad3(x, -1)
	case 5:
		res = fftRad5(x, -1)
	case 6:
		res = fftRad6(x, -1)
	case 7:
		res = fftRad7(x, -1)
	default:
		res = ifftBluestein(x)
		scaled = true
	}

	if !scaled {
		for i := 0; i < n; i++ {
			res[i] *= complex(1/float64(n), 0)
		}
	}

	return res
}

// base returns the smallest integer (up to 7, the max CT radix) such that n = base**m
// for some integer m, or n if there is none.
func base(n int) int {
	if n&(n-1) == 0 {
		return 2
	}

	for i := 3; i <= 7; i++ {
		n2 := n
		for n2%i == 0 {
			n2 /= i
		}
		if n2 == 1 {
			return i
		}
	}

	return n
}

// fftStockham returns the discrete Fourier transform or its unscaled inverse of x for
// flag = 1 or -1 respectively. It is very efficient but limited to len(x) = 2**n.
//
// Stockham autosort algorithm based on
// "Computational Frameworks for the Fast Fourier Transform", Charles Van Loan, SIAM, 1992.
func fftStockham(x []complex128, flag int) []complex128 {
	n := len(x)
	n2 := n >> 1

	tmp := make([]complex128, n)
	y := make([]complex128, n)
	copy(y, x)

	// Iterate log2(n/2) times.
	for r, l := n2, 1; r >= 1; r >>= 1 {
		y, tmp = tmp, y

		// Calculate twiddle factor w = exp(-flag*iπ/l) components (wr, wi).
		wi, wr := math.Sincos(-float64(flag) * math.Pi / float64(l))

		for j, wj := 0, complex(1, 0); j < l; j++ {
			jrs := j * (r << 1)
			for k, m := jrs, jrs>>1; k < jrs+r; k++ {
				t := wj * tmp[k+r]
				y[m] = tmp[k] + t
				y[m+n2] = tmp[k] - t
				m++
			}
			// Increment twiddle factor using w**j+1 = w**j * w.
			wj *= complex(wr, wi)
		}
		l <<= 1
	}
	return y
}

// The following fftRadX functions return the fft or its unscaled inverse for flag = 1 or -1
// respectively using a recursive Cooley-Tukey algorithm with radix-X.
//
// TODO: The functions are O(N logN) but can almost certainly be improved. For example, rewriting
// in a non-recusrive manner would likely be more efficient.

func fftRad3(x []complex128, flag int) []complex128 {
	n := len(x)
	if n < 2 {
		return x
	}

	radix := 3
	np := n / radix

	res := make([]complex128, n)
	for i := 0; i < np; i++ {
		for j := 0; j < radix; j++ {
			res[i+j*np] = x[radix*i+j]
		}
	}

	f0 := fftRad3(res[:np], flag)
	f1 := fftRad3(res[np:2*np], flag)
	f2 := fftRad3(res[2*np:], flag)

	res = make([]complex128, n)

	wexp0 := -2 * math.Pi * float64(flag) / float64(n)
	wi, wr := math.Sincos(wexp0)
	w0 := complex(1, 0)
	for i := 0; i < np; i++ {
		w0 *= complex(wr, wi)
	}

	w := complex(1, 0)
	for j := 0; j < radix; j++ {
		wj := w
		for k := 0; k < np; k++ {
			res[k+j*np] = f0[k] +
				(f1[k]+f2[k]*wj)*wj
			wj *= complex(wr, wi)
		}
		w *= w0
	}
	return res
}

func fftRad5(x []complex128, flag int) []complex128 {
	n := len(x)
	if n < 2 {
		return x
	}

	radix := 5
	np := n / radix

	res := make([]complex128, n)
	for i := 0; i < np; i++ {
		for j := 0; j < radix; j++ {
			res[i+j*np] = x[radix*i+j]
		}
	}

	f0 := fftRad5(res[:np], flag)
	f1 := fftRad5(res[np:2*np], flag)
	f2 := fftRad5(res[2*np:3*np], flag)
	f3 := fftRad5(res[3*np:4*np], flag)
	f4 := fftRad5(res[4*np:], flag)

	res = make([]complex128, n)

	wexp0 := -2 * math.Pi * float64(flag) / float64(n)
	wi, wr := math.Sincos(wexp0)
	w0 := complex(1, 0)
	for i := 0; i < np; i++ {
		w0 *= complex(wr, wi)
	}

	w := complex(1, 0)
	for j := 0; j < radix; j++ {
		wj := w
		for k := 0; k < np; k++ {
			res[k+j*np] = f0[k] +
				(f1[k]+(f2[k]+(f3[k]+f4[k]*wj)*wj)*wj)*wj
			wj *= complex(wr, wi)
		}
		w *= w0
	}
	return res
}

func fftRad6(x []complex128, flag int) []complex128 {
	n := len(x)
	if n < 2 {
		return x
	}

	radix := 6
	np := n / radix

	res := make([]complex128, n)
	for i := 0; i < np; i++ {
		for j := 0; j < radix; j++ {
			res[i+j*np] = x[radix*i+j]
		}
	}

	f0 := fftRad6(res[:np], flag)
	f1 := fftRad6(res[np:2*np], flag)
	f2 := fftRad6(res[2*np:3*np], flag)
	f3 := fftRad6(res[3*np:4*np], flag)
	f4 := fftRad6(res[4*np:5*np], flag)
	f5 := fftRad6(res[5*np:], flag)

	res = make([]complex128, n)

	wexp0 := -2 * math.Pi * float64(flag) / float64(n)
	wi, wr := math.Sincos(wexp0)
	w0 := complex(1, 0)
	for i := 0; i < np; i++ {
		w0 *= complex(wr, wi)
	}

	w := complex(1, 0)
	for j := 0; j < radix; j++ {
		wj := w
		for k := 0; k < np; k++ {
			res[k+j*np] = f0[k] +
				(f1[k]+(f2[k]+(f3[k]+(f4[k]+f5[k]*wj)*wj)*wj)*wj)*wj
			wj *= complex(wr, wi)
		}
		w *= w0
	}
	return res
}

func fftRad7(x []complex128, flag int) []complex128 {
	n := len(x)
	if n < 2 {
		return x
	}

	radix := 7
	np := n / radix

	res := make([]complex128, n)
	for i := 0; i < np; i++ {
		for j := 0; j < radix; j++ {
			res[i+j*np] = x[radix*i+j]
		}
	}

	f0 := fftRad7(res[:np], flag)
	f1 := fftRad7(res[np:2*np], flag)
	f2 := fftRad7(res[2*np:3*np], flag)
	f3 := fftRad7(res[3*np:4*np], flag)
	f4 := fftRad7(res[4*np:5*np], flag)
	f5 := fftRad7(res[5*np:6*np], flag)
	f6 := fftRad7(res[6*np:], flag)

	res = make([]complex128, n)

	wexp0 := -2 * math.Pi * float64(flag) / float64(n)
	wi, wr := math.Sincos(wexp0)
	w0 := complex(1, 0)
	for i := 0; i < np; i++ {
		w0 *= complex(wr, wi)
	}

	w := complex(1, 0)
	for j := 0; j < radix; j++ {
		wj := w
		for k := 0; k < np; k++ {
			res[k+j*np] = f0[k] +
				(f1[k]+(f2[k]+(f3[k]+(f4[k]+(f5[k]+f6[k]*wj)*wj)*wj)*wj)*wj)*wj
			wj *= complex(wr, wi)
		}
		w *= w0
	}
	return res
}

// fftBluestein returns the discrete Fourier transform of x using Bluestein's
// algorith, to employ the circular convolution theorem on slices padded to a
// power of 2 length.
func fftBluestein(x []complex128) []complex128 {
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

	y = fftStockham(y, 1)
	for i, ww := range fftStockham(w, 1) {
		y[i] *= ww
	}
	y = fftStockham(y, -1)

	for i := 0; i < n; i++ {
		y[i] *= complex(real(w[i])/float64(m), -imag(w[i])/float64(m))
	}

	return y[:n]
}

// ifftBluestein returns the scaled inverse discrete Fourier transform of x using
// Bluestein's algorithm.
func ifftBluestein(x []complex128) []complex128 {
	n := len(x)
	y := make([]complex128, n)
	for i, xi := range x {
		y[i] = complex(real(xi), -imag(xi))
	}
	y = fftBluestein(y)
	for i, yi := range y {
		y[i] = complex(real(yi)/float64(n), -imag(yi)/float64(n))
	}
	return y
}
