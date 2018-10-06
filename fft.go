// Copyright (c) 2018, Jack Parkinson. All rights reserved.
// Use of this source code is governed by the BSD 3-Clause
// license that can be found in the LICENSE file.

package fft

const max_radix = 7

// Fft returns a discrete Fourier transform of x.
// It does not check x for NaN or Inf values; these checks should be done separately.
func Fft(x []complex128) []complex128 {
	n := len(x)
	if n < 2 {
		return x
	}
	var res []complex128
	switch r := radix(n); r {
	case 2:
		res = stockham(x, 1)
	case 3:
		res = radix3(x, 1)
	case 5:
		res = radix5(x, 1)
	case 6:
		res = radix6(x, 1)
	case 7:
		res = radix7(x, 1)
	default:
		res = bluesteinFwd(x)
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
	switch r := radix(n); r {
	case 2:
		res = stockham(x, -1)
	case 3:
		res = radix3(x, -1)
	case 5:
		res = radix5(x, -1)
	case 6:
		res = radix6(x, -1)
	case 7:
		res = radix7(x, -1)
	default:
		res = bluesteinBwd(x)
		goto end
	}

	rescale(res, float64(n))

end:
	return res
}

func rescale(x []complex128, scale float64) {
	for i := 0; i < len(x); i++ {
		x[i] *= complex(1/scale, 0)
	}
}

// radix returns the smallest integer r ≤ max_radix such that n = r**p
// for some integer p.
func radix(n int) int {
	if n&(n-1) == 0 {
		return 2
	}

	for i := 3; i <= max_radix; i++ {
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
