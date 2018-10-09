// Copyright (c) 2018, Jack Parkinson. All rights reserved.
// Use of this source code is governed by the BSD 3-Clause
// license that can be found in the LICENSE file.

package fft

const maxRadix = 7

// Fft returns the discrete Fourier transform of x or the (normalised) inverse
// transform if inverse is true. The algorithm is most efficient when the length
// of x is a power of 2, 3, 5, 6 or 7 (in that order).
//
// If the length of x is a power of greater than 7 or is prime, then Bluestein's
// algorithm is used to transform the data and retain O(NlogN) performance.
//
// Fft does not check for NaN or Inf values and will produce erroneous results
// if these are present in x.
func Fft(x []complex128, inverse bool) []complex128 {
	n := len(x)
	if n < 2 {
		return x
	}

	s := 1
	if inverse {
		s = -1
	}

	var res []complex128
	switch r := radix(n); r {
	case 2:
		res = stockham(x, s)
	case 3:
		res = radix3(x, s)
	case 5:
		res = radix5(x, s)
	case 6:
		res = radix6(x, s)
	case 7:
		res = radix7(x, s)
	default:
		if s > 0 {
			res = bluestein(x)
		} else {
			res = bluesteini(x)
		}
		goto end
	}

	// Rescale by n for inverse only.
	if s < 0 {
		for i := 0; i < n; i++ {
			res[i] = complex(real(res[i])/float64(n), imag(res[i])/float64(n))
		}
	}

end:
	return res
}

// radix returns the smallest integer r ≤ max_radix such that n = r**p
// for some integer p.
func radix(n int) int {
	if n&(n-1) == 0 {
		return 2
	}

	for i := 3; i <= maxRadix; i++ {
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
