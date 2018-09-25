// Copyright (c) 2018, Jack Parkinson. All rights reserved.
// Use of this source code is governed by the BSD 3-Clause
// license that can be found in the LICENSE file.

package fft

// Periodogram returns the periodogram of x. It is defined as
// the modulus squared of the discrete Fourier transform of x.
func Periodogram(x []complex128) []float64 {
	n := len(x)
	p := make([]float64, n)
	for i, f := range Fft(x) {
		re, im := real(f), imag(f)
		p[i] = (re*re + im*im) / float64(n)
	}
	return p
}
