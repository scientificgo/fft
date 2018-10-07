// Copyright (c) 2018, Jack Parkinson. All rights reserved.
// Use of this source code is governed by the BSD 3-Clause
// license that can be found in the LICENSE file.

package fft_test

import (
	"fmt"
	"math"
	"math/rand"
	. "scientificgo.org/fft"
	"scientificgo.org/testutils"
	"testing"
)

//
// Setup
//

const TOL = 6 // significant figures

var CASES = [][]int{
	{2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384},
	{3, 9, 27, 81, 243, 729, 2187, 6561, 19683},
	{5, 25, 125, 625, 3125, 15625},
	{6, 36, 216, 1296, 7776},
	{7, 49, 343, 2401, 16807},
	{1, 17, 31, 61, 127, 257, 509, 1021, 2053, 4093},
}

var X []complex128

func init() {
	max := 1
	for _, c := range CASES {
		for _, n := range c {
			if n > max {
				max = n
			}
		}
	}
	X = make([]complex128, max)
	for i := 0; i < max; i++ {
		X[i] = complex(rand.NormFloat64(), rand.NormFloat64())
	}
}

//
// Tests
//

// dft_direct returns the DFT (s=1) or inverse DFT (s=-1) of x using the direct O(n**2) method.
// It is used to validate the values produced by the more efficient FFT routines exposed in the package fft API.
func dft_direct(x []complex128, s int) []complex128 {
	if s*s != 1 {
		return nil
	}
	n := len(x)
	res := make([]complex128, n)
	for i := 0; i < n; i++ {
		for k := 0; k < n; k++ {
			wim, wre := math.Sincos(-2 * math.Pi * float64(s*k*i) / float64(n))
			res[k] += x[i] * complex(wre, wim)
		}
	}
	if s < 0 {
		for i := 0; i < n; i++ {
			res[i] /= complex(float64(n), 0)
		}
	}
	return res
}

func test(t *testing.T, f func([]complex128) []complex128, s int) {
	for _, ns := range CASES {
		for _, n := range ns {
			t.Run(fmt.Sprintf("%v-%v", ns[0], n), func(t *testing.T) {
				direct := dft_direct(X[:n], s)
				res := f(X[:n])
				if !testutils.EqualComplex128s(res, direct, TOL) {
					fname := "Fft"
					if s < 0 {
						fname = "Ifft"
					}
					t.Errorf("%v([%v,...,%v]) = [%v,...,%v], want [%v,...,%v]", fname, X[0], X[n-1], res[0], res[n-1], direct[0], direct[n-1])
				}
			})
		}
	}
}

func TestFft(t *testing.T)  { test(t, Fft, 1) }
func TestIfft(t *testing.T) { test(t, Ifft, -1) }

//
// Benchmarks
//

var GlobalI int

func benchmark(b *testing.B, f func([]complex128) []complex128) {
	for _, ns := range CASES {
		var y []complex128
		for _, n := range ns {
			b.Run(fmt.Sprintf("%v-%v", ns[0], n), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					y = f(X[:n])
				}
				GlobalI = len(y)
			})
		}
	}
}

func BenchmarkFft(b *testing.B)  { benchmark(b, Fft) }
func BenchmarkIfft(b *testing.B) { benchmark(b, Ifft) }
