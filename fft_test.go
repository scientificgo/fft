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

const tol = 8 // significant figures

var cases = [][]int{
	{2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384},
	{3, 9, 27, 81, 243, 729, 2187, 6561, 19683},
	{5, 25, 125, 625, 3125, 15625},
	{6, 36, 216, 1296, 7776},
	{7, 49, 343, 2401, 16807},
	{1, 17, 31, 61, 127, 257, 509, 1021, 2053, 4093},
}

var inputs [][]complex128
var labels []string

func init() {
	max := 1
	ncases := 0
	for _, c := range cases {
		for _, n := range c {
			if n > max {
				max = n
			}
			ncases++
		}
	}
	x := make([]complex128, max)
	for i := 0; i < max; i++ {
		x[i] = complex(rand.NormFloat64(), rand.NormFloat64())
	}

	inputs = make([][]complex128, ncases)
	labels = make([]string, ncases)
	ncase := 0
	for _, c := range cases {
		for _, n := range c {
			labels[ncase] = fmt.Sprintf("%v-%v", c[0], n)
			inputs[ncase] = x[:n]
			ncase++
		}
	}
}

//
// Tests
//

func dftDirect(x []complex128, inverse bool) []complex128 {
	n := len(x)
	s := 1
	if inverse {
		s = -1
	}

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

func test(t *testing.T, inverse bool) {
	var f func(*testing.T, float64, []string, [](func([]complex128) []complex128), [][]complex128)
	testutils.GenerateTest(&f)

	dft := func(x []complex128) []complex128 { return dftDirect(x, inverse) }
	fft := func(x []complex128) []complex128 { return Fft(x, inverse) }

	f(t, tol, labels, [](func([]complex128) []complex128){fft, dft}, inputs)
}

func TestFft(t *testing.T)  { test(t, false) }
func TestIfft(t *testing.T) { test(t, true) }

//
// Benchmarks
//

func benchmark(b *testing.B, f func([]complex128, bool) []complex128, inverse bool) {
	for ncase := 0; ncase < len(labels); ncase++ {
		b.Run(labels[ncase], func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_ = f(inputs[ncase], inverse)
			}
		})
	}
}

func BenchmarkFft(b *testing.B)   { benchmark(b, Fft, false) }
func BenchmarkIffti(b *testing.B) { benchmark(b, Fft, true) }
