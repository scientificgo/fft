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

const acc = 8 // significant figures

var nslices = [][]int{
	{2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384},
	{3, 9, 27, 81, 243, 729, 2187, 6561, 19683},
	{5, 25, 125, 625, 3125, 15625},
	{6, 36, 216, 1296, 7776},
	{7, 49, 343, 2401, 16807},
	{1, 17, 31, 61, 127, 257, 509, 1021, 2053, 4093},
}

var cases = []struct {
	Label string
	Input []complex128
}{}

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

//
// Setup
//

func init() {
	max := 1
	for _, nslice := range nslices {
		for _, n := range nslice {
			if n > max {
				max = n
			}
		}
	}
	x := make([]complex128, max)
	for i := 0; i < max; i++ {
		x[i] = complex(rand.NormFloat64(), rand.NormFloat64())
	}

	for _, nslice := range nslices {
		for _, n := range nslice {
			cases = append(cases, struct {
				Label string
				Input []complex128
			}{fmt.Sprintf("%v-%v", nslice[0], n), x[:n]})
		}
	}
}

//
// Test
//

func test(t *testing.T, inverse bool) {
	dft := func(x []complex128) []complex128 { return dftDirect(x, inverse) }
	fft := func(x []complex128) []complex128 { return Fft(x, inverse) }
	testutils.Test(t, acc, cases, fft, dft)
}

func benchmark(b *testing.B, f func([]complex128, bool) []complex128, inverse bool) {
	for ncase := 0; ncase < len(cases); ncase++ {
		b.Run(cases[ncase].Label, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_ = f(cases[ncase].Input, inverse)
			}
		})
	}
}

func TestFft(t *testing.T)        { test(t, false) }
func TestIfft(t *testing.T)       { test(t, true) }
func BenchmarkFft(b *testing.B)   { benchmark(b, Fft, false) }
func BenchmarkIffti(b *testing.B) { benchmark(b, Fft, true) }
