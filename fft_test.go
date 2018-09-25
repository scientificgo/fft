/*
   SciGo is a scientific library for the Go language.
   Copyright (C) 2018, Jack Parkinson

   This program is free software: you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package fft_test

import (
	"fmt"
	. "github.com/scientificgo/fft"
	"github.com/scientificgo/utils"
	"math"
	"math/rand"
	"testing"
)

const tol = 6

// dftDirect returns the DFT (flag=1) or inverse DFT (flag = -1) of x using the direct O(n**2) method.
// It is used to validate the values produced by the more efficient FFT routines exposed in the package fft API.
func dftDirect(x []complex128, flag int) []complex128 {
	if flag != -1 && flag != 1 {
		return nil
	}

	n := len(x)
	res := make([]complex128, n)
	for i := 0; i < n; i++ {
		for k := 0; k < n; k++ {
			wi, wr := math.Sincos(-2 * math.Pi * float64(flag*k*i) / float64(n))
			res[k] += x[i] * complex(wr, wi)
		}
	}
	return res
}

func periodogramDirect(x []complex128) []float64 {
	n := len(x)
	p := make([]float64, n)
	for i, f := range dftDirect(x, 1) {
		re, im := real(f), imag(f)
		p[i] = (re*re + im*im) / float64(n)
	}
	return p
}

var cases = [][]int{
	{2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096},
	{3, 9, 27, 81, 243, 729, 2187, 6561},
	{5, 25, 125, 625, 3125, 15625},
	{6, 36, 216, 1296, 7776},
	{7, 49, 343, 2401},
	{1, 17, 31, 61, 127, 257, 509, 1021, 2053, 4093},
}

func TestFft(t *testing.T) {
	var x []complex128
	for _, ns := range cases {
		for _, n := range ns {
			for len(x) < n {
				x = append(x, complex(100*rand.NormFloat64(), 100*rand.NormFloat64()))
			}
			xn := x[:n]
			t.Run(fmt.Sprintf("n=%v", n), func(t *testing.T) {
				direct := dftDirect(xn, 1)
				fft := Fft(xn)
				if !utils.EqualComplex128s(fft, direct, tol) {
					t.Errorf("Fft(%v) = %v, want %v", xn, fft, direct)
				}
				xx := Ifft(fft)
				if !utils.EqualComplex128s(xx, xn, tol) {
					t.Errorf("Ifft(%v) = %v, want %v", fft, xx, xn)
				}
				pdirect := periodogramDirect(xn)
				p := Periodogram(xn)
				if !utils.EqualFloat64s(p, pdirect, tol) {
					t.Errorf("PSD(%v) = %v, want %v", xn, p, pdirect)
				}
			})
		}
	}
}

//
// Benchmarks
//

var GlobalI int

func benchmarkFft(f func([]complex128) []complex128, b *testing.B) {
	var x []complex128
	for _, ns := range cases {
		var y []complex128
		for _, n := range ns {
			for len(x) < n {
				x = append(x, complex(100*rand.NormFloat64(), 100*rand.NormFloat64()))
			}
			b.Run(fmt.Sprintf("%v-%v", ns[0], n), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					y = f(x[:n])
				}
				GlobalI = len(y)
			})
		}
	}
}

func BenchmarkFft(b *testing.B)  { benchmarkFft(Fft, b) }
func BenchmarkIfft(b *testing.B) { benchmarkFft(Ifft, b) }
