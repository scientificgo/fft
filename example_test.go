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
	"github.com/scientificgo/scigo/fft"
	"math/cmplx"
)

func Example() {
	// Define some data
	x := []complex128{0, complex(1, -2), complex(2, 8), complex(3, 2)}
	fmt.Printf("x = %v\n", x)

	// Get the transform
	ft := fft.Fft(x)
	fmt.Printf("DFT(x) = %v\n", ft)

	// Get an estimate of the PSD
	psd := fft.Periodogram(x)
	fmt.Printf("PSD(x) ≈ %v\n", psd)

	// Estimate the autocorrelation
	// Convert []float64 psd into a []complex128 slice to pass to Ifft
	ac := make([]complex128, len(psd))
	for i, p := range psd {
		ac[i] = complex(p, 0)
	}
	ac = fft.Ifft(ac)
	fmt.Printf("Corr(x, x) ≈ %v\n", ac)

	// Check that IFFT(FFT(X)) returns original X (to within a specified tolerance,
	// owing to rounding point errors)
	xx := fft.Ifft(ft)
	equal := true
	tolerance := 1e-12
	for i := 0; i < len(xx); i++ {
		if cmplx.Abs(xx[i]-x[i]) > tolerance {
			equal = false
		}
	}

	fmt.Printf("Ifft(Fft(x)) == x = %v", equal)

	// Output:
	// x = [(0+0i) (1-2i) (2+8i) (3+2i)]
	// DFT(x) = [(6+8i) (-6-6i) (-2+8i) (2-10i)]
	// PSD(x) ≈ [25 18 17 26]
	// Corr(x, x) ≈ [(21.5+0i) (1.9999999999999998-2i) (-0.5+0i) (2+2i)]
	// Ifft(Fft(x)) == x = true
}
