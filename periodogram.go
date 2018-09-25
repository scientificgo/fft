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
