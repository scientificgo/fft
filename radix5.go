package fft

import "math"

const (
	w51 = 0.30901699437494742410229341718281905886015458990288143106 - 0.95105651629515357211643933337938214340569863412575022244i  // exp(-2πi*1/5)
	w52 = -0.80901699437494742410229341718281905886015458990288143107 - 0.58778525229247312916870595463907276859765243764314599107i // exp(-2πi*2/5)
	w53 = -0.80901699437494742410229341718281905886015458990288143107 + 0.58778525229247312916870595463907276859765243764314599107i // exp(-2πi*3/5)
	w54 = 0.30901699437494742410229341718281905886015458990288143106 + 0.95105651629515357211643933337938214340569863412575022244i  // exp(-2πi*4/5)
)

func radix5(x []complex128, s int) []complex128 {
	const r = 5
	n := len(x)

	// Copy input into new storage.
	y := make([]complex128, n)
	copy(y, x)

	// Reorder input using base-r digit reversal permutation.
	for i, j := 0, 0; i < n-1; i++ {
		if i < j {
			y[i], y[j] = y[j], y[i]
		}
		k := (r - 1) * n / r
		for k <= j {
			j -= k
			k /= r
		}
		j += k / (r - 1)
	}

	for m := r; m <= n; m *= r {
		// Calculate twiddle factor.
		wim, wre := math.Sincos(-2 * math.Pi * float64(s) / float64(m))
		w := complex(wre, wim)

		mr := m / r
		for i, wi := 0, complex(1, 0); i < mr; i++ {
			for j := 0; j < n; j += m {
				// Retrieve subset of points.
				t0, t1, t2, t3, t4 := y[i+j], y[i+j+mr], y[i+j+2*mr], y[i+j+3*mr], y[i+j+4*mr]

				// Apply twiddle factors w**(i+k) for 1 ≤ k < r.
				t1 *= wi
				t2 *= wi * wi
				t3 *= wi * wi * wi
				t4 *= wi * wi * wi * wi

				// Transform points using r-point DFT.
				y[i+j] += t1 + t2 + t3 + t4
				if s > 0 {
					y[i+j+mr], y[i+j+2*mr], y[i+j+3*mr], y[i+j+4*mr] =
						t0+t1*w51+t2*w52+t3*w53+t4*w54,
						t0+t1*w52+t2*w54+t3*w51+t4*w53,
						t0+t1*w53+t2*w51+t3*w54+t4*w52,
						t0+t1*w54+t2*w53+t3*w52+t4*w51
				} else {
					// 1/w51 = w54, etc.
					y[i+j+mr], y[i+j+2*mr], y[i+j+3*mr], y[i+j+4*mr] =
						t0+t1*w54+t2*w53+t3*w52+t4*w51,
						t0+t1*w53+t2*w51+t3*w54+t4*w52,
						t0+t1*w52+t2*w54+t3*w51+t4*w53,
						t0+t1*w51+t2*w52+t3*w53+t4*w54
				}
			}
			wi *= w
		}
	}
	return y
}
