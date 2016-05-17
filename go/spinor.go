package main

import (
	"math/cmplx"
)

type Spinor struct {
	S [5][5]complex128 // mass-less momenta spinor inner products Gram matrix
	T [5][5]complex128 // mass-less momenta conjugate spinor inner products Gram matrix
}

func initSpinor(spinor *Spinor, evt *Event) {
	var (
		p = [5]PxPyPzE{
			evt.P1,
			evt.P2,
			evt.Pout[0],
			evt.Pout[1],
			evt.Pout[2],
		}
		xx [5]float64
		fx [5]complex128
	)

	for k := 0; k < 5; k++ {
		tx := sqrt(p[k].E() + p[k].Pz())
		xx[k] = tx
		fx[k] = complex(p[k].Px(), p[k].Py()) / complex(tx, 0)
	}

	for j := 0; j < 4; j++ {
		for k := j + 1; k < 5; k++ {
			// spinor product Mangano, Sparke
			cx := fx[j]*complex(xx[k], 0) - fx[k]*complex(xx[j], 0)
			spinor.S[j][k] = +cx
			spinor.S[k][j] = -cx
			spinor.T[k][j] = +cmplx.Conj(cx)
			spinor.T[j][k] = -cmplx.Conj(cx)
		}
	}
	//return spinor
}

func (s *Spinor) APPM(k1, k2, k3 int) complex128 {
	return -RAC8 * s.S[0][1] * s.S[0][k3] * s.S[0][k3] / (s.S[0][k1] * s.S[0][k2] * s.S[1][k1] * s.S[1][k2])
}

func (s *Spinor) APMM(k1, k2, k3 int) complex128 {
	return -RAC8 * s.T[0][1] * s.T[1][k1] * s.T[1][k1] / (s.T[1][k2] * s.T[1][k3] * s.T[0][k2] * s.T[0][k3])
}

func (s *Spinor) BPPM(k1, k2, k3 int) complex128 {
	v := s.T[k1][k2] * s.S[k3][0]
	v2 := v * v
	return -RAC8 * s.T[0][1] * v2
}

func (s *Spinor) BPMM(k1, k2, k3 int) complex128 {
	v := s.T[k1][1] * s.S[k2][k3]
	v2 := v * v
	return -RAC8 * s.S[0][1] * v2
}

func (s *Spinor) BPPP(k1, k2, k3 int) complex128 {
	v123 := s.T[k1][k2] * s.T[k3][1]
	v132 := s.T[k1][k3] * s.T[k2][1]
	v231 := s.T[k2][k3] * s.T[k1][1]
	sum := v123*v123 + v132*v132 + v231*v231
	return -RAC8 * s.S[0][1] * sum
}

func (s *Spinor) BMMM(k1, k2, k3 int) complex128 {
	v123 := s.S[k1][0] * s.S[k2][k3]
	v213 := s.S[k2][0] * s.S[k1][k3]
	v312 := s.S[k3][0] * s.S[k1][k2]
	sum := v123*v123 + v213*v213 + v312*v312
	return -RAC8 * s.T[0][1] * sum
}
