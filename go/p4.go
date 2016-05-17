package main

import (
	"math"
)

type PxPyPzE [4]float64

func (p4 *PxPyPzE) Px() float64 { return p4[0] }
func (p4 *PxPyPzE) Py() float64 { return p4[1] }
func (p4 *PxPyPzE) Pz() float64 { return p4[2] }
func (p4 *PxPyPzE) E() float64  { return p4[3] }

func (p4 *PxPyPzE) M() float64 {
	m2 := p4.M2()
	if m2 < 0.0 {
		return -math.Sqrt(-m2)
	}
	return +math.Sqrt(+m2)
}

func (p4 *PxPyPzE) M2() float64 {
	px := p4.Px()
	py := p4.Py()
	pz := p4.Pz()
	e := p4.E()

	m2 := e*e - (px*px + py*py + pz*pz)
	return m2
}

type byEne []PxPyPzE

func (p byEne) Len() int           { return len(p) }
func (p byEne) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }
func (p byEne) Less(i, j int) bool { return p[i].E() > p[j].E() }
