package main

import (
	"fmt"
	"log"
	"math"
	"sort"
)

type Event struct {
	P1   PxPyPzE      // incoming e- 4-momentum
	P2   PxPyPzE      // incoming e+ 4-momentum
	Pout [100]PxPyPzE // outgoing particles

	z    [100]float64
	mtot float64
	nm   int
}

func NewEvent(etot float64, masses []float64) Event {
	evt := Event{
		P1: PxPyPzE{-0.5 * etot, 0, 0, +0.5 * etot},
		P2: PxPyPzE{+0.5 * etot, 0, 0, +0.5 * etot},
	}
	evt.init(etot, masses)
	return evt
}

func (e *Event) display() {
	for i, p := range e.Pout[:4] {
		fmt.Printf("%d: %v\n", i, p)
	}
}

func (e *Event) sort() {
	sort.Sort(byEne(e.Pout[:4]))
}

func (e *Event) init(etot float64, masses []float64) {
	fmt.Printf("initializing...\n")
	e.z[1] = po2log
	for i := range e.z {
		if i < 2 {
			continue
		}
		e.z[i] = e.z[i-1] + po2log - 2*math.Log(float64(i-1))
	}
	for i, v := range e.z {
		if i < 2 {
			continue
		}
		e.z[i] = v - math.Log(float64(i))
	}

	if len(masses) < 1 || len(masses) > 100 {
		log.Fatal("invalid number of particles: %d\n", len(masses))
	}

	// check whether total energy is sufficient
	e.mtot = 0
	e.nm = 0
	for _, m := range masses {
		if m != 0 {
			e.nm++
		}
		e.mtot += abs(m)
	}
	if e.mtot > etot {
		log.Fatal("mtot > etot (%v>%v)\n", e.mtot, etot)
	}
	fmt.Printf("initializing... [done]\n")
}

func (e *Event) Rambo(n int, etot float64, ms []float64, weight *float64) {
	var (
		q [INP]PxPyPzE
		r PxPyPzE
		b [3]float64
	)

	// generate massless momenta in infinite phase space
	for i := 0; i < n; i++ {
		c := 2*rn() - 1
		s := sqrt(1 - c*c)
		f := twopi * rn()
		e := -math.Log(rn() * rn())
		pz := e * c
		py := e * s * math.Cos(f)
		px := e * s * math.Sin(f)
		q[i] = PxPyPzE{px, py, pz, e}
	}

	// compute the parameters of the conformal transformation
	for i := 0; i < n; i++ {
		v := q[i]
		r[0] += v.Px()
		r[1] += v.Py()
		r[2] += v.Pz()
		r[3] += v.E()
	}
	rmass := r.M()
	invrmass := 1 / rmass
	b[0] = -r[0] * invrmass
	b[1] = -r[1] * invrmass
	b[2] = -r[2] * invrmass
	g := r[3] * invrmass
	a := 1 / (1 + g)
	x := etot * invrmass

	// transorfms the q's conformally into the p's
	for i, qq := range q {
		pout := e.Pout[i]
		bq := b[0]*qq.Px() + b[1]*qq.Py() + b[2]*qq.Pz()
		for k := 0; k < 3; k++ {
			pout[k] = x * (qq[k] + b[k]*(qq.E()+a*bq))
		}
		pout[3] = x * (g*qq.E() + bq)
		e.Pout[i] = pout
	}

	// computes weights
	wt := po2log
	if n != 2 {
		wt = (2*float64(n)-4)*math.Log(etot) + e.z[n-1]
	}
	*weight = math.Exp(wt)
}
