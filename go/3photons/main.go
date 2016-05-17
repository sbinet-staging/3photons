// 3photons is a simple simulation.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"os"
	"runtime/pprof"
	"sort"
	"strings"
	"time"

	"github.com/go-hep/fmom"
)

var (
	sqrt = math.Sqrt
	pow  = math.Pow
	abs  = math.Abs
)

func main() {
	doprof := flag.Bool("prof", false, "enable profiling")
	flag.Parse()
	if *doprof {
		f, err := os.Create("prof.out")
		if err != nil {
			log.Fatal(err)
		}
		defer f.Close()
		err = pprof.StartCPUProfile(f)
		if err != nil {
			log.Fatal(err)
		}
		defer pprof.StopCPUProfile()
	}

	begTime := time.Now()
	input, err := os.Open("values")
	if err != nil {
		log.Fatal(err)
	}
	defer input.Close()

	var (
		itot = 0
		ntot = 0
		str  = ""

		etot    = 0.0
		cut     = Cut{}
		alpha   = 0.0
		alphaz  = 0.0
		convers = 0.0
		param   = Param{}
		sin2w   = 0.0
		brepem  = 0.0
		betap   = 0.0
		betam   = 0.0
		nbin    = 0
	)

	scan := bufio.NewScanner(input)
	scanf := func(format string, val interface{}, comment *string) {
		if !scan.Scan() {
			return
		}
		err := scan.Err()
		if err != nil {
			log.Fatal(err)
		}
		line := scan.Text()
		_, err = fmt.Sscanf(line, format, val)
		if err != nil {
			log.Fatal(err)
		}

		beg := strings.Index(line, "'") + 1
		end := strings.LastIndex(line, "'")
		*comment = strings.Replace(line[beg:end], "''", "'", -1)
	}

	scanf("%d", &itot, &str)
	fmt.Printf("%d | %s\n", itot, str)

	scanf("%f", &etot, &str)
	fmt.Printf("%v | %s\n", etot, str)

	scanf("%f", &cut.A, &str)
	fmt.Printf("%v | %s\n", cut.A, str)

	scanf("%f", &cut.B, &str)
	fmt.Printf("%v | %s\n", cut.B, str)

	scanf("%f", &cut.Emin, &str)
	scanf("%f", &cut.Sin, &str)
	scanf("%f", &alpha, &str)
	scanf("%f", &alphaz, &str)
	scanf("%f", &convers, &str)
	scanf("%f", &param.MZ, &str)
	scanf("%f", &param.GZ, &str)
	scanf("%f", &sin2w, &str)
	scanf("%f", &brepem, &str)
	scanf("%f", &betap, &str)
	scanf("%f", &betam, &str)
	scanf("%d", &nbin, &str)

	fmt.Printf("itot       %v\n", itot)
	fmt.Printf("etot       %v\n", etot)
	fmt.Printf("cut.A      %v\n", cut.A)
	fmt.Printf("cut.B      %v\n", cut.B)
	fmt.Printf("cut.Emin   %v\n", cut.Emin)
	fmt.Printf("alpha      %v\n", alpha)
	fmt.Printf("alphaz     %v\n", alphaz)
	fmt.Printf("convers    %v\n", convers)
	fmt.Printf("param.MZ   %v\n", param.MZ)
	fmt.Printf("param.GZ   %v\n", param.GZ)
	fmt.Printf("sin2w      %v\n", sin2w)
	fmt.Printf("brepem     %v\n", brepem)
	fmt.Printf("beta+      %v\n", betap)
	fmt.Printf("beta-      %v\n", betam)
	fmt.Printf("nbin       %v\n", nbin)
	fmt.Printf("1/alpha    %v\n", 1/alpha)
	fmt.Printf("1/alphaz   %v\n", 1/alphaz)

	// final state particles, flux
	var (
		mfin  = make([]float64, 100)
		flux  = 1 / (2 * etot * etot)
		norm  = pow(twopi, 4-3*INP) / (float64(itot))
		fact  = 1.0 / 6.0 * convers
		e2    = 4 * math.Pi * alpha
		e2z   = 4 * math.Pi * alphaz
		cos2w = 1 - sin2w
		gzr   = param.GZ / param.MZ
	)

	// couplings
	param.GA = -pow(sqrt(e2), 3)
	param.GBP = -sqrt(e2z/(4*cos2w*sin2w)) / pow(param.MZ, 4)
	param.GBM = param.GBP

	// sum over polarizations factors
	param.POLP = 0 - 2*sin2w
	param.POLM = 1 - 2*sin2w
	param.POLP2 = pow(param.POLP, 2)
	param.POLM2 = pow(param.POLM, 2)
	paa := 2.0
	pab := 1 - 4*sin2w
	pbb := 1 - 4*sin2w + 8*pow(sin2w, 2)

	// homogeneity coefficients
	caa := fact * paa
	cab := fact * pab / pow(param.MZ, 2)
	cbb := fact * pbb / pow(param.MZ, 4)

	// weigth of a 3-photon phase space event
	wtev := pow(math.Pi, 2) / 8 * etot * etot
	dzeta := pow((etot / param.MZ), 2)
	epeak := (dzeta - 1) / gzr
	propa := 1 / (pow(epeak, 2) + 1)

	resfin := ResFin{}
	sigma := 0.0
	variance := 0.0
	weight := 0.0

	var (
		result Result
		spinor Spinor
		scalar Scalar
	)
	evt := NewEvent(etot, mfin)
	for i := 0; i < itot; i++ {
		if i%(itot/10) == 0 {
			fmt.Printf("evt %d\n", i)
		}
		evt.Rambo(INP, etot, mfin, &wtev)
		wtev *= norm

		evt.sort()

		initSpinor(&spinor, &evt)
		initScalar(&scalar, &spinor)
		angle := NewAngle(&evt, &scalar)

		if !angle.Cut(&evt, &cut) {
			initResult(&result, &param, &spinor, etot)
			for k := 0; k < NRES; k++ {
				resfin.Spm2Dif[k] = 0
				for l1 := 0; l1 < 2; l1++ {
					for l2 := 0; l2 < 2; l2++ {
						for l3 := 0; l3 < 2; l3++ {
							resfin.Spm2Dif[k] += result.M2[l1][l2][l3][k]
							//	fmt.Printf("%d %d %d %d | %v\n", k, l1, l2, l3, result.M2[l1][l2][l3][k])
						}
					}
				}
				resfin.Spm2[0][k] += resfin.Spm2Dif[k]
				resfin.Var[0][k] += resfin.Spm2Dif[k] * resfin.Spm2Dif[k]
			}
			proba := caa * resfin.Spm2Dif[0]
			proba += cbb * (betap*betap*resfin.Spm2Dif[1] + betam*betam*resfin.Spm2Dif[2]) / (gzr * gzr) * propa
			proba += cab * 2 * betap * (epeak*resfin.Spm2Dif[3] - resfin.Spm2Dif[4]) / gzr * propa
			weight = proba * wtev * 0.25
			sigma += weight
			variance += weight * weight
			ntot++
		} else {
			weight = 0
		}
	}

	dtot := 1 / float64(itot)
	dtot1 := 1 / (float64(itot) - 1)
	// relative uncertainties
	for k := 0; k < NRES; k++ {
		v := (resfin.Var[0][k] - resfin.Spm2[0][k]*resfin.Spm2[0][k]*dtot) * dtot1
		resfin.Var[0][k] = sqrt(v*dtot) / abs(resfin.Spm2[0][k]*dtot)
	}
	// copy for opposite spins
	resfin.Spm2[1] = resfin.Spm2[0]
	resfin.Var[1] = resfin.Var[0]

	// polarizations
	for k := 1; k <= 2; k++ {
		resfin.Spm2[0][k] *= param.POLM2
		resfin.Spm2[1][k] *= param.POLP2
	}
	for k := 3; k <= 4; k++ {
		resfin.Spm2[0][k] *= param.POLM
		resfin.Spm2[1][k] *= param.POLP
	}

	// physical coefficients and Z0 propagator
	for sp := 0; sp < 2; sp++ {
		for k := 0; k < NRES; k++ {
			resfin.Spm2[sp][k] *= fact * flux * wtev
		}
		resfin.Spm2[sp][0] = resfin.Spm2[sp][0]
		resfin.Spm2[sp][1] = resfin.Spm2[sp][1] / (gzr * gzr) / pow(param.MZ, 4) * propa
		resfin.Spm2[sp][2] = resfin.Spm2[sp][2] / (gzr * gzr) / pow(param.MZ, 4) * propa
		resfin.Spm2[sp][3] = resfin.Spm2[sp][3] / gzr / pow(param.MZ, 2) * propa * epeak
		resfin.Spm2[sp][4] = resfin.Spm2[sp][4] / gzr / pow(param.MZ, 2) * propa
	}

	betamin := sqrt((resfin.Spm2[0][0] + resfin.Spm2[1][0]) / (resfin.Spm2[0][1] + resfin.Spm2[1][1]))
	ssp := 0.5 * (resfin.Spm2[0][1] + resfin.Spm2[1][1]) / sqrt(resfin.Spm2[0][0]+resfin.Spm2[1][0])
	ssm := 0.5 * (resfin.Spm2[0][2] + resfin.Spm2[1][2]) / sqrt(resfin.Spm2[0][0]+resfin.Spm2[1][0])
	incssp := sqrt(pow(resfin.Spm2[0][1]*resfin.Var[0][1], 2)+pow(resfin.Spm2[1][1]*resfin.Var[1][1], 2)) / abs(resfin.Spm2[0][1]+resfin.Spm2[1][1])
	incssp += sqrt(pow(resfin.Spm2[0][0]*resfin.Var[0][0], 2)+pow(resfin.Spm2[1][0]*resfin.Var[1][0], 2)) / abs(resfin.Spm2[0][0]+resfin.Spm2[1][0]) / 2.0
	incssm := sqrt(pow(resfin.Spm2[0][2]*resfin.Var[0][2], 2)+pow(resfin.Spm2[1][2]*resfin.Var[1][2], 2)) / abs(resfin.Spm2[0][2]+resfin.Spm2[1][2])
	incssm += sqrt(pow(resfin.Spm2[0][0]*resfin.Var[0][0], 2)+pow(resfin.Spm2[1][0]*resfin.Var[1][0], 2)) / abs(resfin.Spm2[0][0]+resfin.Spm2[1][0]) / 2.0

	variance = (variance - sigma*sigma/float64(itot)) / float64(itot-1)
	prec := sqrt(variance/float64(itot)) / abs(sigma/float64(itot))
	sigma *= flux

	delta := time.Since(begTime)

	// store numerical results
	out, err := os.Create("res.dat")
	if err != nil {
		log.Fatal(err)
	}
	defer out.Close()

	fprint := func(format string, args ...interface{}) {
		fmt.Fprintf(out, format, args...)
	}
	fprint(" %v\n\n", time.Now())
	fprint(" Number of events:       %16d\n", itot)
	fprint(" ... after cuts:         %16d\n", ntot)
	fprint(" energy in CM:           %g\n", etot)
	fprint(" cut-cos(photon/beam):   %g\n", cut.A)
	fprint(" cut-cos(photon,photon): %g\n", cut.B)
	fprint(" cut-sin(norm,beam):     %g\n", cut.Sin)
	fprint(" cut-energy:             %g\n", cut.Emin)
	fprint(" 1/alpha:                %g\n", 1/alpha)
	fprint(" 1/alpha(Z):             %g\n", 1/alphaz)
	fprint(" conv-factor GeV-2/pb:   %g\n", convers)
	fprint(" M(Z0):                  %g\n", param.MZ)
	fprint(" Width(Z0):              %g\n", param.GZ)
	fprint(" Sin^2(ThetaW):          %g\n", sin2w)
	fprint(" Br(Z->e+e-):            %g\n", brepem)
	fprint(" Beta+:                  %g\n", betap)
	fprint(" Beta-:                  %g\n", betam)
	fprint(" ------------------------------------\n")
	fprint(" x-section (pb):         %g\n", sigma)
	fprint(" std.dev   (pb):         %g\n", sigma*prec)
	fprint(" rel.precision:          %g\n", prec)
	fprint(" ------------------------------------\n")
	fprint(" Beta min:               %g\n", betamin)
	fprint(" Stat. Signf Beta+:      %g\n", ssp)
	fprint(" Error:                  %g\n", ssp*incssp)
	fprint(" Stat. Signf Beta-:      %g\n", ssm)
	fprint(" Error:                  %g\n", ssm*incssm)
	fprint(" Time:                   %v\n", delta)

	fmt.Printf("delta: %v\n", delta)
}

type Cut struct {
	A    float64 // cut on maximum cosine of (beam, photons) angle
	B    float64 // cut on maximum cosine of (photon, photon) angle
	Emin float64 // cut on minimum photon energy
	Sin  float64 // cut on minimum sine of (beam, normal to the photon plane) angle
}

type Param struct {
	MZ    float64 // Z0 mass (GeV/c2)
	GZ    float64 // Z0 width
	GA    float64 // SM EM coupling contribution
	GBP   float64 // beta+ anomalous contribution EW coupling
	GBM   float64 // beta- anomalous contribution EW coupling
	POLP  float64 // EW polarization factors for beta+ anomalous contribution
	POLM  float64 // EW polarization factors for beta- anomalous contribution
	POLP2 float64 // EW polarization factors for beta+ anomalous contribution
	POLM2 float64 // EW polarization factors for beta- anomalous contribution
	IMPR  bool    // predicate controlling dump of result
}

type Event struct {
	P1   fmom.PxPyPzE      // incoming e- 4-momentum
	P2   fmom.PxPyPzE      // incoming e+ 4-momentum
	Pout [100]fmom.PxPyPzE // outgoing particles

	z    [100]float64
	mtot float64
	nm   int
}

func NewEvent(etot float64, masses []float64) Event {
	evt := Event{
		P1: fmom.PxPyPzE{-0.5 * etot, 0, 0, +0.5 * etot},
		P2: fmom.PxPyPzE{+0.5 * etot, 0, 0, +0.5 * etot},
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
		q [INP]fmom.PxPyPzE
		r fmom.PxPyPzE
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
		q[i] = fmom.PxPyPzE{px, py, pz, e}
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

type byEne []fmom.PxPyPzE

func (p byEne) Len() int           { return len(p) }
func (p byEne) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }
func (p byEne) Less(i, j int) bool { return p[i].E() > p[j].E() }

const twopi = 2 * math.Pi

var po2log = math.Log(0.5 * math.Pi)

const (
	NRES = 8 // number of type of results
	INP  = 3 // number of momenta being generated
	RAC8 = 2 * math.Sqrt2
)

type ResFin struct {
	Spm2Dif [NRES]float64    // element of x-section
	Spm2    [2][NRES]float64 // total z-section
	Var     [2][NRES]float64 // variance of the sum
}

type Spinor struct {
	S [5][5]complex128 // mass-less momenta spinor inner products Gram matrix
	T [5][5]complex128 // mass-less momenta conjugate spinor inner products Gram matrix
}

func initSpinor(spinor *Spinor, evt *Event) {
	var (
		p = [5]fmom.PxPyPzE{
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

type Scalar [5][5]float64

func initScalar(s *Scalar, spinor *Spinor) {
	var cx complex128
	// Lorentz scalar product
	for j := 0; j < 4; j++ {
		for k := j + 1; k < 5; k++ {
			cx = spinor.S[j][k]
			s[j][k] = 0.5 * real(cx*cmplx.Conj(cx))
			s[k][j] = s[j][k]
		}
	}
	//return s
}

type Angle struct {
	cos1p1 float64 // cosine between highest energy photon and electron beam
	cos2p1 float64 // cosine between middle  energy photon and electron beam
	cos3p1 float64 // cosine between lowest  energy photon and electron beam
	cos12  float64 // cosine between highest energy photon and middle energy photon
	cos13  float64 // cosine between highest energy photon and lowest energy photon
	cos23  float64 // cosine between middle  energy photon and lowest energy photon
	cosn   float64 // cosine between outgoing photons plane and elextron beam
	cosac  float64 // cosine between outgoing photons plane perpendicular momenta
}

func NewAngle(evt *Event, scalar *Scalar) Angle {
	p := [5]fmom.PxPyPzE{
		evt.P1,
		evt.P2,
		evt.Pout[0],
		evt.Pout[1],
		evt.Pout[2],
	}

	nx := evt.Pout[0].Py()*evt.Pout[1].Pz() - evt.Pout[1].Py()*evt.Pout[0].Pz()
	ny := evt.Pout[0].Pz()*evt.Pout[1].Px() - evt.Pout[1].Pz()*evt.Pout[0].Px()
	nz := evt.Pout[0].Px()*evt.Pout[1].Py() - evt.Pout[1].Px()*evt.Pout[0].Py()
	nn := sqrt(nx*nx + ny*ny + nz*nz)

	angle := Angle{
		cos1p1: 1 - (*scalar)[0][2]/(p[0][3]*p[2][3]),
		cos2p1: 1 - (*scalar)[0][3]/(p[0][3]*p[3][3]),
		cos3p1: 1 - (*scalar)[0][4]/(p[0][3]*p[4][3]),
		cos12:  1 - (*scalar)[2][3]/(evt.Pout[0].E()*evt.Pout[1].E()),
		cos13:  1 - (*scalar)[2][4]/(evt.Pout[0].E()*evt.Pout[2].E()),
		cos23:  1 - (*scalar)[3][4]/(evt.Pout[1].E()*evt.Pout[2].E()),
		cosn:   (nx*p[0][0] + ny*p[0][1] + nz*p[0][2]) / p[0][3] / nn,
		cosac: (p[2][1]*p[3][1] + p[2][2]*p[3][2]) / sqrt(
			(p[2][1]*p[2][1]+p[2][2]*p[2][2])*(p[3][1]*p[3][1]+p[3][2]*p[3][2])),
	}

	return angle
}

func (a *Angle) Cut(evt *Event, cut *Cut) bool {
	var o bool
	for i := 0; i < INP; i++ {
		o = o || (evt.Pout[i].E() < cut.Emin)
	}
	o = (o ||
		(abs(a.cos1p1) > cut.A) ||
		(abs(a.cos2p1) > cut.A) ||
		(abs(a.cos3p1) > cut.A) ||
		(a.cos12 > cut.B) ||
		(a.cos13 > cut.B) ||
		(a.cos23 > cut.B) ||
		(abs(a.cosn) < cut.Sin))

	return o
}

// Result is the array of squared matrix elements contribution with detail oh helicities
type Result struct {
	M2 [2][2][2][NRES]float64
}

func initResult(res *Result, param *Param, spinor *Spinor, etot float64) {
	var (
		//res = new(Result)
		a  complex128
		bp complex128
		bm complex128
	)

	// minus signs comes from outgoing photons
	for ll1 := 0; ll1 < 2; ll1++ {
		l1 := -(2*ll1 - 1)
		for ll2 := 0; ll2 < 2; ll2++ {
			l2 := -(2*ll2 - 1)
			for ll3 := 0; ll3 < 2; ll3++ {
				l3 := -(2*ll3 - 1)
				// helicity amplitudes
				switch {
				case l1 == +1 && l2 == +1 && l3 == +1:
					a = 0
					bp = 0
					bm = spinor.BPPP(2, 3, 4)
				//	fmt.Printf("+++ %v %v %v\n", a, bp, bm)

				case l1 == -1 && l2 == -1 && l3 == -1:
					a = 0
					bp = 0
					bm = spinor.BMMM(2, 3, 4)
				//	fmt.Printf("--- %v %v %v\n", a, bp, bm)

				case l1 == +1 && l2 == +1 && l3 == -1:
					a = spinor.APPM(2, 3, 4)
					bp = spinor.BPPM(2, 3, 4)
					bm = 0.0
				//	fmt.Printf("++- %v %v %v\n", a, bp, bm)

				case l1 == +1 && l2 == -1 && l3 == +1:
					a = spinor.APPM(4, 2, 3)
					bp = spinor.BPPM(4, 2, 3)
					bm = 0.0
				//	fmt.Printf("+-+ %v %v %v\n", a, bp, bm)

				case l1 == -1 && l2 == +1 && l3 == +1:
					a = spinor.APPM(3, 4, 2)
					bp = spinor.BPPM(3, 4, 2)
					bm = 0.0
				//	fmt.Printf("-++ %v %v %v\n", a, bp, bm)

				case l1 == +1 && l2 == -1 && l3 == -1:
					a = spinor.APMM(2, 3, 4)
					bp = spinor.BPMM(2, 3, 4)
					bm = 0.0
				//	fmt.Printf("+-- %v %v %v\n", a, bp, bm)

				case l1 == -1 && l2 == -1 && l3 == +1:
					a = spinor.APMM(4, 2, 3)
					bp = spinor.BPMM(4, 2, 3)
					bm = 0.0
					//	fmt.Printf("--+ %v %v %v\n", a, bp, bm)

				case l1 == -1 && l2 == +1 && l3 == -1:
					a = spinor.APMM(3, 2, 4)
					bp = spinor.BPMM(3, 2, 4)
					bm = 0.0
					//	fmt.Printf("-+- %v %v %v\n", a, bp, bm)
				}

				// couplings
				a *= complex(param.GA, 0)
				bp *= complex(param.GBP, 0)
				bm *= complex(param.GBM, 0)

				// fmt.Printf("== %v %v %v\n", a, bp, bm)

				aabs := cmplx.Abs(a)
				a2 := aabs * aabs
				bpabs := cmplx.Abs(bp)
				bp2 := bpabs * bpabs
				bmabs := cmplx.Abs(bm)
				bm2 := bmabs * bmabs

				abp := a * cmplx.Conj(bp)

				//squared matrix elements
				res.M2[ll1][ll2][ll3][0] = a2
				res.M2[ll1][ll2][ll3][1] = bp2
				res.M2[ll1][ll2][ll3][2] = bm2
				res.M2[ll1][ll2][ll3][3] = 2.0 * real(abp)
				res.M2[ll1][ll2][ll3][4] = 2.0 * imag(abp)
				// res.M2[ll1][ll2][ll3][5] = 0.0
				// res.M2[ll1][ll2][ll3][6] = 0.0
				// res.M2[ll1][ll2][ll3][7] = 0.0

				// fmt.Printf("l1-l2-l3-0 %v\n", res.M2[ll1][ll2][ll3][0])
			}
		}
	}
	//return res
}

func b2i(b bool) int {
	if b {
		return 1
	}
	return 0
}

const modulo = 1000000000

var (
	ncall = 0
	mcall = 55
	ia    [56]int64
)

// rn returns a random number between 0 and 1.
// Rand implements a random number function taken from Knuth's RANF
// (Semi numerical algorithms).
//
// Method is X(N)=MOD(X(N-55)-X(N-24), 1/FMODUL)
func rn() float64 {
	const fmodul = 1e-9
	if ncall == 0 {
		in55(&ia, int64(234612947))
		ncall = 1
	}
	if mcall == 0 {
		irn55(&ia)
		mcall = 55
	}
	mcall -= 1
	return float64(ia[mcall+1]) * fmodul
}

func in55(ia *[56]int64, ix int64) {
	(*ia)[55] = ix
	j := ix
	k := int64(1)
	for i := 1; i <= 54; i++ {
		ii := (21 * i) % 55
		(*ia)[ii] = k
		k = j - k
		if k < 0 {
			k += modulo
		}
		j = (*ia)[ii]
	}
	for i := 0; i < 10; i++ {
		irn55(ia)
	}
}

func irn55(ia *[56]int64) {
	var j int64
	for i := 1; i <= 24; i++ {
		j = (*ia)[i] - (*ia)[i+31]
		if j < 0 {
			j += modulo
		}
		(*ia)[i] = j
	}
	for i := 25; i <= 55; i++ {
		j = (*ia)[i] - (*ia)[i-24]
		if j < 0 {
			j += modulo
		}
		(*ia)[i] = j
	}
}
