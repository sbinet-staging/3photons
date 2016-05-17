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
	"strings"
	"time"
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
	p := [5]PxPyPzE{
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

				case l1 == -1 && l2 == -1 && l3 == -1:
					a = 0
					bp = 0
					bm = spinor.BMMM(2, 3, 4)

				case l1 == +1 && l2 == +1 && l3 == -1:
					a = spinor.APPM(2, 3, 4)
					bp = spinor.BPPM(2, 3, 4)
					bm = 0.0

				case l1 == +1 && l2 == -1 && l3 == +1:
					a = spinor.APPM(4, 2, 3)
					bp = spinor.BPPM(4, 2, 3)
					bm = 0.0

				case l1 == -1 && l2 == +1 && l3 == +1:
					a = spinor.APPM(3, 4, 2)
					bp = spinor.BPPM(3, 4, 2)
					bm = 0.0

				case l1 == +1 && l2 == -1 && l3 == -1:
					a = spinor.APMM(2, 3, 4)
					bp = spinor.BPMM(2, 3, 4)
					bm = 0.0

				case l1 == -1 && l2 == -1 && l3 == +1:
					a = spinor.APMM(4, 2, 3)
					bp = spinor.BPMM(4, 2, 3)
					bm = 0.0

				case l1 == -1 && l2 == +1 && l3 == -1:
					a = spinor.APMM(3, 2, 4)
					bp = spinor.BPMM(3, 2, 4)
					bm = 0.0
				}

				// couplings
				a *= complex(param.GA, 0)
				bp *= complex(param.GBP, 0)
				bm *= complex(param.GBM, 0)

				aabs := cmplx.Abs(a)
				a2 := aabs * aabs
				bpabs := cmplx.Abs(bp)
				bp2 := bpabs * bpabs
				bmabs := cmplx.Abs(bm)
				bm2 := bmabs * bmabs

				abp := a * cmplx.Conj(bp)

				//squared matrix elements
				m2 := res.M2[ll1][ll2][ll3][:]
				m2[0] = a2
				m2[1] = bp2
				m2[2] = bm2
				m2[3] = 2.0 * real(abp)
				m2[4] = 2.0 * imag(abp)
				m2[5] = 0.0
				m2[6] = 0.0
				m2[7] = 0.0
			}
		}
	}
}
