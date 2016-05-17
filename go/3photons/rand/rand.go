package rand

const modulo = 1000000000

var (
	ncall = 0
	mcall = 55
	ia    [56]int64
)

// Rand returns a random number between 0 and 1.
// Rand implements a random number function taken from Knuth's RANF
// (Semi numerical algorithms).
//
// Method is X(N)=MOD(X(N-55)-X(N-24), 1/FMODUL)
func Rand() float64 {
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
