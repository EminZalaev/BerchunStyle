package main

import (
	"fmt"
	"os"
)

func main() {
	N := 48
	Tc := 368
	Ts := 8
	lambda := 1 / float64(Tc)
	mu := 1.0 / float64(Ts)
	//K := make([]float64, N)

	MfromNwFile, _ := os.Create("MfromNw.txt")
	defer MfromNwFile.Close()

	MQFile, _ := os.Create("MQ.txt")
	defer MQFile.Close()

	PQFile, _ := os.Create("PQ.txt")
	defer PQFile.Close()

	MnFile, _ := os.Create("Mn.txt")
	defer MnFile.Close()

	KznFile, _ := os.Create("Kzn.txt")
	defer KznFile.Close()

	for n := 1; n < 10; n++ {
		Kcalc := Kfunc(n, N, lambda, mu)
		nString := fmt.Sprintf("%d", n)

		fmt.Println("K: ", Kcalc)

		Po := PoFunc(Kcalc, N)
		fmt.Println("Po: ", Po)

		MNw := MfromNw(Kcalc, n, N, Po)
		fmt.Println("M(Nw): ", MNw)

		MNwString := fmt.Sprintf("%f", MNw)
		MfromNwFile.WriteString(nString + "	" + MNwString + "\n")

		MQ := MfromQ(Kcalc, n, N, Po)
		fmt.Println("M(Q): ", MQ)

		MQString := fmt.Sprintf("%f", MQ)
		MQFile.WriteString(nString + "	" + MQString + "\n")

		PQ := PfromQ(Kcalc, n, N, Po)
		fmt.Println("P(Q): ", PQ)

		PQString := fmt.Sprintf("%f", PQ)
		PQFile.WriteString(nString + "	" + PQString + "\n")

		MN := MfromN(Kcalc, n, N, Po)
		fmt.Println("M(N): ", MN)

		MNString := fmt.Sprintf("%f", MN)
		MnFile.WriteString(nString + "	" + MNString + "\n")

		Kzn := MN / float64(n)
		fmt.Println("Kzn: ", Kzn)

		KznString := fmt.Sprintf("%f", Kzn)
		KznFile.WriteString(nString + "	" + KznString + "\n")
	}
}

func Kfunc(n int, N int, lambda float64, mu float64) []float64 {
	K := make([]float64, N)
	K[0] = (float64(N) * lambda) / mu

	for i := 1; i < N; i++ {
		if i < n {
			K[i] = K[i-1] * ((float64(N-i) * lambda) / (float64(i+1) * mu))
		} else {
			K[i] = K[i-1] * ((float64(N-i) * lambda) / (float64(n) * mu))
		}
	}
	return K
}

func PoFunc(K []float64, N int) float64 {
	var Po float64
	sum := 1.0
	for i := 0; i < N; i++ {
		sum += K[i]
	}
	Po = float64(1) / sum

	return Po
}

func MfromNw(K []float64, n int, N int, Po float64) float64 {
	var M_N float64

	for i := 0; i < N; i++ {
		M_N += float64(i+1) * K[i]
	}
	M_N = Po * M_N

	return M_N
}

func MfromQ(K []float64, n int, N int, Po float64) float64 {
	var M_N float64
	j := 1

	for i := n; i < N; i++ {
		M_N += float64(j) * K[i]
		j++
	}
	M_N = Po * M_N

	return M_N
}

func PfromQ(K []float64, n int, N int, Po float64) float64 {
	var M_N float64

	for i := n; i < N; i++ {
		M_N += K[i]
	}
	M_N = Po * M_N

	return M_N
}

func MfromN(K []float64, n int, N int, Po float64) float64 {
	var M_N float64

	for i := 0; i < N; i++ {
		if i < n {
			M_N += float64(i+1) * K[i]
		} else {
			M_N += float64(n) * K[i]
		}
	}
	M_N = Po * M_N

	return M_N
}
