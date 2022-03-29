package main

import (
	"fmt"
	"math"
	"os"
)

func main() {
	Tc := 49
	Ts := 213
	Tw := 587
	//n := 4

	//lambda := 3.0
	//mu := 1.0

	lambda := 1 / float64(Tc)
	mu := 1 / float64(Ts)
	v := 1 / float64(Tw)
	//part 1.1
	//fmt.Println("\nlambda =", lambda, "\nmu = ", mu)
	//
	//PotkazFile, _ := os.Create("Potkaz.txt")
	//defer PotkazFile.Close()
	//
	//MfromNFile, _ := os.Create("MfromN.txt")
	//defer MfromNFile.Close()
	//
	//KzFile, _ := os.Create("Kz.txt")
	//defer KzFile.Close()
	//
	//for n := 1; n <= 10; n++ {
	//	K := make([]float64, n)
	//
	//	Kcomp := Kfunc(lambda, mu, n, K)
	//	fmt.Println("\nK1: ", Kcomp)
	//
	//	Po := PoCalc(K, n)
	//	fmt.Println("Po: ", Po)
	//
	//	Potkaz := PoCalc(K, n) * K[n-1]
	//	fmt.Println("Pотказ: ", Potkaz)
	//
	//	M_N := MfromN(K, Po, n)
	//	fmt.Println("M(N): ", M_N)
	//
	//	Kz := M_N / float64(n)
	//	fmt.Println("Kz: ", Kz)
	//
	//	PotkazString := fmt.Sprintf("%f", Potkaz)
	//	nString := fmt.Sprintf("%d", n)
	//	PotkazFile.WriteString(nString + " " + PotkazString + "\n")
	//
	//	MfromNString := fmt.Sprintf("%f", M_N)
	//	MfromNFile.WriteString(nString + " " + MfromNString + "\n")
	//
	//	KzString := fmt.Sprintf("%f", Kz)
	//	KzFile.WriteString(nString + " " + KzString + "\n")
	//}
	//
	////part 1.2
	//fmt.Println("\nlambda =", lambda, "\nmu = ", mu, "\nn =", n, "\nQ =", Q)

	//Potkaz2File, _ := os.Create("Potkaz2.txt")
	//defer Potkaz2File.Close()
	//
	//MfromN2File, _ := os.Create("MfromN2.txt")
	//defer MfromN2File.Close()
	//
	//Kz2File, _ := os.Create("Kz2.txt")
	//defer Kz2File.Close()
	//
	//PQ2File, _ := os.Create("PQ2.txt")
	//defer PQ2File.Close()
	//
	//MQ2File, _ := os.Create("MQ2.txt")
	//defer MQ2File.Close()
	//
	//KzQ2File, _ := os.Create("KzQ2.txt")
	//defer KzQ2File.Close()
	//
	//for n := 1; n <= 5; n++ {
	//	QStirng := fmt.Sprintf("%d", n)
	//
	//	Potkaz2File.WriteString(QStirng + "	")
	//	MfromN2File.WriteString(QStirng + "	")
	//	Kz2File.WriteString(QStirng + "	")
	//	PQ2File.WriteString(QStirng + "	")
	//	MQ2File.WriteString(QStirng + "	")
	//	KzQ2File.WriteString(QStirng + "	")
	//
	//	for Q := 1; Q <= 17; Q++ {
	//		K := make([]float64, n+Q)
	//
	//		K2 := Kfunc2(n, Q, lambda, mu, K)
	//		fmt.Println("\nK2: ", K2)
	//
	//		Po2 := PoCalc(K2, n+Q)
	//		fmt.Println("Po2: ", Po2)
	//
	//		Potkaz2 := PoCalc(K2, n+Q) * K2[n+Q-1]
	//		fmt.Println("Pотказ2: ", Potkaz2)
	//
	//		Potkaz2String := fmt.Sprintf("%f", Potkaz2)
	//		Potkaz2File.WriteString(Potkaz2String + "	")
	//
	//		M_N2 := MfromN2(K2, Po2, n, Q)
	//		fmt.Println("M(N)2: ", M_N2)
	//
	//		MfromN2String := fmt.Sprintf("%f", M_N2)
	//		MfromN2File.WriteString(MfromN2String + "	")
	//
	//		Kz2 := MfromN2(K2, Po2, n, Q) / float64(n)
	//		fmt.Println("Kz2: ", Kz2)
	//
	//		Kz2String := fmt.Sprintf("%f", Kz2)
	//		Kz2File.WriteString(Kz2String + "	")
	//
	//		PQ2 := PfromQ2(Po2, K2, lambda, mu, n, Q)
	//		fmt.Println("P(Q)2: ", PQ2)
	//
	//		PQ2String := fmt.Sprintf("%f", PQ2)
	//		PQ2File.WriteString(PQ2String + "	")
	//
	//		M_Q2 := MfromQ(Po2, K2, n, Q)
	//		fmt.Println("M(Q)2: ", M_Q2)
	//
	//		MQ2String := fmt.Sprintf("%f", M_Q2)
	//		MQ2File.WriteString(MQ2String + "	")
	//
	//		K_ZQ2 := MfromQ(Po2, K2, n, Q) / float64(n)
	//		fmt.Println("Kzq2: ", K_ZQ2)
	//
	//		KzQ2String := fmt.Sprintf("%f", K_ZQ2)
	//		KzQ2File.WriteString(KzQ2String + "	")
	//
	//	}
	//	Potkaz2File.WriteString("\n")
	//	MfromN2File.WriteString("\n")
	//	Kz2File.WriteString("\n")
	//	PQ2File.WriteString("\n")
	//	MQ2File.WriteString("\n")
	//	KzQ2File.WriteString("\n")
	//}

	// part 1.3
	//
	//MfromN3File, _ := os.Create("MfromN3.txt")
	//defer MfromN3File.Close()
	//
	//Kz3File, _ := os.Create("Kz3.txt")
	//defer Kz3File.Close()
	//
	//PQ3File, _ := os.Create("PQ3.txt")
	//defer PQ3File.Close()
	//
	//MQ3File, _ := os.Create("MQ3.txt")
	//defer MQ3File.Close()
	//
	//for n := 1; n <= 19; n++ {
	//	if geta(lambda, n, mu) >= 1 {
	//		continue
	//	}
	//
	//	nString := fmt.Sprintf("%d", n)
	//
	//	K3 := make([]float64, n)
	//	fmt.Println("\nlambda =", lambda, "\nmu = ", mu)
	//
	//	Kcomp3 := Kfunc(lambda, mu, n, K3)
	//	fmt.Println("\nK3: ", Kcomp3)
	//
	//	Po3 := PoCalc3(Kcomp3, n, lambda, mu)
	//	fmt.Println("Po3: ", Po3)
	//
	//	M_N3 := MformN3(Kcomp3, n, Po3, lambda, mu)
	//	fmt.Println("M(N)3: ", M_N3)
	//
	//	MfromN3String := fmt.Sprintf("%f", M_N3)
	//	MfromN3File.WriteString(nString + "	" + MfromN3String + "\n")
	//
	//	Kz3 := M_N3 / float64(n)
	//	fmt.Println("Kz3: ", Kz3)
	//
	//	Kz3String := fmt.Sprintf("%f", Kz3)
	//	Kz3File.WriteString(nString + "	" + Kz3String + "\n")
	//
	//	PQ3 := PfromQ3(K3, n, lambda, mu, Po3)
	//	fmt.Println("P(Q)3: ", PQ3)
	//
	//	PQ3String := fmt.Sprintf("%f", PQ3)
	//	PQ3File.WriteString(nString + "	" + PQ3String + "\n")
	//
	//	MQ3 := MfromQ3(K3, n, lambda, mu, Po3)
	//	fmt.Println("M(Q)3: ", MQ3)
	//
	//	MQ3String := fmt.Sprintf("%f", MQ3)
	//	MQ3File.WriteString(nString + "	" + MQ3String + "\n")
	//
	//	K_ZQ3 := MQ3 / float64(n)
	//	fmt.Println("Kzq3: ", K_ZQ3)
	//}

	//part 4

	MfromN4File, _ := os.Create("MfromN4.txt")
	defer MfromN4File.Close()

	Kz4File, _ := os.Create("Kz4.txt")
	defer Kz4File.Close()

	PQ4File, _ := os.Create("PQ4.txt")
	defer PQ4File.Close()

	MQ4File, _ := os.Create("MQ4.txt")
	defer MQ4File.Close()

	for n := 1; n <= 12; n++ {
		K4 := make([]float64, n)
		nString := fmt.Sprintf("%d", n)

		fmt.Println("\nlambda =", lambda, "\nmu = ", mu, "\nQ =", "\nV =", v)

		Kcomp4 := Kfunc(lambda, mu, n, K4)
		fmt.Println("\nK4: ", Kcomp4)

		Po4, ii := PoCalc4(K4, n, lambda, mu, v)
		fmt.Println("Po4: ", Po4, "\nii: ", ii)

		M_N4 := MfromN4(K4, n, Po4, lambda, mu, v, ii)
		fmt.Println("M(N)4: ", M_N4)

		MfromN4String := fmt.Sprintf("%f", M_N4)
		MfromN4File.WriteString(nString + "	" + MfromN4String + "\n")

		Kz4 := M_N4 / float64(n)
		fmt.Println("K(Z)4: ", Kz4)

		Kz4String := fmt.Sprintf("%f", Kz4)
		Kz4File.WriteString(nString + "	" + Kz4String + "\n")

		PQ4 := PfromQ4(K4, n, Po4, lambda, mu, v, ii)
		fmt.Println("P(Q)4: ", PQ4)

		PQ4String := fmt.Sprintf("%f", PQ4)
		PQ4File.WriteString(nString + "	" + PQ4String + "\n")

		MQ4 := MfromQ4(K4, n, Po4, lambda, mu, v, ii)
		fmt.Println("M(Q)4: ", MQ4)

		MQ4String := fmt.Sprintf("%f", MQ4)
		MQ4File.WriteString(nString + "	" + MQ4String + "\n")

		Kzq4 := MQ4 / float64(n)
		fmt.Println("Kzq4: ", Kzq4)
		fmt.Println()
	}
}

func Kfunc(lambda float64, mu float64, n int, K []float64) []float64 {
	K[0] = lambda / mu

	for i := 1; i < n; i++ {
		K[i] = K[i-1] * lambda / ((float64(i + 1)) * mu)
	}

	return K
}

func Kfunc2(n int, Q int, lambda float64, mu float64, K []float64) []float64 {
	var i int
	K[0] = lambda / mu

	for i = 1; i < n; i++ {
		K[i] = K[i-1] * lambda / ((float64(i + 1)) * mu)
	}

	for i = n; i < (n + Q); i++ {
		a := math.Pow(lambda, float64(i-n+1))
		b := math.Pow(float64(n)*mu, float64(i-n+1))
		A := a / b
		K[i] = K[n-1] * A
	}

	return K
}

func PoCalc(K []float64, n int) float64 {
	var sum float64

	for i := 0; i < n; i++ {
		sum += K[i]
	}

	Po := 1 / (1 + sum)

	return Po
}

func PoCalc3(K []float64, n int, lambda float64, mu float64) float64 {
	sum := 1.0

	for i := 0; i < n; i++ {
		sum += K[i]
	}

	return float64(1) / (sum + (K[n-1] * (geta(lambda, n, mu) / (float64(1) - geta(lambda, n, mu)))))
}

func PoCalc4(K []float64, n int, lambda float64, mu float64, v float64) (float64, int) {
	sum := 1.0
	var A float64
	var Poo float64
	var ii int
	var i int
	Po := 0.0
	accur := 0.000000001

	for i = 0; i < n; i++ {
		sum += K[i]
	}

	for i = 1; i > 0; i++ {
		if i == 1 {
			A = getA(lambda, mu, n, i, v)
		} else {
			A = A * getA(lambda, mu, n, i, v)
		}

		sum += K[n-1] * A
		Poo = Po
		Po = float64(1) / sum

		if math.Abs(Poo-Po) <= accur {
			ii = i
			return Po, ii
		}
	}

	return Po, ii
}

func MfromN(K []float64, Po float64, n int) float64 {
	var sum float64

	for i := 0; i < n; i++ {
		sum += K[i] * float64(i+1)
	}

	return Po * sum
}

func MfromN2(K []float64, Po float64, n int, Q int) float64 {
	var sum float64

	for i := 0; i < n+Q; i++ {
		if i < n {
			sum += K[i] * float64(i+1)
		} else {
			sum += K[i] * float64(n)
		}
	}

	return Po * sum
}

func MformN3(K []float64, n int, Po float64, lambda float64, mu float64) float64 {
	var sum float64

	for i := 0; i < n; i++ {
		sum += float64(i+1) * K[i]
	}

	return Po * (sum + float64(n)*(K[n-1]*(geta(lambda, n, mu)/(float64(1)-(geta(lambda, n, mu))))))
}

func MfromN4(K []float64, n int, Po float64, lambda float64, mu float64, v float64, ii int) float64 {
	var M_N, A, MN float64
	var i int

	for i = 0; i < n; i++ {
		M_N += float64(i+1) * K[i]
	}

	for i = 1; i < ii; i++ {
		if i == 1 {
			A = getA(lambda, mu, n, i, v)
		} else {
			A = A * getA(lambda, mu, n, i, v)
		}

		M_N += float64(n) * K[n-1] * A
		MN = Po * M_N
	}

	return MN
}

func PfromQ2(Po float64, K []float64, lambda float64, mu float64, n int, Q int) float64 {
	i := n
	var sum float64

	for i = n; i < n+Q; i++ {
		sum += math.Pow(geta(lambda, n, mu), float64(i-n+1))
	}

	return K[n-1] * Po * sum
}

func PfromQ3(K []float64, n int, lambda float64, mu float64, Po float64) float64 {
	return K[n-1] * Po * (geta(lambda, n, mu) / (float64(1) - geta(lambda, n, mu)))
}

func PfromQ4(K []float64, n int, Po float64, lambda float64, mu float64, v float64, ii int) float64 {
	var P_Q, PQ, A float64

	for i := 1; i < ii; i++ {
		if i == 1 {
			A = getA(lambda, mu, n, i, v)
		} else {
			A = A * getA(lambda, mu, n, i, v)
		}

		P_Q += A
		PQ = Po * K[n-1] * P_Q
	}

	return PQ
}

func MfromQ(Po float64, K []float64, n int, Q int) float64 {
	i := n
	var sum float64

	for i = n; i < n+Q; i++ {
		sum += float64(i-n+1) * K[i]
	}

	return Po * sum
}

func MfromQ3(K []float64, n int, lambda float64, mu float64, Po float64) float64 {
	return K[n-1] * Po * (geta(lambda, n, mu) / math.Pow((float64(1)-geta(lambda, n, mu)), 2.0))
}

func MfromQ4(K []float64, n int, Po float64, lambda float64, mu float64, v float64, ii int) float64 {
	var M_Q, MQ, A float64

	for i := 1; i < ii; i++ {
		if i == 1 {
			A = getA(lambda, mu, n, i, v)
		} else {
			A = A * getA(lambda, mu, n, i, v)
		}

		M_Q += float64(i) * A
		MQ = Po * K[n-1] * M_Q
	}

	return MQ
}

func getA(lambda float64, mu float64, n int, i int, v float64) float64 {
	return lambda / (float64(n)*mu + float64(i)*v)
}

func geta(lambda float64, n int, mu float64) float64 {
	return lambda / (float64(n) * mu)
}
