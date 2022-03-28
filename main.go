package main

import (
	"fmt"
	"math"
)

func main() {
	//Tc := 49
	//Ts := 213
	//Tw := 587
	//n := 4

	lambda := 3.0
	mu := 1.0
	//part 1.1
	for n := 1; n <= 4; n++ {
		K := make([]float64, n)
		Kcomp := Kfunc(lambda, mu, n, K)
		fmt.Println("K1: ", Kcomp)
		Po := PoCalc(K, n)
		fmt.Println("Po: ", Po)
		Potkaz := PoCalc(K, n) * K[n-1]
		fmt.Println("Pотказ: ", Potkaz)
		M_N := MfromN(K, Po, n)
		fmt.Println("M(N): ", M_N)
		K_Z := M_N / float64(n)
		fmt.Println("Kz: ", K_Z)
	}

	//part 1.2
	Q := 2
	n := 4
	K := make([]float64, n+Q)
	K2 := Kfunc2(n, Q, lambda, mu, K)
	fmt.Println("K2: ", K2)
	Po2 := PoCalc(K2, n+Q)
	fmt.Println("Po2", Po2)
	Potkaz2 := PoCalc(K2, n+Q) * K2[n+Q-1]
	fmt.Println("Pотказ2: ", Potkaz2)
	M_N2 := MfromN2(K2, Po2, n, Q)
	fmt.Println("M(N)2: ", M_N2)
	K_Z2 := MfromN2(K2, Po2, n, Q) / float64(n)
	fmt.Println("Kz2: ", K_Z2)
	PQ2 := P_Q2(Po2, K2, lambda, mu, n, Q)
	fmt.Println("P(Q)2: ", PQ2)
	M_Q2 := MfromQ(Po2, K2, n, Q)
	fmt.Println("M(Q)2: ", M_Q2)
	K_ZQ2 := MfromQ(Po2, K2, n, Q) / float64(n)
	fmt.Println("Kzq2: ", K_ZQ2)

	// part 1.3
	K3 := make([]float64, n)
	Kcomp3 := Kfunc(lambda, mu, n, K3)
	fmt.Println(Kcomp3)
	Po3 := PoCalc3(K3, n, lambda, mu)
	fmt.Println(Po3)
	M_N3 := MformN3(K3, n, Po3, lambda, mu)
	fmt.Println(M_N3)
	K_Z3 := M_N3 / float64(n)
	fmt.Println(K_Z3)
	PQ3 := PQ_3(K3, n, lambda, mu, Po3)
	fmt.Println(PQ3)
	MQ3 := MQ_3(K3, n, lambda, mu, Po3)
	fmt.Println(MQ3)
	K_ZQ3 := MQ3 / float64(n)
	fmt.Println(K_ZQ3)
}

func Kfunc(lambda float64, mu float64, n int, K []float64) []float64 {
	K[0] = lambda / mu
	for i := 1; i < n; i++ {
		K[i] = K[i-1] * lambda / ((float64(i + 1)) * mu)
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

func MfromN(K []float64, Po float64, n int) float64 {
	var sum float64
	for i := 0; i < n; i++ {
		sum += K[i] * float64(i+1)
	}
	return Po * sum
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

func P_Q2(Po float64, K []float64, lambda float64, mu float64, n int, Q int) float64 {
	i := n
	var sum float64
	for i = n; i < n+Q; i++ {
		sum += math.Pow((lambda / (float64(n) * mu)), float64(i-n+1))
	}

	return K[n-1] * Po * sum
}

func MfromQ(Po float64, K []float64, n int, Q int) float64 {
	i := n
	var sum float64
	for i = n; i < n+Q; i++ {
		sum += float64(i-n+1) * K[i]
	}
	return Po * sum
}

func PoCalc3(K []float64, n int, lambda float64, mu float64) float64 {
	var sum float64
	for i := 0; i < n; i++ {
		sum += K[i]
	}
	return float64(1) / (sum + (K[n-1] * (lambda / (float64(n) * mu) / (float64(1) - (lambda / float64(n) * mu)))))
}

func MformN3(K []float64, n int, Po float64, lambda float64, mu float64) float64 {
	var sum float64
	for i := 0; i < n; i++ {
		sum += float64(i+1) * K[i]
	}
	return Po * (sum + float64(n)*(K[n-1])*(lambda/(float64(n)*mu)/(float64(1)-(lambda/float64(n)*mu))))
}

func PQ_3(K []float64, n int, lambda float64, mu float64, Po float64) float64 {
	return K[n-1] * Po * (lambda / (float64(n) * mu) / (float64(1) - (lambda / float64(n) * mu)))
}

func MQ_3(K []float64, n int, lambda float64, mu float64, Po float64) float64 {
	return K[n-1] * Po * K[n-1] * Po * (lambda / (float64(n) * mu) / math.Pow((float64(1)-(lambda/float64(n)*mu)), 2.0))
}
