#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits>

using namespace std;

double inf = std::numeric_limits<double>::infinity();

float fact(int N)
{
    if(N<0)
        return 0;
    if (N==0)
        return 1;
    else
        return N*fact(N-1);
}

float degree_number(float number, float degree)
{
    return (pow(number, degree));
}

//массив коэффициентов для пунктов 1, 3 и 4.
float *Pk_1(int n, float L, float mu)
{
    float *P;
    int i=1;

	//с факториалом.
    /*for(int i=0; i<n; i++)
    {
        P[i] = (pow(L, i+1) / (fact(i+1) * pow(mu, i+1)));
        cout << P[i] << "\t";
    }*/

	//без факториала.
    P[0] = L/mu;
    cout << "P" << i << "=" << P[0] << "\t";
    for(i; i<n; i++)
    {
        P[i] = P[i-1] * (L/(float(i+1)*mu));
        cout << "P" << i+1 << "=" << P[i] << "\t";
    }

    cout << "\n";
    return P;
}

//массив коэффициентов для пункта 2 (с определенной очередью).
float *Pk_2(int n, int m, float L, float mu)
{
    float *P, A, a, b;
    int i=1, k;

    P[0] = L/mu;
    cout << "P" << i << "=" << P[0] << "\t";
    for(i; i<n; i++)
    {
        P[i] = P[i-1] * (L/(float(i+1)*mu));
        cout << "P" << i+1 << "=" << P[i] << "\t";
    }

    k = i;
    for(i; i<(n+m); i++)
    {
        a = degree_number(L, float(i-k+1));
        b = degree_number((float(k)*mu), float(i-k+1));
        A = float(a)/float(b);
        P[i] = P[int(k)-1] * A;
        cout << "P" << k << "+" << i-k+1 << "=" << P[i] << "\t";
    }

    cout << "\n";
    return P;
}

//Ро для пунктов 1 и 2.
float Po(float *Pk, int n)
{
    float sum=1, Po;
    int i=0;

    for(i; i<n; i++)
        sum += Pk[i];
    Po = float(1)/sum;

    return Po;
}

float Po_3(float *Pk, int n, float L, float mu)
{
	if((L/(float(n)*mu)>1))
		return float(0);

    float sum=1, Po, a, A;
    int i=0;

    for(i; i<n; i++)
        sum += Pk[i];
    a = (L/(float(n)*mu));
    A = Pk[n-1] * (a/(float(1)-a));
    Po = float(1)/(sum + A);

    return Po;
}

//счетчик итераций достижения точности.
int ii;

float Po_4(float *Pk, int n, float L, float mu, float v)
{
    float sum=1, Po=0.0, accur=0.000000001, Poo, A;
    int i=0;

    for(i; i<n; i++)
        sum += Pk[i];

    for(i=1; i>0; i++)
	{
        if(i==1)
			A = (L/(float(n)*mu + float(i)*v));
        else
			A = A * (L/(float(n)*mu + float(i)*v));

        sum += Pk[n-1]*A;
        Poo = Po;
        Po = float(1)/sum;

        if(fabs(Poo - Po) <= accur)
		{
			ii = i;
            cout << "i=" << ii << " " << Po << endl;
            return Po;
        }
    }

    return Po;
}

float P_otk_1(float *Pk, int n, float Po, float L, float mu)
{
    float Potk = Po * Pk[n-1];
    return Potk;
}

float P_otk_2(float *Pk, int n, int m, float Po, float L, float mu)
{
    float Potk = Po * Pk[n+m-1];
    return Potk;
}

float M_N_1(float *Pk, int n, float Po)
{
    float M_N=0;

    for(int i=0; i<n; i++)
        M_N += float(i+1)*Pk[i];
    M_N = Po*M_N;

    return M_N;
}

float M_N_2(float *Pk, int n, int m, float Po, float L, float mu)
{
    float M_N=0;

    for(int i=0; i<(n+m); i++)
	{
        if(i<n)
			M_N += float(i+1)*Pk[i];
        else
            M_N += float(n)*Pk[i];
    }
    M_N = Po*M_N;

    return M_N;
}

float M_N_3(float *Pk, int n, float Po, float L, float mu)
{
	if((L/(float(n)*mu)>1))
		return float(n);

    float M_N=0, a, A;

    for(int i=0; i<n; i++)
        M_N += float(i+1)*Pk[i];
    a = (L/(float(n)*mu));
    A = Pk[n-1] * (a/(float(1)-a));
    M_N = Po*(M_N + float(n)*A);

    return M_N;
}

float M_N_4(float *Pk, int n, float Po, float L, float mu, float v, int j)
{
    float M_N=0, A, MN=0, MN_;
	int i;

    for(i=0; i<n; i++)
        M_N += float(i+1)*Pk[i];

    for(i=1; i<j; i++)
	{
        if(i==1)
			A = (L/(float(n)*mu + float(i)*v));
        else
			A = A * (L/(float(n)*mu + float(i)*v));

        M_N += float(n)*Pk[n-1]*A;
        MN_ = MN;
        MN = Po * M_N;
    }

    return MN;
}

float P_Q_2(float *Pk, int n, int m, float L, float mu, float Po)
{
    float sum=0, PQ, a, b, A;
    int i=n;

	//способ 1
    /*for(i; i<(n+m); i++)
	{
        sum += Pk[i];
        //cout << Pk[i] << endl;
    }
    //cout << sum << endl;
    PQ = Po * sum;*/

	//способ 2
    for(i; i<(n+m); i++)
	{
        a = degree_number(L, float(i-n+1));
        b = degree_number((float(n)*mu), float(i-n+1));
        A = float(a)/float(b);
        sum += A;
    }
    PQ = Pk[n-1] * Po * sum;

    return PQ;
}

float P_Q_3(float *Pk, int n, float L, float mu, float Po)
{
	if((L/(float(n)*mu)>1))
		return float(1);

    float sum=0, PQ, a, b, A;

    a = (L/(float(n)*mu));
    PQ = Pk[n-1] * Po * (a/(float(1)-a));

    return PQ;
}

float P_Q_4(float *Pk, int n, float Po, float L, float mu, float v, int j)
{
    float P_Q=0, A, PQ=0, PQ_;
	int i;

    for(i=1; i<j; i++)
	{
        if(i==1)
            A = (L/(float(n)*mu + float(i)*v));
        else
			A = A * (L/(float(n)*mu + float(i)*v));

        P_Q += A;
        PQ_ = PQ;
        PQ = Po * Pk[n-1] * P_Q;
    }

    return PQ;
}

float M_Q_2(float *Pk, int n, int m, float L, float mu, float Po)
{
    float sum=0, MQ, a, b, A;
    int i=n;

    for(i; i<(n+m); i++)
        sum += (float(i-n+1) * Pk[i]);
    MQ = Po * sum;

    return MQ;
}

float M_Q_3(float *Pk, int n, float L, float mu, float Po)
{
	if((L/(float(n)*mu)>1))
		return inf;

    float sum=0, MQ, a, b, A;

    a = (L/(float(n)*mu));
    MQ = Pk[n-1] * Po * (a/((float(1)-a)*(float(1)-a)));

    return MQ;
}

float M_Q_4(float *Pk, int n, float Po, float L, float mu, float v, int j)
{
    float M_Q=0, A, MQ=0, MQ_;
	int i;

    for(i=1; i<j; i++)
	{
        if(i==1)
            A = (L/(float(n)*mu + float(i)*v));
        else
			A = A * (L/(float(n)*mu + float(i)*v));

        M_Q += float(i) * A;
        MQ_ = MQ;
        MQ = Po * Pk[n-1] * M_Q;
    }

    return MQ;
}

int main()
{
    int n, m,
		Tc=30, Ts=217, Tw=538;

    float *Pk_v1, *Pk_v2,   L, mu,   a1, a2,   Kz1, Kz2,
		PQ2, MQ2,  Kzq2,   *Pk_v4, Po_v4,  v,   PQ4, MQ4,   a3, P_otk_v3,
		Po_v1, P_otk_v1, M_N_v1,   Po_v2, P_otk_v2, M_N_v2,
		*Pk_v3, Po_v3,  M_N_v3, Kz3,  PQ3, MQ3, Kzq3,   M_N_v4, Kz4, Kzq4;

    L=1/float(Tc);
	mu=1/float(Ts);
	v=1/float(Tw);
    //L=float(6); mu=float(1);

	//1.1.
    cout << "PART 1" << endl << endl;
    std::ofstream fout("1.1.Kz.txt");
    for(n=1; n<=14; n++)
    {
        //n=8;
        cout << "L=" << L << " " << "mu=" << mu << endl;
        cout << "N=" << n << endl;
        fout << std::setw(10) << std::left << n;
        Pk_v1 = Pk_1(n, L, mu);

        Po_v1 = Po(Pk_v1, n);
        cout << "Po_1=" << Po_v1 << endl;
        P_otk_v1 = P_otk_1(Pk_v1, n, Po_v1, L, mu);
        cout << "P_otk_1=" << P_otk_v1 << endl;
        M_N_v1 = M_N_1(Pk_v1, n, Po_v1);
        cout << "M_N_1=" << M_N_v1 << endl;
        a1 = L/(float(n)*mu);
        cout << "A_1=" << a1 << endl;
        Kz1 = M_N_v1/float(n);
        cout << "Kz_1=" << Kz1 << endl << endl;

        //fout << std::setw(10) << std::left << P_otk_v1 << endl;
        //fout << std::setw(10) << std::left << M_N_v1 << endl;
        fout << std::setw(10) << std::left << Kz1 << endl;
    }
    fout.close();

	//1.2.
    cout << "PART 2" << endl << endl;
    //m=10; n=8;
    std::ofstream fout("1.2.1.Potk.txt");
    //for(n=1; n<=9; n++)
	for(m=1; m<=19; m++)
    {
        //n=8; m=10;
		fout << std::setw(12) << std::left << m;
		//for(m=1; m<=10; m++)
		for(n=1; n<=9; n++)
		{
			cout << "L=" << L << " " << "mu=" << mu << endl;
			cout << "N=" << n << endl;
			cout << "M=" << m << endl;
			Pk_v2 = Pk_2(n, m, L, mu);

			Po_v2 = Po(Pk_v2, n+m);
			cout << "Po_2=" << Po_v2 << endl;
			P_otk_v2 = P_otk_2(Pk_v2, n, m, Po_v2, L, mu);
			cout << "P_otk_2=" << P_otk_v2 << endl;
			M_N_v2 = M_N_2(Pk_v2, n, m, Po_v2, L, mu);
			cout << "M_N_2=" << M_N_v2 << endl;
			a2 = L/(float(n)*mu);
			cout << "A_2=" << a2 << endl;
			Kz2 = M_N_v2/float(n);
			cout << "Kz_2=" << Kz2 << endl;
			PQ2 = P_Q_2(Pk_v2, n, m, L, mu, Po_v2);
			cout << "PQ_2=" << PQ2 << endl;
			MQ2 = M_Q_2(Pk_v2, n, m, L, mu, Po_v2);
			cout << "MQ_2=" << MQ2 << endl;
			Kzq2 = MQ2/float(n);
			cout << "Kzq_2=" << Kzq2 << endl << endl;

			//fout << std::setw(12) << std::left << P_otk_v2;
			//fout << std::setw(12) << std::left << M_N_v2;
			//fout << std::setw(12) << std::left << Kz2;
			//fout << std::setw(12) << std::left << PQ2;
			//fout << std::setw(12) << std::left << MQ2;
			fout << std::setw(12) << std::left << Kzq2;
		}
		fout << endl;
    }
    fout.close();

	//1.3.
    cout << "PART 3" << endl << endl;
    std::ofstream fout("1.3.MQ.txt");
    for(n=1; n<=14; n++)
    {
        //n=8;
        cout << "L=" << L << " " << "mu=" << mu << endl;
        cout << "N=" << n << endl;
        fout << std::setw(10) << std::left << n;
        Pk_v3 = Pk_1(n, L, mu);

        Po_v3 = Po_3(Pk_v3, n, L, mu);
        cout << "Po_3=" << Po_v3 << endl;
        M_N_v3 = M_N_3(Pk_v3, n, Po_v3, L, mu);
		cout << "M_N_3=" << M_N_v3 << endl;
		Kz3 = M_N_v3/float(n);
		cout << "Kz_3=" << Kz3 << endl;
        PQ3 = P_Q_3(Pk_v3, n, L, mu, Po_v3);
        cout << "PQ_3=" << PQ3 << endl;
        MQ3 = M_Q_3(Pk_v3, n, L, mu, Po_v3);
		cout << "MQ_3=" << MQ3 << endl;
		Kzq3 = MQ3/float(n);
		cout << "Kzq_3=" << Kzq3 << endl << endl;

        //fout << std::setw(10) << std::left << M_N_v3 << endl;
        //fout << std::setw(10) << std::left << Kz3 << endl;
        //fout << std::setw(10) << std::left << PQ3 << endl;
		fout << std::setw(10) << std::left << MQ3 << endl;
    }
    fout.close();

	//1.4.
    cout << "PART 4" << endl << endl;
    std::ofstream fout("1.4.MQ.txt");
    for(n=1; n<=14; n++)
    {
        //n=8; //v=1;
        cout << "L=" << L << " " << "mu=" << mu << endl;
        cout << "N=" << n << endl;
        cout << "V=" << v << endl;
        fout << std::setw(10) << std::left << n;
        Pk_v4 = Pk_1(n, L, mu);

        Po_v4 = Po_4(Pk_v4, n, L, mu, v);
        cout << "Po_4=" << Po_v4 << endl;
        M_N_v4 = M_N_4(Pk_v4, n, Po_v4, L, mu, v, ii);
        cout << "M_N_4=" << M_N_v4 << endl;
        Kz4 = M_N_v4/float(n);
        cout << "Kz_4=" << Kz4 << endl;
        PQ4 = P_Q_4(Pk_v4, n, Po_v4, L, mu, v, ii);
        cout << "PQ_4=" << PQ4 << endl;
        MQ4 = M_Q_4(Pk_v4, n, Po_v4, L, mu, v, ii);
        cout << "MQ_4=" << MQ4 << endl;
        Kzq4 = MQ4/float(n);
        cout << "Kzq_4=" << Kzq4 << endl << endl;

        //fout << std::setw(10) << std::left << M_N_v4 << endl;
        //fout << std::setw(10) << std::left << Kz4 << endl;
        //fout << std::setw(10) << std::left << PQ4 << endl;
		fout << std::setw(10) << std::left << MQ4 << endl;
    }
    fout.close();

    return 0;
}
