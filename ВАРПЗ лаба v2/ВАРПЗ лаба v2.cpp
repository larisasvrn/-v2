#include "Header.h"

double f(double x, double u, double v);
double u_right(double x);
double v_right(double x);

int main() {
	/*
	* input:
	* a b
	* n
	* Eps
	* K
	* alpha0 alpha1
	* beta0 beta1
	* gamma0 gamma1
	* A=v(a), B=u(b)
	*/
	ifstream input("input2.txt");
	/*
	* output:
	* IER(0 - ошибок нет, выводим решение. 1 - превышенно число итераций, решения нет. 2 - деление на ноль, решение нет
	* alpha - найденное значение u(a)
	* beta - найденное значение u((a+b)/2)
	* alpha - найденное значение v((a+b)/2)
	* a u(a) |u(a)-u_right(a)| v(a) |v(a)-v_right(a)|
	* a+h u(a+h) |u(a+h)-u_right(a+h)| v(a+h) |v(a+h)-v_right(a+h)|
	* ...
	* a+i*h u(a+i*h) |u(a+i*h)-u_right(a+i*h)| v(a+i*h) |v(a+i*h)-v_right(a+i*h)|	xi=a+i*h - очередная точка сетки, u(xi) - найденное значение функции, |u(xi - u_right(xi)| - погрешность функции в точке xi, v(xi) - найденное значение производной функции, |v(xi - v_right(xi)| - погрешность производной функции в точке xi
	* ...
	* b u(b) |u(b)-u_right(b)| v(b) |v(b)-v_right(b)|
	*/
	ofstream output("output.txt");
	double a, b, Eps, alpha0, alpha1, beta0, beta1, gamma0, gamma1, A, B, * u, * v;
	int n, L, K, IER;
	input >> a >> b >> n >> Eps >> K >> alpha0 >> alpha1 >> beta0 >> beta1 >> gamma0 >> gamma1 >> A >> B;
	u = new double[n + 1];
	v = new double[n + 1];
	shooting(a, b, n, u, v, f, A, B, alpha0, alpha1, beta0, beta1, gamma0, gamma1, Eps, K, L, IER);
	if (IER == 0) {
		output << IER << endl
			<< L << endl
			<< alpha0 << '\t' << beta0 << '\t' << gamma0 << endl;
		double h = (b - a) / n;
		for (int i = 0; i <= n; i++) {
			double x = a + i * h;
			output << x << '\t' << u[i] << '\t' << abs(u[i] - u_right(x)) << '\t' << v[i] << '\t' << abs(v[i] - v_right(x)) << endl;
		}
	}
	else {
		output << IER;
	}
	input.close();
	output.close();
}

/*
double f(double x, double u, double v) {
	return 3 * u + v - 3 * x * x * x * x - 7 * x * x * x + 12 * x * x + 11 * x - 4;
}
double u_right(double x) {
	return x * x * x * x + x * x * x - x * x - x + 1;
}
double v_right(double x) {
	return 4 * x * x * x + 3 * x * x - 2 * x - 1;
}*/

double f(double x, double u, double v) {
	return 6 * u + v + v * v - 4 * exp(-4 * x);
}
double u_right(double x) {
	return exp(-2 * x);
}
double v_right(double x) {
	return -2 * exp(-2 * x);
}