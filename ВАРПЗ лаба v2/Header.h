#include <cstdlib>
#include <fstream>
using namespace std;

// Максим лох

/*
* Метод Рунге-кутта
* K1 = h * f(xi, yi)
* K2 = h * f(xi + h/2, yi + K1/2)
* K3 = h * f(xi + h/2, yi + K2/2)
* K4 = h * f(xi + h, yi + K3)
* yi+1 = yi + (K1 + 2 * K2 + 2 * K3 + K4)/6
*
* сделан для решения системы дифуров с двумия переменными
*/

void RK_solve(double a, double b, int n, double* u, double* v, double(*f1)(double, double, double), double(*f2)(double, double, double)) {
	double h = (b - a) / n;
	for (int i = 0; i < n; i++) {
		double
			x = a + i * h,
			K1u = h * f1(x, u[i], v[i]),
			K1v = h * f2(x, u[i], v[i]),
			K2u = h * f1(x + h * 0.5, u[i] + K1u * 0.5, v[i] + K1v * 0.5),
			K2v = h * f2(x + h * 0.5, u[i] + K1u * 0.5, v[i] + K1v * 0.5),
			K3u = h * f1(x + h * 0.5, u[i] + K2u * 0.5, v[i] + K2v * 0.5),
			K3v = h * f2(x + h * 0.5, u[i] + K2u * 0.5, v[i] + K2v * 0.5),
			K4u = h * f1(x + h, u[i] + K3u, v[i] + K3v),
			K4v = h * f2(x + h, u[i] + K3u, v[i] + K3v);
		u[i + 1] = u[i] + (K1u + 2 * K2u + 2 * K3u + K4u) / 6;
		v[i + 1] = v[i] + (K1v + 2 * K2v + 2 * K3v + K4v) / 6;
	}
}

void shooting(double a, double b, int n, double*& u, double*& v, double (*f)(double, double, double), double A, double B, double& alpha0, double alpha1, double& beta0, double beta1, double& gamma0, double gamma1, double Eps, int K, int& L, int& IER) {
	double c = (a + b) * 0.5;
	double
		* u1_0 = new double[n / 2 + 1],
		* v1_0 = new double[n / 2 + 1],
		* u2_0 = new double[n / 2 + 1],
		* v2_0 = new double[n / 2 + 1];
	u = new double[n];
	v = new double[n];

	//заполняем начальные значения для u1. u1(a) = alpha0, v1(a) = A
	u1_0[0] = alpha0;
	v1_0[0] = A;

	//заполняем начальные значения для u2. u2(c) = beta0, v2(c) = gamma0
	u2_0[0] = beta0;
	v2_0[0] = gamma0;

	IER = 0;

	//решаем задачу Коши для u1 на отрезке [a, c] и для на отрезке [c, b]
	RK_solve(a, c, n / 2, u1_0, v1_0, [](double x, double u, double v) {return v; }, f);
	RK_solve(c, b, n / 2, u2_0, v2_0, [](double x, double u, double v) {return v; }, f);

	L = 0;

	/*
	* Проверяем условие
	* u1(c) = beta
	* u1(c) = gamma
	* u2(b) = B
	* с точностью до Eps
	*/
	if (abs(u1_0[n / 2] - beta0) <= Eps && abs(v1_0[n / 2] - gamma0) <= Eps && abs(u2_0[n / 2] - B) <= Eps) {
		/*
		* Если условие выполнено, переписываем значение u1, u2 в u
		* u(x) = u1(x) a <= x < c
		* u(c) = (u1(c) + u2(c))/2
		* u(x) = u2(x) c < x <=b
		*/
		for (int i = 0; i < n / 2; i++) {
			u[i] = u1_0[i];
			v[i] = v1_0[i];
		}
		u[n / 2] = (u1_0[n / 2] + u2_0[0]) * 0.5;
		v[n / 2] = (v1_0[n / 2] + v2_0[0]) * 0.5;
		for (int i = 1; i <= n / 2; i++) {
			u[i + n / 2] = u2_0[i];
			v[i + n / 2] = v2_0[i];
		}
		return;
	}
	double
		* u1_1 = new double[n / 2 + 1],
		* v1_1 = new double[n / 2 + 1],
		* u2_1 = new double[n / 2 + 1],
		* v2_1 = new double[n / 2 + 1];

	//заполняем начальные значения для u1: u1(a) = alpha1, v1(a) = A, для u2: u2(c) = beta1, v2(c) = gamma1
	u1_1[0] = alpha1;
	v1_1[0] = A;
	u2_1[0] = beta1;
	v2_1[0] = gamma1;
	RK_solve(a, c, n / 2, u1_1, v1_1, [](double x, double u, double v) {return v; }, f);
	RK_solve(c, b, n / 2, u2_1, v2_1, [](double x, double u, double v) {return v; }, f);
	L = 1;

	//проверяем условие при начальных значениях alpha1, beta1, gamma1
	if (abs(u1_1[n / 2] - beta1) <= Eps && abs(v1_1[n / 2] - gamma1) <= Eps && abs(u2_1[n / 2] - B) <= Eps) {
		for (int i = 0; i < n / 2; i++) {
			u[i] = u1_1[i];
			v[i] = v1_1[i];
		}
		u[n / 2] = (u1_1[n / 2] + u2_1[0]) * 0.5;
		v[n / 2] = (v1_1[n / 2] + v2_1[0]) * 0.5;
		for (int i = 1; i <= n / 2; i++) {
			u[i + n / 2] = u2_1[i];
			v[i + n / 2] = v2_1[i];
		}
		alpha0 = alpha1;
		beta0 = beta1;
		gamma0 = gamma1;
		return;
	}

	//пока число итераций не превысильно максимальное заданное
	while (L < K) {
		double
			* u2_01 = new double[n / 2 + 1],
			* v2_01 = new double[n / 2 + 1],
			* u2_10 = new double[n / 2 + 1],
			* v2_10 = new double[n / 2 + 1];

		//заполняем смешанные начальный условия для u2. (u1 зависит только от alpha, поэтому для него нет смешанного начального условия)
		u2_01[0] = beta0;
		v2_01[0] = gamma1;
		RK_solve(c, b, n / 2, u2_01, v2_01, [](double x, double u, double v) {return v; }, f);

		//проверяем u2(beta0, gamma1) = B
		if (abs(u2_01[n / 2] - B) <= Eps) {

			/*
			* проверяем условия
			* u1(alpha0) = beta0
			* v1(alpha0) = gamma1
			*/
			if (abs(u1_0[n / 2] - beta0) <= Eps && abs(v1_0[n / 2] - gamma1) <= Eps) {

				//если условие выполнено, то первую часть переписываем из решения для alpha0, вторую переписываем из решения для смешанного условия
				for (int i = 0; i < n / 2; i++) {
					u[i] = u1_0[i];
					v[i] = v1_0[i];
				}
				u[n / 2] = (u1_0[n / 2] + u2_01[0]) * 0.5;
				v[n / 2] = (v1_0[n / 2] + v2_01[0]) * 0.5;
				for (int i = 1; i <= n / 2; i++) {
					u[i + n / 2] = u2_01[i];
					v[i + n / 2] = v2_01[i];
				}

				/*
				* найденные значения переписываются в alpha0, beta0, gamma0
				* alpha = alpha0
				* beta = beta0
				* gamma = gamma1
				*/
				gamma0 = gamma1;
				return;
			}

			/*
			* проверяем условия
			* u1(alpha1) = beta0
			* v1(alpha1) = gamma1
			*/
			else if (abs(u1_1[n / 2] - beta0) <= Eps && abs(v1_1[n / 2] - gamma1) <= Eps) {
				for (int i = 0; i < n / 2; i++) {
					u[i] = u1_1[i];
					v[i] = v1_1[i];
				}
				u[n / 2] = (u1_1[n / 2] + u2_01[0]) * 0.5;
				v[n / 2] = (v1_1[n / 2] + v2_01[0]) * 0.5;
				for (int i = 1; i <= n / 2; i++) {
					u[i + n / 2] = u2_01[i];
					v[i + n / 2] = v2_01[i];
				}

				/*
				* alpha = alpha1
				* beta = beta0
				* gamma = gamma1
				*/
				alpha0 = alpha1;
				gamma0 = gamma1;
				return;
			}
		}

		//заполняем смешанное условие для u2 вторым способом
		u2_10[0] = beta1;
		v2_10[0] = gamma0;
		RK_solve(c, b, n / 2, u2_10, v2_10, [](double x, double u, double v) {return v; }, f);

		//проверяем u2(beta1, gamma0) = B
		if (abs(u2_10[n / 2] - B) <= Eps) {


			/*
			* проверяем условия
			* u1(alpha0) = beta1
			* v1(alpha0) = gamma0
			*/
			if (abs(u1_0[n / 2] - beta1) <= Eps && abs(v1_0[n / 2] - gamma0) <= Eps) {
				for (int i = 0; i < n / 2; i++) {
					u[i] = u1_0[i];
					v[i] = v1_0[i];
				}
				u[n / 2] = (u1_0[n / 2] + u2_10[0]) * 0.5;
				v[n / 2] = (v1_0[n / 2] + v2_10[0]) * 0.5;
				for (int i = 1; i <= n / 2; i++) {
					u[i + n / 2] = u2_10[i];
					v[i + n / 2] = v2_10[i];
				}

				/*
				* alpha = alpha0
				* beta = beta1
				* gamma = gamma0
				*/
				beta0 = beta1;
				return;
			}

			/*
			* проверяем условия
			* u1(alpha0) = beta0
			* v1(alpha0) = gamma1
			*/
			else if (abs(u1_1[n / 2] - beta1) <= Eps && abs(v1_1[n / 2] - gamma0) <= Eps) {
				for (int i = 0; i < n / 2; i++) {
					u[i] = u1_1[i];
					v[i] = v1_1[i];
				}
				u[n / 2] = (u1_1[n / 2] + u2_10[0]) * 0.5;
				v[n / 2] = (v1_1[n / 2] + v2_10[0]) * 0.5;
				for (int i = 1; i <= n / 2; i++) {
					u[i + n / 2] = u2_10[i];
					v[i + n / 2] = v2_10[i];
				}

				/*
				* alpha = alpha1
				* beta = beta1
				* gamma = gamma0
				*/
				alpha0 = alpha1;
				beta0 = beta1;
				return;
			}
		}
		if (abs(alpha1 - alpha0) < 1.e-16 || abs(beta1 - beta0) < 1.e-16 || abs(gamma1 - gamma0) < 1.e-16) {
			IER = 2;
			return;
		}
		double
			//выписываем производные функций phi по alpha, beta, gamma (для остальных 5 из 9 можно найти точное значение)
			dphi1_alpha = (u1_1[n / 2] - u1_0[n / 2]) / (alpha1 - alpha0),
			dphi2_alpha = (v1_1[n / 2] - v1_0[n / 2]) / (alpha1 - alpha0),
			dphi3_beta = (u2_1[n / 2] - u2_01[n / 2]) / (beta1 - beta0),
			dphi3_gamma = (u2_1[n / 2] - u2_10[n / 2]) / (gamma1 - gamma0),

			//считаем определитель левой части системы уравнений относительно delta_alpha, delta_beta, delta_gamma
			delta = dphi1_alpha * dphi3_beta + dphi2_alpha * dphi3_gamma;
		if (abs(delta) < 1.e-16) {
			IER = 2;
			return;
		}

		//считаем коэффициенты для метода Крамера
		double
			delta1 = -(u2_1[n / 2] - B) - (u1_1[n / 2] - beta1) * dphi3_beta - (v1_1[n / 2] - gamma1) * dphi3_gamma,
			delta2 = -(v1_1[n / 2] - gamma1) * dphi1_alpha * dphi3_gamma - (u2_1[n / 2] - B) * dphi1_alpha + (u1_1[n / 2] - beta1) * dphi2_alpha * dphi3_gamma,
			delta3 = -(u1_1[n / 2] - beta1) * dphi2_alpha * dphi3_beta + (v1_1[n / 2] - gamma1) * dphi1_alpha * dphi3_beta - (u2_1[n / 2] - B) * dphi2_alpha;

		/*
		* переприсваиваем
		* (alpha0, beta0, gamma0) = (alpha1, beta1, gamma1)
		* значения из соотвествующих массивов, чтобы не считать заново решение для alpha0, beta0, gamma0
		*/
		delete[] u1_0;
		delete[] v1_0;
		delete[] u2_0;
		delete[] v2_0;
		delete[] u2_01;
		delete[] v2_01;
		delete[] u2_10;
		delete[] v2_10;
		u1_0 = u1_1;
		v1_0 = v1_1;
		u2_0 = u2_1;
		v2_0 = v2_1;
		u1_1 = new double[n / 2 + 1];
		v1_1 = new double[n / 2 + 1];
		u2_1 = new double[n / 2 + 1];
		v2_1 = new double[n / 2 + 1];
		alpha0 = alpha1;
		beta0 = beta1;
		gamma0 = gamma1;

		//считаем новые alpha1, beta1, gamma1
		alpha1 = alpha0 + delta1 / delta;
		beta1 = beta0 + delta2 / delta;
		gamma1 = gamma0 + delta3 / delta;

		//заполняем новые начальные условия и решаем задачи Коши
		u1_1[0] = alpha1;
		v1_1[0] = A;
		u2_1[0] = beta1;
		v2_1[0] = gamma1;
		RK_solve(a, c, n / 2, u1_1, v1_1, [](double x, double u, double v) {return v; }, f);
		RK_solve(c, b, n / 2, u2_1, v2_1, [](double x, double u, double v) {return v; }, f);
		L++;
		if (abs(u1_1[n / 2] - beta1) <= Eps && abs(v1_1[n / 2] - gamma1) <= Eps && abs(u2_1[n / 2] - B) <= Eps) {
			for (int i = 0; i < n / 2; i++) {
				u[i] = u1_1[i];
				v[i] = v1_1[i];
			}
			u[n / 2] = (u1_1[n / 2] + u2_1[0]) * 0.5;
			v[n / 2] = (v1_1[n / 2] + v2_1[0]) * 0.5;
			for (int i = 1; i <= n / 2; i++) {
				u[i + n / 2] = u2_1[i];
				v[i + n / 2] = v2_1[i];
			}
			alpha0 = alpha1;
			beta0 = beta1;
			gamma0 = gamma1;
		}
	}
	IER = 1;
}