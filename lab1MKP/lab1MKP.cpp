#include <iostream>
#include <cmath>
#include <limits>
#define M_PI 3.14159265358979323846

double solver(double meanAnomaly, double e, double epsilon = 1e-9, int maxIterations = 1000) {
    //Начальное приближение E = M
    double E = meanAnomaly;

    for (int i = 0; i < maxIterations; ++i) {
        //Решение уравнения Кеплера
        double deltaE = (meanAnomaly - (E - e * sin(E))) /
            (1 - e * cos(E));
        E += deltaE;

        //Проверка сходимости
        if (fabs(deltaE) < epsilon) {
            break;
        }
    }

    return E;
}

double solver2(double meanAnomaly, double e, double epsilon = 1e-9, int maxIterations = 1000) {
    //Начальное приближение E = M
    double a = meanAnomaly - M_PI; //Любое значение < 0
    double b = meanAnomaly + M_PI; //Любое значение > 0

    double fa = a - e * sin(a) - meanAnomaly;
    double fb = b - e * sin(b) - meanAnomaly;

    if (fa * fb > 0) {
        std::cerr << "Ошибка: Начальные значения выбраны неправильно." << std::endl;
        return std::numeric_limits<double>::quiet_NaN(); //Возвращаем nan в случае ошибки
    }

    double c, fc;

    for (int i = 0; i < maxIterations; ++i) {
        c = (a + b) / 2.0;
        fc = c - e * sin(c) - meanAnomaly;

        if (fabs(b - a) < 2 * epsilon || fabs(fc) < epsilon) {
            break;
        }

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        }
        else {
            a = c;
            fa = fc;
        }
    }

    return c;
}

double solver3(double meanAnomaly, double e, double epsilon = 1e-9, int maxIterations = 1000) {

    const double chi = (1 + sqrt(5)) / 2;

    //Начальные значения
    double a = meanAnomaly - M_PI; //Любое значение < 0
    double b = meanAnomaly + M_PI; //Любое значение > 0

    double fa = a - e * sin(a) - meanAnomaly;
    double fb = b - e * sin(b) - meanAnomaly;

    if (fa * fb > 0) {
        std::cerr << "Ошибка: Начальные значения выбраны неправильно." << std::endl;
        return std::numeric_limits<double>::quiet_NaN(); //Возвращаем nan в случае ошибки
    }

    double c, fc;

    for (int i = 0; i < maxIterations; ++i) {
        //Разделение отрезка в соотношении золотого сечения
        c = b - (b - a) / chi;
        fc = c - e * sin(c) - meanAnomaly;

        if (fabs(b - a) < 2 * epsilon || fabs(fc) < epsilon) {
            break;
        }

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        }
        else {
            a = c;
            fa = fc;
        }
    }

    return c;
}

double solver4(double meanAnomaly, double e, double epsilon = 1e-9, int maxIterations = 1000) {
    double E = meanAnomaly; //Начальное приближение E = M

    for (int i = 0; i < maxIterations; ++i) {
        double f = E - e * sin(E) - meanAnomaly;
        double f_d = 1 - e * cos(E);

        //Метод Ньютона
        double deltaE = f / f_d;
        E -= deltaE;

        //Проверка сходимости
        if (fabs(deltaE) < epsilon) {
            break;
        }
    }

    return E;
}

int main() {
    setlocale(LC_ALL, "RU");
    double mass_mercury = 0.33e24;
    double r1 = 10300.0; // апоцентр
    double r2 = 200.0;   // перицентр
    double t = 6;        // время 

    //Вычисление периода обращения (T)
    double T = 2 * M_PI * sqrt(pow(r1, 3) / (6.67430e-11 * mass_mercury));

    //Вычисление средней аномалии (M)
    double M = 2 * M_PI * t / T;

    //Вычисление большой полуоси (a)
    double a = (r1 + r2) / 2.0;

    //Вычисление эксцентриситета (e)
    double e = (r1 - r2) / (r1 + r2);

    //Вычисление малой полуоси (b)
    double b = a * sqrt(1 - e * e);

    //Выбор метода
    int method;
    std::cout << "Выберите метод (1 - Итерации, 2 - Половинного деления , 3 - Золотого сечения, 4 - ньютона): ";
    std::cin >> method;

    double E;
    switch (method) {
    case 1:
        E = solver(M, e);
        break;
    case 2:
        E = solver2(M, e);
        break;
    case 3:
        E = solver3(M, e);
        break;
    case 4:
        E = solver4(M, e);
        break;
    default:
        std::cerr << "Ошибка: Некорректный метод." << std::endl;
        return 1;
    }

    if (!std::isnan(E)) {
        std::cout << "Средняя аномалия (M): " << M << " радиан" << std::endl;
        std::cout << "Эксцентриситет (e): " << e << std::endl;
        std::cout << "Эксцентрическая аномалия (E): " << E << " радиан" << std::endl;
    }

    return 0;
}
