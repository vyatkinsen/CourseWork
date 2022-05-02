#include <iostream>
#include <cmath>
#include "zeroin.h"
#include "matrix.h"
#include "quanc8.h"

using namespace std;

double helperFunctionForQUANC8(double x) {
    return exp(-x) / pow(x, 3);
}

double calculateEpsilon() {
    double RESULT, FLAG, ERREST, RELERR = 0, ABSERR = 1.0e-10; //вспомогательные переменные для QUANC8
    int NOFUN;
    double downLimit = 0.1;
    double upLimit = 0.2;
    double constanta = 0.03045946;
    quanc8(helperFunctionForQUANC8, downLimit, upLimit,
           ABSERR, RELERR, &RESULT, &ERREST, &NOFUN, &FLAG);
//    Раскомментируйте следующую строку, если нужно округлить значение корня до 2-ух знаков после запятой
//    RESULT = floor(RESULT * 100) / 100;
//    Раскомментируйте следующую строку, если нужно округлить значение константы до 2-ух знаков после запятой
//    constanta = floor(constanta * 100) / 100;
//    Раскомментируйте следующую строку, если нужно округлить вниз значение интеграла
    RESULT = floor(RESULT);
//    Раскомментируйте следующую строку, если нужно округлить вверх значение интеграла
//    RESULT = ceil(RESULT);
    return constanta * RESULT;
}

double helperFunctionForZEROIN(double x) {
    return sin(x) - pow(cos(x), 2) + 0.25;
}

double calculateL() {
    double a = 0.0, b = 100.0, tol = 1e-8; //вспомогательные переменные для ZEROIN
    double minPositiveRoot = Zeroin(helperFunctionForZEROIN, a, b, tol);
    double constanta = 1.909859;
//    Раскомментируйте следующую строку, если нужно округлить значение корня до 2-ух знаков после запятой
//    minPositiveRoot = floor(minPositiveRoot * 100) / 100;
//    Раскомментируйте следующую строку, если нужно округлить значение константы до 2-ух знаков после запятой
//    constanta = floor(constanta * 100) / 100;
    return minPositiveRoot * constanta;
}

void makeAllCalculations(int n, double eps, double h, double **values) {
    MATRIX(system);
    VECTOR(answers, n - 1);
    int ipvt[n - 1]; //вспомогательный массив для SOLVE
    double cond[n - 1], work[n - 1]; //вспомогательные массивы для DECOMP
    double T[n + 1];

    for (int i = 0; i < n - 1; i++) for (int j = 0; j < n - 1; j++) system[i][j] = 0.0;
    system[0][0] = -2 * eps - pow(h, 3);
    system[0][1] = eps + pow(h, 2) / 2;
    for (int i = 1; i < n - 2; i++) {
        system[i][i - 1] = eps - (pow(h, 2) * (i + 1) / 2);
        system[i][i] = -2 * eps - pow(h, 3) * (i + 1);
        system[i][i + 1] = eps + pow(h, 2) * (i + 1) / 2;
    }
    system[n - 2][n - 3] = eps - pow(h, 2) * (n - 1) / 2;
    system[n - 2][n - 2] = -2 * eps - pow(h, 3) * (n - 1);

    for (int i = 0; i < n; i++) answers[i] = 0;
    answers[n - 2] = -(eps + pow(h, 2) * (n - 1) / 2);

    decomp(n - 1, system, cond, ipvt, work);
    solve(n - 1, system, answers, ipvt);

    for (int i = 0;i < n; i++) T[i + 1] = answers[i];
    T[0] = 0;
    T[n] = 1;
    *values = T;
}

void printResultForEach(double *values, double eps, double l, double h, int n) {
    cout << "-----------------------------------------------";
    printf("\neps = %.2f; l = %.2f; h = %.2f; n = %.2d\n", eps, l, h, n);
    double counter = 0.0;
    int i = 0;
    l = l + 0.5 - (l < 0);
    int s = (int) l;
    printf("\nx\t\tSOLVE T(x)\n");
    while (counter <= s) {
        printf("%.2f\t%.5f\n", counter, values[i++]);
        counter += h;
    }
}

void printResultForAll(double *firstValues, double *secondValues, double eps, double l, double h) {
    cout << "-----------------------------------------------";
    printf("\neps = %.2f;\t\tl = %.2f;\t\th = %.2f;\n", eps, l, h);
    double counter = 0.0;
    int i = 0, y = 0;
    l = l + 0.5 - (l < 0);
    int s = (int) l;
    printf("\nx\t\tT(x)_n=10\tT(x)_n=20\tDIFF\n");
    while (counter <= s) {
        printf("%.2f\t%.5f\t\t%.5f\t\t%.5f\n", counter, firstValues[i],
               secondValues[y], firstValues[i] - secondValues[y]);
        i++;
        y += 2;
        counter += h;
    }
}


int main() {
    double eps = calculateEpsilon();
    double l = calculateL();
    int n_1 = 10;
    int n_2 = 20;
    double h_1 = l / n_1;
    double h_2 = l / n_2;

    double firstArray[n_1 + 1];
    double secondArray[n_2 + 1];
    double *firstArrayPointer;
    double *secondArrayPointer;

    makeAllCalculations(n_1, eps, h_1, &firstArrayPointer);
    for (int x = 0; x <= n_1; x++) firstArray[x] = firstArrayPointer[x];
    printResultForEach(firstArray, eps, l, h_1, n_1);

    makeAllCalculations(n_2, eps, h_2, &secondArrayPointer);
    for (int x = 0; x <= n_2; x++) secondArray[x] = secondArrayPointer[x];
    printResultForEach(secondArray, eps, l, h_2, n_2);

    printResultForAll(firstArray, secondArray, eps, l, h_1);
    return 0;
}