/*
BSD 3-Clause License
Copyright (c) 2020, MasaharuYoshihara
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <math.h>
#include <random>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <omp.h>

using namespace std;


const int x = 20;
const int y = 20;
const double dt = 0.0001;
const double betaN = 100;
const double betaD = 10;
const int betaR = 1000000;
const int kcinv = 10;
const int ktinv = 1;
const int gamma = 1;
//const int gammaR = 1;
const int kRS = 300000;
const int m = 1;
const int n = 3;
const int N = 4; //position of PV
const int dist = 3;//disturbance exponential
const int PV = pow(10, dist);
const int paramax = 100;

const int numIterMax = 2000000;

template <typename ArrayType>
void set_array2d(std::vector<std::vector<ArrayType>>* array, int d1, int d2) {

    (*array).resize(d1);
    for (int i = 0; i < d1; i++)
        (*array)[i].resize(d2);
}

template <typename ArrayType>
void set_array1d(std::vector<ArrayType>* array, int d3) {

    (*array).resize(d3);

}

void get_ave_sd(
    std::vector<std::vector<double>>* beforeR, //INPUT:2D-array of R
    vector <vector <int>>* skip,               //INPUT:Skip index of calculating ave and sd.
    double* ave,                               //OUTPUT:Average
    double* sd                                 //OUTPUT:SD
) {

    *ave = 0.0;
    for (int i = 0; i < (*beforeR).size(); i++)
        for (int j = 0; j < (*beforeR)[i].size(); j++)
            *ave += (*beforeR)[i][j];

    double askip = 0;
    for (int s = 0; s < (*skip).size(); s++)
        askip += (*beforeR)[(*skip)[s][0]][(*skip)[s][1]];

    *ave = (*ave - askip) / double(x * y - (*skip).size());

    double dev = 0;
    for (int i = 0; i < (*beforeR).size(); i++)
        for (int j = 0; j < (*beforeR)[i].size(); j++)
            dev += ((*beforeR)[i][j] - *ave) * ((*beforeR)[i][j] - *ave);

    double dskip = 0;
    for (int s = 0; s < (*skip).size(); s++)
        dskip += ((*beforeR)[(*skip)[s][0]][(*skip)[s][1]] - *ave) * ((*beforeR)[(*skip)[s][0]][(*skip)[s][1]] - *ave);

    double dev2 = dev - dskip;
    *sd = sqrt(dev2 / double(x * y - (*skip).size()));
}

int main() {
    FILE* fp;

    fopen_s(&fp, "SensitivityAnalysis_gammaR_1.txt", "w");

    if (fp == NULL) {
        printf("failed\n");
        return -1;
    }
    else {
        printf("creating file, please wait\n");
    }

    std::vector<double>
        diff31, diff32, diff33, diff34, diff, logdiff;
    set_array1d(&diff31, paramax);
    set_array1d(&diff32, paramax);
    set_array1d(&diff33, paramax);
    set_array1d(&diff34, paramax);
    set_array1d(&diff, paramax);
    set_array1d(&logdiff, paramax);

    std::vector<int> numIter;
    set_array1d(&numIter, paramax);

#pragma omp parallel for
    for (int k = 0; k < paramax + 1; k++) {

        int i;
        int j;

        double gammaR = pow(10, (0.1 * k) - 10);

        std::vector<std::vector<double>>
            beforeD, beforeN, nextD, nextN, aveD, aveN, deltaD, deltaN, beforeR, deltaR, nextR, aveR, pow1, pow2;

        std::vector<std::vector<int>>
            firstD, firstN, firstR;

        set_array2d(&beforeD, x, y);
        set_array2d(&beforeN, x, y);
        set_array2d(&beforeR, x, y);
        set_array2d(&nextD, x, y);
        set_array2d(&nextN, x, y);
        set_array2d(&nextR, x, y);
        set_array2d(&aveD, x, y);
        set_array2d(&aveN, x, y);
        set_array2d(&aveR, x, y);
        set_array2d(&firstD, x, y);
        set_array2d(&firstN, x, y);
        set_array2d(&firstR, x, y);
        set_array2d(&deltaD, x, y);
        set_array2d(&deltaN, x, y);
        set_array2d(&deltaR, x, y);
        set_array2d(&pow1, x, y);
        set_array2d(&pow2, x, y);


        /*
         for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                beforeD[i][j] = 0;
                beforeN[i][j] = 0;
                beforeR[i][j] = 0;
                nextD[i][j] = 0;
                nextN[i][j] = 0;
                nextR[i][j] = 0;
            }
        }
        */

        std::random_device seed_gen;
        std::mt19937 engine(100);
        std::normal_distribution<> dist(1.0, 0.1);
        for (i = 0; i < x; i++) {
            for (j = 0; j < y; j++) {
                beforeD[i][j] = dist(engine);
                beforeN[i][j] = dist(engine);
                beforeR[i][j] = 0.0;
            }
        }

        numIter[k] = 0;
        for (; ; ) {

            int sumfirstD = 0, sumfirstN = 0, sumfirstR = 0;

            for (i = 0; i < x; i++) {
                for (j = 0; j < y; j++) {
                    beforeD[N][N] = 1 * PV;

                    int ii1 = (i - 1 + x) % x;
                    int ii2 = (i + 1) % x;
                    int jj1 = (j - 1 + y) % y;
                    int jj2 = (j + 1) % y;

                    aveD[i][j] = (beforeD[i][jj2] + beforeD[ii2][j] + beforeD[i][jj1] + beforeD[ii1][j]) / 4;
                    aveN[i][j] = (beforeN[i][jj2] + beforeN[ii2][j] + beforeN[i][jj1] + beforeN[ii1][j]) / 4;
                    aveR[i][j] = (beforeR[i][jj2] + beforeR[ii2][j] + beforeR[i][jj1] + beforeR[ii1][j]) / 4;

                    pow1[i][j] = pow(beforeN[i][j] * aveD[i][j], n);
                    pow2[i][j] = pow(beforeR[i][j], m);

                    deltaN[i][j] = betaN - gamma * beforeN[i][j] - beforeN[i][j] * aveD[i][j] * ktinv - beforeN[i][j] * beforeD[i][j] * kcinv;
                    deltaR[i][j] = betaR * pow1[i][j] / (kRS + pow1[i][j]) - gammaR * beforeR[i][j];
                    deltaD[i][j] = betaD / (1 + pow2[i][j]) - gamma * beforeD[i][j] - beforeD[i][j] * aveN[i][j] * ktinv - beforeN[i][j] * beforeD[i][j] * kcinv;

                    nextD[i][j] = beforeD[i][j] + deltaD[i][j] * dt;
                    nextN[i][j] = beforeN[i][j] + deltaN[i][j] * dt;
                    nextR[i][j] = beforeR[i][j] + deltaR[i][j] * dt;

                    if (nextD[i][j] < 0.0)nextD[i][j] = 0.0;
                    if (nextN[i][j] < 0.0)nextN[i][j] = 0.0;
                    if (nextR[i][j] < 0.0)nextR[i][j] = 0.0;

                    firstD[i][j] = 0;
                    firstN[i][j] = 0;
                    firstR[i][j] = 0;

                    if (fabs(deltaD[i][j]) < 1e-3)firstD[i][j] = 1;
                    if (fabs(deltaN[i][j]) < 1e-3)firstN[i][j] = 1;
                    if (fabs(deltaR[i][j]) < 1e-3)firstR[i][j] = 1;

                    firstD[N][N] = 0;
                    firstN[N][N] = 0;
                    firstR[N][N] = 0;

                    sumfirstD += firstD[i][j];
                    sumfirstN += firstN[i][j];
                    sumfirstR += firstR[i][j];
                }
            }

            for (i = 0; i < x; i++) {
                for (j = 0; j < y; j++) {
                    beforeD[i][j] = nextD[i][j];
                    beforeN[i][j] = nextN[i][j];
                    beforeR[i][j] = nextR[i][j];
                }
            }

            numIter[k]++;

            //if (numIter % 10000 == 0 || numIter == 1)
            //fprintf(stdout, "Iter%10d deltaD=%10.2e deltaN=%10.2e deltaR=%10.2e\n", numIter, deltaD[N + 3][N + 3], deltaN[N + 3][N + 3], deltaR[N + 3][N + 3]);

            if ((sumfirstD == 399 && sumfirstN == 399 && sumfirstR == 399) || numIter[k] == numIterMax) {
                //fprintf(stdout, "Iter%10d deltaD=%10.2e deltaN=%10.2e deltaR=%10.2e\n", numIter, deltaD[N + 3][N + 3], deltaN[N + 3][N + 3], deltaR[N + 3][N + 3]);
                break;
            }
        }

        vector <vector <int>> skip3 = { { N, N }, { N - 1, N }, { N + 1, N }, { N, N - 1 }, { N, N + 1 } };
        double ave3, sd3;
        get_ave_sd(&beforeR, &skip3, &ave3, &sd3);

        diff31[k] = (beforeR[N - 1][N] - ave3) / sd3;
        diff32[k] = (beforeR[N + 1][N] - ave3) / sd3;
        diff33[k] = (beforeR[N][N - 1] - ave3) / sd3;
        diff34[k] = (beforeR[N][N + 1] - ave3) / sd3;

        diff[k] = (diff31[k] + diff32[k] + diff33[k] + diff34[k]) / 4;
        logdiff[k] = log2(diff[k]);

        fprintf(stdout, "%d %.2f %7d \n", k, logdiff[k], numIter[k]);//changed by tnagata. numIterも各スレッドの結果を格納する配列としました

    }

    for (int k = 0; k < paramax + 1; k++)
        fprintf(fp, "%d %.2f %7d \n", k, logdiff[k], numIter[k]);


    fclose(fp);
    return 0;
}