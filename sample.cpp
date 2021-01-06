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

using namespace std;

int i;
int j;
const int x = 20;
const int y = 20;
const double dt = 0.0001;
const int betaN = 100;
const int betaD = 10;
const int betaR = 1000000;
const int kcinv = 10;
const int ktinv = 1;
const int gamma = 1;
const int gammaR = 1;
const int kRS = 300000;
const int m = 1;
const int n = 3;
const int N = 4; //position of PV
const int dist = 3;//disturbance exponential

std::vector<std::vector<double>>
beforeD, beforeN, nextD, nextN, aveD, aveN, deltaD, deltaN, beforeR, deltaR, nextR, aveR, pow1, pow2;

std::vector<std::vector<int>>
firstD, firstN, firstR;

void calc_ave(int i, int j) {
	int ii1 = (i - 1 + x) % x;
	int ii2 = (i + 1) % x;
	int jj1 = (j - 1 + y) % y;
	int jj2 = (j + 1) % y;

	aveD[i][j] = (beforeD[i][jj2] + beforeD[ii2][j] + beforeD[i][jj1] + beforeD[ii1][j]) / 4;
	aveN[i][j] = (beforeN[i][jj2] + beforeN[ii2][j] + beforeN[i][jj1] + beforeN[ii1][j]) / 4;
	aveR[i][j] = (beforeR[i][jj2] + beforeR[ii2][j] + beforeR[i][jj1] + beforeR[ii1][j]) / 4;
}

int main() {
	FILE* fp;

	fopen_s(&fp, "sample.txt", "w");

	if (fp == NULL) {
		printf("failed\n");
		return -1;
	}
	else {
		printf("creating file, please wait\n");
	}

	beforeD.resize(x);
	beforeN.resize(x);
	nextD.resize(x);
	nextN.resize(x);
	aveD.resize(x);
	aveN.resize(x);
	deltaD.resize(x);
	deltaN.resize(x);
	beforeR.resize(x);
	deltaR.resize(x);
	nextR.resize(x);
	aveR.resize(x);
	pow1.resize(x);
	pow2.resize(x);
	firstD.resize(x);
	firstN.resize(x);
	firstR.resize(x);
	for (i = 0; i < x; i++) {
		beforeD[i].resize(y);
		beforeN[i].resize(y);
		nextD[i].resize(y);
		nextN[i].resize(y);
		aveD[i].resize(y);
		aveN[i].resize(y);
		deltaD[i].resize(y);
		deltaN[i].resize(y);
		beforeR[i].resize(y);
		deltaR[i].resize(y);
		nextR[i].resize(y);
		aveR[i].resize(y);
		pow1[i].resize(y);
		pow2[i].resize(y);
		firstD[i].resize(y);
		firstN[i].resize(y);
		firstR[i].resize(y);
	}

	/*
	std::random_device seed_gen;
	std::mt19937 engine(100);
	std::normal_distribution<> dist(1000, 100);
	for (i = 0; i < x; i++) {
		for (j = 0; j < y; j++) {
			beforeD[i][j] = dist(engine);
			beforeN[i][j] = dist(engine);
			beforeR[i][j] = 0.0;
		}
	}
	*/

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

	int numIter = 0;
	while (1)
	{
		int sumfirstD = 0, sumfirstN = 0, sumfirstR = 0;
		for (i = 0; i < x; i++) {
			for (j = 0; j < y; j++) {
				double PV = pow(10, dist);
				beforeD[N][N] = 1 * PV;

				calc_ave(i, j);

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

		numIter++;
		if (numIter % 10000 == 0 || numIter == 1)
			fprintf(stdout, "Iter%10d deltaD=%10.2e deltaN=%10.2e deltaR=%10.2e\n", numIter, deltaD[N + 3][N + 3], deltaN[N + 3][N + 3], deltaR[N + 3][N + 3]);

		if (sumfirstD == 399 && sumfirstN == 399 && sumfirstR == 399) {
			fprintf(stdout, "Iter%10d deltaD=%10.2e deltaN=%10.2e deltaR=%10.2e\n", numIter, deltaD[N + 3][N + 3], deltaN[N + 3][N + 3], deltaR[N + 3][N + 3]);

			break;
		}
	}

	double ave1 = 0;
	for (i = 0; i < x; i++) {
		for (j = 0; j < y; j++) {
			ave1 += beforeR[i][j];
		}
	}

	double ave2 = (ave1 - beforeR[N][N]) / 399;

	double dev1 = 0;
	for (i = 0; i < x; i++) {
		for (j = 0; j < y; j++) {
			dev1 += (beforeR[i][j] - ave2) * (beforeR[i][j] - ave2);
		}
	}

	double dev2 = dev1 - (beforeR[N][N] - ave2) * (beforeR[N][N] - ave2);

	double sd = sqrt(dev2 / 399);

	double diff = (beforeR[N + 1][N] - ave2) / sd;

	fprintf(fp, "Condition = sample\n");
	fprintf(fp, "\n");

	fprintf(fp, "Iteration Number= %8d\n", numIter);
	fprintf(fp, "\n");

	fprintf(fp, "diff (SD) = %.2f\n", diff);
	fprintf(fp, "\n");

	fprintf(fp, "Delta\n");

	for (i = 0; i < x; i++) {
		for (j = 0; j < y; j++) {
			fprintf(fp, " %.2f", beforeD[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fprintf(fp, "Notch\n");

	for (i = 0; i < x; i++) {
		for (j = 0; j < y; j++) {
			fprintf(fp, " %.2f", beforeN[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fprintf(fp, "Reporter\n");

	for (i = 0; i < x; i++) {
		for (j = 0; j < y; j++) {
			fprintf(fp, " %.2f", beforeR[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fclose(fp);
	return 0;
}