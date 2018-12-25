#ifndef ADF_H
#define ADF_H

double ADF(mat R, int n, double delta, vec moments);
double adfCov(int i, int j, int k, int h, mat R, vec moments);
int findpos (int i, int j, int k, int h);
double FRHO (int i, int j, int k, int h, vec M);

#endif
