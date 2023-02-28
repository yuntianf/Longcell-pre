#pragma once
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include<math.h>
#include<cmath>
#include <stdlib.h>
#include <utility>
#include <iostream>
#include <fstream>

using namespace std;

vector<int> kmer_include(string seq, set<string> dic);
vector<vector<int> > barcodes_cos_vector(vector<string> barcodes, set<string> dic);
double cos_sim(vector<int> a, vector<int> b);
bool cos_sim_comp(pair <int, double> cos1, pair <int, double> cos2);
vector<string> barcode_cand_cos(string seq, vector<string> barcodes, set<string> dic, vector<vector<int> > index,
    int top = 8, double thresh = 0.25);
vector<int> pos_filter(vector<double> start, vector<int> edit);
void barcodeMatch(vector<string> seq, vector<string> barcodes, vector<string> readnames,
    double mu = 20, double sigma = 10, double sigma_start = 10, int k = 8, int batch = 100,
    int top = 8, double cos_thresh = 0.25, double alpha = 0.05, int edit_thresh = 6,
    string outname = "integration_result.txt");
