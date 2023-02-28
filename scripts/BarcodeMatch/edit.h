#pragma once
#include<string>
#include <algorithm>
#include<vector>
#include<set>
#include<math.h>
#include<cmath>
#include<fstream>
#include <stdlib.h>
#include <utility>

using namespace std;

vector<string> reverse(vector<string> s);
set<string> kmer(vector<string> s, int k);
vector<string> kmer(string s, int k);
int editDist(string word1, string word2);
pair<int, int> minEditDist(string seq, string barcode);