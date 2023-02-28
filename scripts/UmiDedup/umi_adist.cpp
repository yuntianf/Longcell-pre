#include <Rcpp.h>
#include<string>
#include <algorithm>
#include<vector>
#include<set>
#include<math.h>
#include<cmath>

using namespace std;
using namespace Rcpp;
vector<string> kmer(string s, int k);
int editDist(string word1, string word2);

int minEditDis(string seq1, string seq2, int k);
bool kmer_overlap(string seq1, string seq2, int k);

// [[Rcpp::export]]
NumericMatrix umi_adist(vector<string> umi_flank, int umi_len = 10, int k = 6){
  int umi_size = umi_flank.size();
  NumericVector umi1,umi2,umi_dis;
  int temp_dis = 0;
  
  for(int i =0;i < umi_size;i++){
    for(int j = i+1;j < umi_size;j++){
      if(kmer_overlap(umi_flank[i],umi_flank[j],k)){
        temp_dis = minEditDis(umi_flank[i],umi_flank[j],umi_len);
        umi1.push_back(i+1);
        umi2.push_back(j+1);
        umi_dis.push_back(temp_dis);
      }
    }
  }
  
  NumericMatrix dis(umi1.size(),3);
  dis(_,0) = umi1;
  dis(_,1) = umi2;
  dis(_,2) = umi_dis;
  return(dis);
}

vector<string> kmer(string s, int k) {
  int s_len = s.length();
  vector<string> out;
  
  for (int i = 0; i < s_len - k + 1; i++) {
    out.push_back(s.substr(i, k));
  }
  return(out);
}


int editDist(string word1, string word2) {
  vector<vector<int>> dp = vector<vector<int>>(word1.size() + 1, vector<int>(word2.size() + 1, 0));
  
  for (int i = 0; i <= word1.size(); i++) {
    dp[i][0] = i;
  }
  
  for (int j = 1; j <= word2.size(); j++) {
    dp[0][j] = j;
  }
  
  for (int i = 0; i < word1.size(); i++) {
    for (int j = 0; j < word2.size(); j++) {
      if (word1[i] == word2[j]) {
        dp[i+1][j+1] = dp[i][j];
      }
      else {
        dp[i+1][j+1] = min(min(dp[i][j+1], dp[i+1][j]), dp[i][j]) + 1;
      }
    }
  }
  
  return dp[word1.size()][word2.size()];
}

int minEditDis(string seq1, string seq2, int k = 10) {
  vector<string> seq1_kmer = kmer(seq1, k);
  vector<string> seq2_kmer = kmer(seq2, k);
  
  int kmer_size = seq1_kmer.size();
  
  int edit = seq1.size();
  for(int i = 0;i < kmer_size;i++){
    for(int j = 0;j < kmer_size;j++){
      int dis = editDist(seq1_kmer[i],seq2_kmer[j]);
      if(edit > dis){
        edit = dis;
        }
      if(dis == 0){
        break;
        }
      }
    }
  return(edit);
}

bool kmer_overlap(string seq1, string seq2, int k = 6){
  vector<string> seq1_kmer = kmer(seq1, k);
  vector<string> seq2_kmer = kmer(seq2, k);
  
  int kmer_size = seq1_kmer.size();
  
  for(int i = 0;i < kmer_size;i++){
    for(int j = 0;j < kmer_size;j++){
      if(seq1_kmer[i] == seq2_kmer[j]){
        return(true);
        }
      }
  }
  return(false);
}

#/*** R
#editDist("TGCGCTTGATAATT","TGTGTGATAATTTC")
#*/
