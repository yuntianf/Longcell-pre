#include <Rcpp.h>
#include<string>
#include <algorithm>
#include<vector>
#include<set>
#include<math.h>
#include<cmath>
#include <stack>

using namespace std;
using namespace Rcpp;

int UMI_NS[10000][10000] = {};
int UMI_CORRES[50000] = {};
static int none = -100;

// ##### isoform cluster #####
vector<string> flatten(vector<vector<string> > const &vec){
  vector<string> flattened;
  for (auto const &v: vec) {
    flattened.insert(flattened.end(), v.begin(), v.end());
  }
  return flattened;
}


inline vector<string> str_split(string s, string split){
  vector<string> out;
  
  size_t pos = 0;
  string token;
  while ((pos = s.find(split)) != string::npos) {
    token = s.substr(0, pos);
    out.push_back(token);
    s.erase(0, pos + split.length());
  }
  if(s.size() > 0){
    out.push_back(s);
  }
  
  return(out);
}


inline vector<int> isoform2sites(string iso, string split = "|", string sep = ","){
  vector<string> exons = str_split(iso,split);
  
  vector<vector<string> > bins;
  int len = exons.size();
  for(int i = 0;i < len;i++){
    bins.push_back(str_split(exons[i],sep));
  }
  vector<string> sites_chr = flatten(bins);
  
  vector<int> sites_int;
  transform(sites_chr.begin(), sites_chr.end(), back_inserter(sites_int),
            [](const string& str) { return stoi(str); });
  
  return(sites_int);
}


inline int iso_len(vector<int> sites){
  int sites_size = sites.size();
  
  if(sites_size % 2 > 0){
    cout << "the length of sites should be an even" << endl;
    exit(1);
  }
  int len = 0;
  for(int i = 0;i < sites_size/2;i++){
    int sites_s = sites[i*2];
    int sites_e = sites[i*2 + 1];
    
    len= len+sites_e-sites_s+1;
  } 
  return(len);
}


inline int bin2_intersect(int a_start,int a_end, 
                          int b_start,int b_end){
  int len = 0;
  if(a_start >= b_start && a_start <= b_end){
    len = min({a_end,b_end}) - a_start  + 1;
  }
  else if(b_start >= a_start && b_start <= a_end){
    len = min({a_end,b_end}) - b_start  + 1;
  }
  else{
    len = 0;
  }
  return(len);
}


inline vector<int> sites_chop(vector<int> sites, int start,int end){
  int sites_size = sites.size();
  
  if(sites_size % 2 > 0){
    cout << "the length of sites should be an even" << endl;
    exit(1);
  }
  
  vector<int> chop;
  for(int i = 0;i < sites_size/2;i++){
    int sites_s = sites[i*2];
    int sites_e = sites[i*2 + 1];
    if(sites_e <= start || sites_s >= end){
      continue;
    }
    else{
      chop.push_back(max({sites_s,start}));
      chop.push_back(min({sites_e,end}));
    }
  }
  return(chop);
}


inline int iso2_dis(string a,string b,
                        string split = "|",string sep = ","){
  if(a == b){
    return(0);
  }
  vector<int> sites_a = isoform2sites(a, split,sep);
  vector<int> sites_b = isoform2sites(b, split,sep);
  
  int a_size = sites_a.size();
  int b_size = sites_b.size();
  
  int start = max({sites_a[0],sites_b[0]});
  int end = min({sites_a[a_size - 1],sites_b[b_size - 1]});
  
  if(start >= end){
    return(-1);
  }
  
  vector<int> chop_a = sites_chop(sites_a, start, end);
  vector<int> chop_b = sites_chop(sites_b, start, end);
  
  int chop_a_size = chop_a.size();
  int chop_b_size = chop_b.size();
  
  int intersect = 0;
  int next_start = 0;
  for(int i =0;i < chop_a_size/2;i++){
    int j = next_start;
    bool flag = false;
    
    while(j < chop_b_size/2){
      int a_start = chop_a[2*i];
      int a_end = chop_a[2*i + 1];
      int b_start = chop_b[2*j];
      int b_end = chop_b[2*j + 1];
      
      int inter = bin2_intersect(a_start,a_end,b_start,b_end);
      
      if(inter > 0){
        if(!flag){
          next_start = j;
          flag = true;
        }
        intersect += inter;
      }
      else{
        if(flag){
          break;
        }
      }
      j++;
    }
  }
  
  int dis = iso_len(chop_a) + iso_len(chop_b) - 2*intersect;
  return(dis);
}

NumericMatrix iso_dist(vector<string> iso_set,int thresh = 80,
                       string split = "|",string sep = ","){
  int iso_size = iso_set.size();
  
  NumericMatrix iso_dis ((iso_size+1)*iso_size/2,2);
  int id = 0;
  
  for(int i = 0;i < iso_size;i++){
    for(int j = i;j < iso_size;j++){
      if(i == j){
        iso_dis(id,0) = i;
        iso_dis(id,1) = j;
        id = id + 1;
      }
      else{
        int dis = iso2_dis(iso_set[i],iso_set[j],
                           split,sep);
        //cout << i << "-" << j<<":" << dis << endl;
        if(dis >=0 && dis <= thresh){
          iso_dis(id,0) = i;
          iso_dis(id,1) = j;    
          id = id+ 1;
        }
      }
    }
  }
  return(iso_dis(Range(0,id-1),_));
}

// ##### UMI distance #####
// ### needleman ###
inline int needle(string A,string B,
                        int match_score = 1, 
                        int mismatch_score = -1,
                        int gap_score = -1){
  int n = A.size();
  int m = B.size();
  
  NumericMatrix dp(n+1,m+1);
  for (int i=0;i<=n;i++) dp(i,0) = dp(0,i) = i * gap_score;
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=m;j++)
    {
      int S = (A[i-1] == B[j-1]) ? match_score : mismatch_score;
      dp(i,j) = max({dp(i-1,j-1) + S, max({dp(i-1,j) + gap_score, dp(i,j-1) + gap_score})});
    }
  }
  
  return dp(n,m);
}

// ### edit distance ###
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


// ### index ###
inline int pair2id(int x, int y){
  int m = max({x,y});
  int n = min({x,y});
  int loc = (m+1)*m/2;
  int id = loc+n;
  return(id);
}

inline List index(vector<string> data, vector<string> uniq){
  List L;
  int uniq_size = uniq.size();
  int data_size = data.size();
  
  for(int i =0;i < uniq_size;i++){
    vector<int> sub_index;
    for(int j = 0;j <data_size;j++){
      if(data[j] == uniq[i]){
        sub_index.push_back(j);
      }
    }
    L.push_back(sub_index);
  }
  
  return(L);
}

inline vector<string> vec_extract(vector<string> data, vector<int> index){
  int len = index.size();
  
  vector<string> extract;
  for(int i = 0;i < len;i++){
    extract.push_back(data[index[i]]);
  }
  
  return(extract);
}


vector<vector<int> > umi_needle(vector<string> umi, vector<int> count,
                                vector<int> umi_id1, vector<int> umi_id2,
                                int thresh = 5) {
  int umi_set1_size = umi_id1.size();
  int umi_set2_size = umi_id2.size();
  
  vector<vector<int> > umi_ns;
  
  for (int i = 0; i < umi_set1_size; i++) {
    for (int j = 0; j < umi_set2_size; j++) {
      int ns = none;
      if(UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]] != none){
        ns = UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]];
      }
      else {
        ns = needle(umi[umi_id1[i]], umi[umi_id2[j]]);
        UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]] = ns;
        UMI_NS[UMI_CORRES[umi_id2[j]]][UMI_CORRES[umi_id1[i]]] = ns;
      }
      if (ns >= thresh) {
        int score = 1;
        if(umi_id1[i] == umi_id2[j]){
          score = count[umi_id1[i]]*(count[umi_id1[i]]-1)/2;
        }
        else{
          score = count[umi_id1[i]]*count[umi_id2[j]];
        }
        vector<int> sub_ns = { umi_id1[i],umi_id2[j],score*ns,score };
        umi_ns.push_back(sub_ns);
      }
    }
  }
  return(umi_ns);
}

vector<vector<int> > umi_edit(vector<string> umi,
                                vector<int> umi_id1, vector<int> umi_id2,
                                int thresh = 2, int k = 10) {
  int umi_set1_size = umi_id1.size();
  int umi_set2_size = umi_id2.size();
  
  vector<vector<int> > umi_ns;
  
  for (int i = 0; i < umi_set1_size; i++) {
    for (int j = 0; j < umi_set2_size; j++) {
      int ns = none;
      if(UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]] != none){
        ns = UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]];
      }
      else {
        ns = minEditDis(umi[umi_id1[i]], umi[umi_id2[j]],k);
        UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]] = ns;
        UMI_NS[UMI_CORRES[umi_id2[j]]][UMI_CORRES[umi_id1[i]]] = ns;
      }
      if (ns <= thresh) {
        vector<int> sub_ns = { umi_id1[i],umi_id2[j],ns };
        umi_ns.push_back(sub_ns);
      }
    }
  }
  return(umi_ns);
}


inline void initialize(vector<string> umi){
  set<string> umi_set(umi.begin(), umi.end());
  vector<string> umi_uniq;
  umi_uniq.assign(umi_set.begin(), umi_set.end());
  
  int umi_size = umi.size();
  int umi_uniq_size = umi_uniq.size();
  int umi_len = umi[0].size();
  
  for(int i = 0;i < umi_uniq_size;i++){
    for(int j = 0;j < umi_uniq_size;j++){
      if(i == j){
        UMI_NS[i][j] = umi_len;
      }
      else{
        UMI_NS[i][j] = none;
      }
    }
  }
  
  for(int i = 0;i < umi_size;i++){
    for(int j = 0;j < umi_uniq_size;j++){
      if(umi[i] == umi_uniq[j]){
        UMI_CORRES[i] = j;
      }
    }
  }
}


// [[Rcpp::export]]
vector<vector<int> > umi_graph_table(vector<string> umi,
                              vector<string> isoform,
                              vector<int> count,
                              int sim_thresh = 5,int iso_thresh = 80,
                              string split = "|",string sep = ","){
  int umi_len = umi.size();
  int iso_len = isoform.size();
  int count_len = count.size();
  
  if(iso_len != umi_len || iso_len != count_len || umi_len != count_len){
    cout << "The size of isoforms and umi don't match!" << endl;
  }
  
  set<string> iso_set(isoform.begin(), isoform.end());
  vector<string> iso_uniq;
  iso_uniq.assign(iso_set.begin(), iso_set.end());
  
  NumericMatrix iso_src = iso_dist(iso_uniq,iso_thresh,
                                   split = split,sep = sep);
  //return(iso_src);
  List iso_index = index(isoform,iso_uniq);
  initialize(umi);
  
  int iso_pair_count = iso_src.nrow();
  vector<vector<int> > umi_ns;
  for(int i = 0;i < iso_pair_count;i++){
    vector<vector<int> > temp_ns = umi_needle(umi,count,
                                       iso_index[iso_src(i,0)],
                                       iso_index[iso_src(i,1)],
                                       sim_thresh);
    
    if(temp_ns.size() > 0){
      umi_ns.insert(umi_ns.end(),temp_ns.begin(),temp_ns.end());
    }
  }
  return(umi_ns);
}

// [[Rcpp::export]]
vector<vector<int> > umi_edit_table(vector<string> umi,
                                     vector<string> isoform,
                                     int edit_thresh = 2,
                                     int k = 10,
                                     int iso_thresh = 80,
                                     string split = "|",string sep = ","){
  int umi_len = umi.size();
  int iso_len = isoform.size();
  
  if(iso_len != umi_len){
    cout << "The size of isoforms and umi don't match!" << endl;
  }
  
  set<string> iso_set(isoform.begin(), isoform.end());
  vector<string> iso_uniq;
  iso_uniq.assign(iso_set.begin(), iso_set.end());
  
  NumericMatrix iso_src = iso_dist(iso_uniq,iso_thresh,
                                   split = split,sep = sep);
  //return(iso_src);
  List iso_index = index(isoform,iso_uniq);
  initialize(umi);
  
  int iso_pair_count = iso_src.nrow();
  vector<vector<int> > umi_ns;
  for(int i = 0;i < iso_pair_count;i++){
    vector<vector<int> > temp_ns = umi_edit(umi,
                                            iso_index[iso_src(i,0)],
                                            iso_index[iso_src(i,1)],
                                            edit_thresh,k);
    
    if(temp_ns.size() > 0){
      umi_ns.insert(umi_ns.end(),temp_ns.begin(),temp_ns.end());
    }
  }
  return(umi_ns);
}