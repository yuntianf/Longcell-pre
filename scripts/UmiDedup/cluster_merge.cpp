#include <Rcpp.h>
#include<string>
#include <algorithm>
#include<vector>
#include<set>
#include<math.h>
#include<cmath>

using namespace std;
using namespace Rcpp;


int cluster_disance(NumericMatrix cluster, NumericMatrix dis,int cluster1, int cluster2);

// [[Rcpp::export]]

NumericMatrix cluster_merge(NumericMatrix cluster, NumericMatrix dis, int thresh = 2){
  NumericVector umi_id = cluster(_,0);
  NumericVector cluster_id = cluster(_,1);
    
  NumericVector cluster_type = unique(cluster_id).sort(); 
  int cluster_num = cluster_type.size();
  if(cluster_num == 1){
    return(cluster);
    }
  
  NumericMatrix cluster_dis(cluster_num,cluster_num);
  cluster_dis.fill(-1);
  
  for(int i = 0;i < cluster_num;i++){
    for(int j = i+1; j < cluster_num;j++){
      int temp = cluster_disance(cluster,dis,cluster_type[i],cluster_type[j]);
      if(temp >= 0 & temp <= thresh){
        cluster_id[cluster_id == cluster_type[j]] = cluster_type[i];
        NumericMatrix update_cluster(cluster_id.size(),2);
        update_cluster.column(0) = umi_id;
        update_cluster.column(1) = cluster_id;
        return(cluster_merge(update_cluster,dis,thresh));
        }
      }
    }
  
  return(cluster);
  }

int cluster_disance(NumericMatrix cluster, NumericMatrix dis,int cluster1, int cluster2){
  NumericVector umi_id = cluster(_,0);
  NumericVector cluster_id = cluster(_,1);
  
  NumericVector cluster1_id = umi_id[cluster_id == cluster1];
  NumericVector cluster2_id = umi_id[cluster_id == cluster2];

  
  double cluster1_size = cluster1_id.size();
  double cluster2_size = cluster2_id.size();
  double dis_sum = 0;
  
  double dis_size = cluster1_size*cluster2_size;
  LogicalVector loc = in(dis(_,0),cluster1_id) & in(dis(_,1),cluster2_id);
  NumericVector dis_sum_vec = dis(_,2);
  dis_sum_vec = dis_sum_vec[loc];
  dis_sum = sum(dis_sum_vec);
  
  /*
  for(int i = 0;i < cluster1_size;i++){
    for(int j = 0; j < cluster2_size;j++){
      int temp = dis(cluster1_id[i]-1,cluster2_id[j]-1);
      if(temp == -1){
        inf_num++;
        }
      else{
        dis_sum += temp;
        }
      }
  }*/
  
  if(dis_sum_vec.size() < 4*dis_size/5){
    return(-1);
    }
  else{
    return(dis_sum/dis_sum_vec.size());
    }
  }