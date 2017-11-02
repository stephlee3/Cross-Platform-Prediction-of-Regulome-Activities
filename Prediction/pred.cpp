#include <iostream>
#define MAXN 1200000
using namespace std;

double coef_all[MAXN][20];
int predictor_idx[MAXN][20];
double micro_test_mean[2000][20];
double exon_test_mean[2000][20];
double DH_exon[MAXN][20];
double DH_micro[MAXN][20];
double DNase[MAXN][5];
//coef_all:1108603*7
//predictor_idx:1108603*7
//micro_test_mean:1000*10
//exon_test_mean:1000*10
	
//DH.pred.exon 1108603*10
//DH.pred.micro 1108603*10
//DNase_mean_sd 1108603*2
int n=1108603;
int m=500;
int q=8;
int ct=7;

int main(){
	// read the data
	freopen("regress_coef.txt","r",stdin);
	for(int i=1;i<=n;i++){
		for(int j=1;j<=q;j++){
			cin>>coef_all[i][j];
		}
	}
	fclose(stdin); 
	
	freopen("regress_predictor.txt","r",stdin);
	for(int i=1;i<=n;i++){
		for(int j=1;j<=q;j++){
			cin>>predictor_idx[i][j];
		}
	}
	fclose(stdin);
	
	freopen("exon_test_mean.txt","r",stdin);
	for(int i=1;i<=m;i++){
		for(int j=1;j<=ct;j++){
			cin>>exon_test_mean[i][j];
		}
	}
	fclose(stdin);
	
	freopen("micro_test_mean.txt","r",stdin);
	for(int i=1;i<=m;i++){
		for(int j=1;j<=ct;j++){
			cin>>micro_test_mean[i][j];
		}
	}
	fclose(stdin);
/*	
	for(int i=1;i<=n;i++){
		for(int j=1;j<=7;j++){
			cout<<coef_all[i][j]<<" "<<predictor_idx[i][j];
			cout<<endl; 
		}
	}
	
	for(int i=1;i<=m;i++){
		for(int j=1;j<=10;j++){
			cout<<exon_test_mean[i][j]<<" "<<micro_test_mean[i][j]; 
			cout<<endl;
		}
	}
*/	
	// prediction
	for(int i=1;i<=n;i++){
		double coef_tmp[20];
		int  pred_idx[20];
		for(int j=1;j<=q;j++){
			coef_tmp[j]=coef_all[i][j];
			pred_idx[j]=predictor_idx[i][j];
		}
		double expr_val[20];
		double expr_val2[20];
		for(int j=1;j<=ct;j++){
			for(int k=1;k<=q;k++){
				expr_val[k]=micro_test_mean[pred_idx[k]][j];
				expr_val2[k]=exon_test_mean[pred_idx[k]][j];
			}
			for(int k=1;k<=q;k++){
				DH_micro[i][j]+=expr_val[k]*coef_tmp[k];
				DH_exon[i][j]+=expr_val2[k]*coef_tmp[k];
			}
		}
	}
	//scale
	//double DNase[MAXN][5];
	freopen("DNase_mean_sd.txt","r",stdin);
	for(int i=1;i<=n;i++){
		cin>>DNase[i][1]>>DNase[i][2];
	}
	fclose(stdin);
	
	for(int i=1;i<=n;i++){
		for(int j=1;j<=ct;j++){
			DH_micro[i][j]=DH_micro[i][j]*DNase[i][2]+DNase[i][1];
			DH_exon[i][j]=DH_exon[i][j]*DNase[i][2]+DNase[i][1];
		}
	}
	
	freopen("DH_pred_micro.txt","w",stdout);
	for(int i=1;i<=n;i++){
		for(int j=1;j<=ct;j++){
			cout<<DH_micro[i][j]<<" ";
		}
		cout<<endl;
	}
	fclose(stdout);
//	cout<<"*************";
	
	freopen("DH_pred_exon.txt","w",stdout);
	for(int i=1;i<=n;i++){
		for(int j=1;j<=ct;j++){
			cout<<DH_exon[i][j]<<" ";
		}
		cout<<endl;
	}
	fclose(stdout);
	
	return 0;
}
