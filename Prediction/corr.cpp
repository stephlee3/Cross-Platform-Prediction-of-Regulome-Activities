#include<iostream>
#include <math.h>
using namespace std;
double exon[1200000][15];
double micro[1200000][15];
double dnase[1200000][15];
double corr[1200000];

double a[10];
double b[10];
int n=1108603;
int k=7;

double cor(double x[],double y[],int m){
	
	double x_sum=0,y_sum=0;
	double x_mean=0,y_mean=0;
	for(int i=1;i<=m;i++){
		x_sum+=x[i];
		y_sum+=y[i];
	}
	x_mean= x_sum/m;
	y_mean= y_sum/m;
	double c=0,d=0,e=0;
	for(int i=1;i<=m;i++){
		c+=(x[i]-x_mean)*(y[i]-y_mean);
	}
	for(int i=1;i<=m;i++){
		d+=(x[i]-x_mean)*(x[i]-x_mean);
	}
	for(int i=1;i<=m;i++){
		e+=(y[i]-y_mean)*(y[i]-y_mean);
	}
	double r=0;
	if(d==0||e==0) return 0;
	r=c/(sqrt(d)*sqrt(e)); 
	return r;
	
}
double pred_err(double x[1200000][15],double y[1200000][15],int l,int m){
	double ans=0.0;
	double a=0,b=0;
	double y_mean=0,y_sum=0;
	for(int i=1;i<=l;i++){
		for(int j=1;j<=m;j++){
			y_sum+=y[i][j];
		}
	}
	y_mean=y_sum/(l*m);
	
	for(int i=1;i<=l;i++){
		for(int j=1;j<=m;j++){
			a+=(x[i][j]-y[i][j])*(x[i][j]-y[i][j]);
			b+=(y[i][j]-y_mean)*(y[i][j]-y_mean);
		}
	}
	ans=a/b;
	return ans;
}

int main(){
	
	freopen("DH_pred_exon.txt","r",stdin);
	for(int i=1;i<=n;i++){
		for(int j=1;j<=k;j++){
			cin>>exon[i][j];
		}
	}
	freopen("DH_pred_micro.txt","r",stdin);
	for(int i=1;i<=n;i++){
		for(int j=1;j<=k;j++){
			cin>>micro[i][j];
		}
	}
	freopen("DNase_c.txt","r",stdin);
	for(int i=1;i<=n;i++){
		for(int j=1;j<=k;j++){
			cin>>dnase[i][j];
		}
	}
	
	freopen("cor_micro.txt","w",stdout);
	for(int i=1;i<=n;i++){
		corr[i]=cor(micro[i],dnase[i],k);
		cout<<corr[i]<<endl;
	}
	freopen("cor_exon.txt","w",stdout);
	for(int i=1;i<=n;i++){
		corr[i]=cor(exon[i],dnase[i],k);
		cout<<corr[i]<<endl;
	}
	
	double t1,t2=0;
	t1=pred_err(exon,dnase,n,k);
	t2=pred_err(micro,dnase,n,k);
	freopen("pred_err.txt","w",stdout);
	cout<<"exon:"<<t1<<endl;
	cout<<"micro:"<<t2<<endl; 
	
  return 0;
} 
