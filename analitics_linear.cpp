//SFS with nonconstant Ic CVC 29-09-2021
//Overdamped JJ
#include<iostream>
#include<cmath>
#include<cstdio>
#include<stdlib.h>

int main()
	{
	using namespace std;

	double
        wF=0.5,
        G=0.1,
        r=0.1,
        alpha=0.1,
        delta_w=0.01,
        w_max=1,
        my_max,
        w;

	//Creation result's files
	//Создания файлов для результатов
	FILE
	*f;

	f=fopen("my_max_w.dat","w");


    w=0;
    do{
        my_max=pow(wF,2)*G*r/sqrt(pow((pow(wF,2)-pow(w,2)),2)+pow(2*alpha*wF*w,2));

		fprintf(f,"%f\t%.16f\n",w,my_max);
        w=w+delta_w;
        cout<<w<<endl;
    }while(w<=w_max);
	fclose(f);

return(0);
}
