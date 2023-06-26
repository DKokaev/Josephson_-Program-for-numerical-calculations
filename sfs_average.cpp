#include<iostream>
#include<cmath>
#include<cstdio>

//######################################################################################
// declearing the function of dmx/dt
double fmx(	double mx, double my, double mz,
			double ph,
			double G, double r, double alpha, double omega_F)
{
  double Heffx ,Heffy, Heffz;

  Heffx =0;
  Heffy= G*r*sin(ph-r*my);
  Heffz= mz;

  double result;
    result=-(omega_F/(1+pow(alpha, 2)))*(
				my*Heffz-mz*Heffy
                +alpha*(mx*(mx*Heffx+my*Heffy+mz*Heffz)-Heffx));
  return  result;
}
// declearing the function of dmy/dt
double fmy(	double mx, double my, double mz,
			double ph,
			double G, double r, double alpha, double omega_F)
{
  double Heffx ,Heffy, Heffz;

  Heffx =0;
  Heffy= G*r*sin(ph-r*my);
  Heffz= mz;

  double result;
    result=-(omega_F/(1+pow(alpha, 2)))*(
				mz*Heffx-mx*Heffz
                +alpha*(my*(mx*Heffx+my*Heffy+mz*Heffz)-Heffy));
  return  result;
}
// declearing the function of dmz/dt
double fmz(	double mx, double my, double mz,
			double ph,
			double G, double r, double alpha, double omega_F)
{
  double Heffx ,Heffy, Heffz;

  Heffx =0;
  Heffy= G*r*sin(ph-r*my);
  Heffz= mz;

  double result;
    result=-(omega_F/(1+pow(alpha, 2)))*(
				mx*Heffy-my*Heffx
                +alpha*(mz*(mx*Heffx+my*Heffy+mz*Heffz)-Heffz));
  return  result;
}
// declearing the function of dV/dt
double fV (	double V,	double ph,
			double mx, 	double my, 	double mz,
			double I, 	double beta_c,
			double G, 	double r, 	double alpha, double omega_F)
{
  double result;
    result=(1/beta_c)*(I
                   -(V-r*fmy(mx, my, mz, ph, G, r, alpha, omega_F))
                   -sin(ph-r*my)
                );
  return  result;
}
// declearing the function of dph/dt
double fph (double  V)
{
  double result;
    result=V;
  return  result;
}


int main()
	{
	using namespace std;
	int
		t,	//index of time
			//индекс по времени
		ti,	//number of time cycles from which of we begin integration
			//номер цикла по времени с которого начнается интегрирование
		a=1,//parameter for direction of current
			//параметр для определения направлении тока
		tN,	//total number of time cycles
			//количество точек (циклов) по времени

		breakcurrant=0,
		backcurrent=0;
		//######################################################################################
	double
		Ti=50,	//Initial value of time averaging
				//Начальная значения времени усреднения
		Tmax=1000,	//Time domain for one step of current
				//Временной домен для одного шага тока
		delta_t=0.01,//Step of time
				//Шаг по времени
		delta_I=0.005,	//Step of current outside current interval [Imore1,Imore2]
                        //Шаг по току вне интервала [Imore1,Imore2]

		I0=0.1,	//Initial value of bias current
				//Начальное значения тока
		Imax=1.2,	//Maximal value of bias current
					//Максимальное значения тока
		Imin=0,	//Минимальное значения тока
		Ibreak=0,//Ток при которой вычисление останавливается

		beta_c=25,	//McCumber parameter
					//Параметр МакКамбера
		r=0.5,	//Spin orbit coupling parameter
				//Параметр спин-орбитальной связи
		G=0.1,	//Отношение магнитной энергии к джозефсоновской
		alpha=0.1,//Гильбертовское затухание
		omega_F=0.5,//Частота ферромагнитного резонанса
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Runge-Kutta coefficients
		//Коэффициенты Рунге-Кутта
        mx, mx1, mx2, mx3, mx4,
        my, my1, my2, my3, my4,
        mz, mz1, mz2, mz3, mz4,
		ph, P1,P2,P3,P4,
		V, V1,V2,V3,V4,

        mx_max,    mx_min,
        my_max,    my_min,
        mz_max,    mz_min,
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		intV,//Интеграл напряжении по времени
		V_av,//Среднее напряжение
		I,	//Bias current
			//Входной ток
		time,	//время
		Is,		intIs,		Is_av,
		Iqp,	intIqp,		Iqp_av,
		Idisp,	intIdisp,	Idisp_av,
        V0;
//#######################################################################################
	//Creation result's files
	//Создания файлов для результатов
	FILE
	*f,
	*fIs,	*fIqp,	*fIdisp,

	*fmxmaxmin_I_up,	*fmymaxmin_I_up,	*fmzmaxmin_I_up,
	*fmxmax_V_up,		*fmymax_V_up,		*fmzmax_V_up,

	*fmxmaxmin_I_down,	*fmymaxmin_I_down,	*fmzmaxmin_I_down,
	*fmxmax_V_down,		*fmymax_V_down,		*fmzmax_V_down,

	*fVmaxmin_I_up,
	*fVmaxmin_I_down;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	f=fopen("Voltage.dat","w");

	fmxmaxmin_I_up=fopen("mx_maxmin_I_up.dat","w");
	fmymaxmin_I_up=fopen("my_maxmin_I_up.dat","w");
	fmzmaxmin_I_up=fopen("mz_maxmin_I_up.dat","w");

	fmxmax_V_up=fopen("mx_max_V_up.dat","w");
	fmymax_V_up=fopen("my_max_V_up.dat","w");
	fmzmax_V_up=fopen("mz_max_V_up.dat","w");


	fmxmaxmin_I_down=fopen("mx_maxmin_I_down.dat","w");
	fmymaxmin_I_down=fopen("my_maxmin_I_down.dat","w");
	fmzmaxmin_I_down=fopen("mz_maxmin_I_down.dat","w");

	fmxmax_V_down=fopen("mx_max_V_down.dat","w");
	fmymax_V_down=fopen("my_max_V_down.dat","w");
	fmzmax_V_down=fopen("mz_max_V_down.dat","w");

	fIs=fopen("Is_I.dat","w");
	fIqp=fopen("Iqp_I.dat","w");
	fIdisp=fopen("Idisp_I.dat","w");


	tN=(Tmax/delta_t);
	ti=tN*Ti/Tmax;
	//***************************************
	//initial condition
	//Начальные условия
	ph=0;
	intV=0;

	intIs=0;
	intIqp=0;
	intIdisp=0;

	V=0;

	mx=0;
	my=0;
	mz=1;

	I=I0;
//#######################################################################
	do{	// Here the cycle of current begins
		//Здесь начинается цикл по току
		mx_max=mx;	mx_min=mx;
		my_max=my;	my_min=my;
		mz_max=mz;	mz_min=mz;

		t=0;
		do{
			V0=V;

            time=t*delta_t;

            V1=delta_t*fV(   V,
                        ph,
                        mx,my, mz,
                        I, beta_c,
                        G, r, alpha, omega_F);
            P1=delta_t*fph(  V);
            mx1=delta_t*fmx( mx, my, mz,
                        ph,
                        G, r, alpha, omega_F
                        );
            my1=delta_t*fmy( mx, my, mz,
                        ph,
                        G, r, alpha, omega_F
                        );
            mz1=delta_t*fmz(   mx, my, mz,
                          ph,
                          G, r, alpha, omega_F
                          );
            //Вторые коэффициенты Рунге Кутта
            V2=delta_t*fV(   V+V1/2,
                        ph+P1/2,
                        mx+mx1/2,my+my1/2, mz+mz1/2,
                        I, beta_c,
                        G, r, alpha, omega_F);
            P2=delta_t*fph(  V+V1/2);
            mx2=delta_t*fmx( mx+mx1/2,my+my1/2, mz+mz1/2,
                        ph+P1/2,
                        G, r, alpha, omega_F
                        );
            my2=delta_t*fmy(   mx+mx1/2,my+my1/2, mz+mz1/2,
                          ph+P1/2,
                          G, r, alpha, omega_F
                          );
            mz2=delta_t*fmz( mx+mx1/2,my+my1/2, mz+mz1/2,
                        ph+P1/2,
                        G, r, alpha, omega_F
                        );
            //Третье коэффициенты Рунге-Кутта
            V3=delta_t*fV(
                        V+V2/2,
                        ph+P1/2,
                        mx+mx2/2,my+my2/2, mz+mz2/2,
                        I, beta_c,
                        G, r, alpha, omega_F);
            P3=delta_t*fph(V+V2/2);
            mx3=delta_t*fmx( mx+mx2/2,my+my2/2, mz+mz2/2,
                        ph+P2/2,
                        G, r, alpha, omega_F
                        );
            my3=delta_t*fmy( mx+mx2/2,my+my2/2, mz+mz2/2,
                        ph+P2/2,
                        G, r, alpha, omega_F
                        );
            mz3=delta_t*fmz( mx+mx2/2,my+my2/2, mz+mz2/2,
                        ph+P2/2,
                        G, r, alpha, omega_F
                    );
            //Четвертые коэффициенты Рунге-Кутта
            V4=delta_t*fV(   V+V3,
                        ph+P3,
                        mx+mx3,my+my3, mz+mz3,
                        I, beta_c,
                        G, r, alpha, omega_F);
            P4=delta_t*fph(V+V3);
            mx4=delta_t*fmx(   mx+mx3,my+my3, mz+mz3,
                          ph+P3,
                          G, r, alpha, omega_F
                          );
            my4=delta_t*fmy(   mx+mx3,my+my3, mz+mz3,
                          ph+P3,
                          G, r, alpha, omega_F
						  );
            mz4=delta_t*fmz(   mx+mx3,my+my3, mz+mz3,
                          ph+P3,
                          G, r, alpha, omega_F
						  );
			V=V+(1/6.)*(V1+2*V2+2*V3+V4);

			ph=ph+(1/6.)*(P1+2*P2+2*P3+P4);

			mx=mx+(1/6.)*(mx1+2*mx2+2*mx3+mx4);
			my=my+(1/6.)*(my1+2*my2+2*my3+my4);
			mz=mz+(1/6.)*(mz1+2*mz2+2*mz3+mz4);


			Is=sin(ph-r*my);
			Iqp=V;
			Idisp=beta_c*(V-V0)/delta_t;

            if (t>ti){
                intV=V*delta_t+intV;

				intIs=Is*delta_t+intIs;
				intIqp=Iqp*delta_t+intIqp;
				intIdisp=Idisp*delta_t+intIdisp;
            };
	//Ending of Runge-Kutta method
	//Конец метода Рунге-Кутта
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            if(mx>mx_max){mx_max=mx;};
			if(mx<mx_min){mx_min=mx;};
            if(my>my_max){my_max=my;};
			if(my<my_min){my_min=my;};
			if(mz>mz_max){mz_max=mz;};
			if(mz<mz_min){mz_min=mz;};
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//Here the cycle of time finishs
			//Здесь цикл по времени закончивается
		}while(++t<tN+1);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Voltage averaging
		//Усреднение напряжения
		V_av=intV/(Tmax-Ti);
		intV=0;

		Is_av=intIs/(Tmax-Ti);
		intIs=0;
		Iqp_av=intIqp/(Tmax-Ti);
		intIqp=0;
		Idisp_av=intIdisp/(Tmax-Ti);
		intIdisp=0;

		fprintf(fIs,"%f\t%e\n",I,Is_av);
		fprintf(fIqp,"%f\t%e\n",I,Iqp_av);
		fprintf(fIdisp,"%f\t%e\n",I,Idisp_av);

		if(a==1){
			fprintf(fmxmaxmin_I_up,"%f\t%e\t%e\n",I,mx_max,mx_min);
			fprintf(fmymaxmin_I_up,"%f\t%e\t%e\n",I,my_max,my_min);
			fprintf(fmzmaxmin_I_up,"%f\t%e\t%e\n",I,mz_max,mz_min);

			fprintf(fmxmax_V_up,"%f\t%e\n",V_av,mx_max);
			fprintf(fmymax_V_up,"%f\t%e\n",V_av,my_max);
			fprintf(fmzmax_V_up,"%f\t%e\n",V_av,mz_max);
		}

		if(a==-1){
			fprintf(fmxmaxmin_I_down,"%f\t%e\t%e\n",I,mx_max,mx_min);
			fprintf(fmymaxmin_I_down,"%f\t%e\t%e\n",I,my_max,my_min);
			fprintf(fmzmaxmin_I_down,"%f\t%e\t%e\n",I,mz_max,mz_min);

			fprintf(fmxmax_V_down,"%f\t%e\n",V_av,mx_max);
			fprintf(fmymax_V_down,"%f\t%e\n",V_av,my_max);
			fprintf(fmzmax_V_down,"%f\t%e\n",V_av,mz_max);
		}
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//showing results on the screen
		//вывод результатов на экран
		cout<<"\n\nI:"<<I<<"\tV:"<<V_av;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Here will be recorded bias Current and total voltage
		//Здесь записывается входной ток и суммарное напряжение
		fprintf(f,"%f\t%e\n",I,V_av);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Условие для возврата тока
		if ((I<=Imin+delta_I/2)&&(backcurrent==0)&&(a==-1)){
			backcurrent=1;
			breakcurrant=1;
			a=1;
		};
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Condition for changing current dirrection
		//Условие для изменения направлении тока
		if (I>Imax){
			cout<<"here is max of current so we  go back"<<endl;
			a=-1;};
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Here bias current will be changed
		//Здесь изменяется направление тока
		I=I+a*delta_I;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if((backcurrent==1)&&(I>=Ibreak-delta_I/2)&&(a==1))
			{
				break;
			};
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Here cycle of current finishs
	//Здесь закончивается цикл по току
	}while(I<=100);
	//**********************
	// Here all files will be closed
	//Здесь все файлы закрываются
	fclose(f);
    fclose(fmxmaxmin_I_up);
	fclose(fmymaxmin_I_up);
	fclose(fmzmaxmin_I_up);

	fclose(fmxmax_V_up);
	fclose(fmymax_V_up);
	fclose(fmzmax_V_up);


	fclose(fmxmaxmin_I_down);
	fclose(fmymaxmin_I_down);
	fclose(fmzmaxmin_I_down);

	fclose(fmxmax_V_down);
	fclose(fmymax_V_down);
	fclose(fmzmax_V_down);

	fclose(fIs);
	fclose(fIqp);
	fclose(fIdisp);
return(0);
}
