//  ******************************************************************
//  LASLETT Model from Laslett, Eveson, and Polacheck (2002).
//  Created by Martell on 2012-05-22.
//  Copyright (c) 2012. All rights reserved.
//  Comments: 
//   -Acknowledge Paige Eveson for help with the bivariate normal
//  ******************************************************************


DATA_SECTION
	!! ad_comm::change_datafile_name("HBCdata.dat");
	int iyr;
	LOC_CALCS
		int on,opt;
		if((on = option_match(ad_comm::argc,ad_comm::argv,"-iyr",opt))>-1)
		{
			iyr = atoi(ad_comm::argv[on+1]);
		}
		else
		{
			iyr = 1989;
		}
	END_CALCS
	init_int nobs;
	int i1;
	int i2;
	//#YEAR RIVER L1 L2 DT DL
	init_matrix data(1,nobs,1,6);
	LOC_CALCS
		int i;
		i1=0;
		i2=0;
		for(i=1;i<=nobs;i++)
		{
			if(data(i,1)==iyr && i1==0) i1 = i;
			if(data(i,1)==iyr+1) break;
			i2 = i;
		}
	END_CALCS
	ivector ii(i1,i2);
	!! ii.fill_seqadd(i1,1);
	vector l1(i1,i2);
	vector l2(i1,i2);
	vector dt(i1,i2);
	LOC_CALCS
		for(i=i1;i<=i2;i++)
		{
			l1(i) = data(i,3);
			l2(i) = data(i,4);
			dt(i) = data(i,5)/365.25;
		}
	END_CALCS

INITIALIZATION_SECTION
 k 0.2;

PARAMETER_SECTION
	init_bounded_number mu_linf(50,550,1);
	init_bounded_number sd_linf(0,100,2);
	init_bounded_number k(0,10,1);
	init_bounded_number sig(0,150,1);
	init_bounded_number muA(0,50,1);
	init_bounded_number sigA(0,25,2);
	//init_bounded_dev_vector A(1,nobs,-5,5,2);
	random_effects_vector A(i1,i2,2);
	objective_function_value f;
	
	number sig2;
	sdreport_number sd2_linf;
	
	vector age(i1,i2);
	vector l1_hat(i1,i2);
	vector l2_hat(i1,i2);
	vector linf(i1,i2);
	vector f1(i1,i2);
	vector f2(i1,i2);
	vector rho(i1,i2);
	vector nll(i1,i2);
	
PROCEDURE_SECTION

	calc_relative_age();
	calc_predicted_length();
	calc_objective_function();
	
FUNCTION calc_relative_age
	int i;
	dvariable d,tmp;
	for(i=i1;i<=i2;i++)
	{
		//d      = mfexp(-k*dt(i));
		//tmp    = (l2(i)-l1(i))/(l2(i)-d*l1(i));
		//age(i) = -1./k*log(tmp)*exp(A(i)*sigA);  //add random effect here
		//tmp      = -log(fabs(mu_linf-l1(i))/mu_linf) / k;
		//age(i)   = tmp*exp(A(i)*sigA);
		age(i) = muA*exp(A(i)*sigA);
	}

FUNCTION calc_predicted_length
	int i;
	dvariable t1,t3;
	sig2     = square(sig);
	sd2_linf = square(sd_linf);
	for(i=i1;i<=i2;i++)
	{
		f1(i)      = (1.-exp(-k*age(i)));
		f2(i)      = (1.-exp(-k*(age(i)+dt(i))));

		t1         = sd2_linf/(sig2+sd2_linf*(square(f1(i))+square(f2(i))));
		t3         = f1(i)*(l1(i)-mu_linf*f1(i))+f2(i)*(l2(i)-mu_linf*f2(i));
		linf(i)    = mu_linf+t1*t3;  
		l1_hat(i)  = mu_linf*f1(i);
		l2_hat(i)  = mu_linf*f2(i);
	}
	
FUNCTION calc_objective_function
	int i;
	dvariable s1,s2,t1,t2,t3,t9;
	dvariable r1,r2;
	for(i=i1;i<=i2;i++)
	{
		s1     = sqrt(sd2_linf*square(f1(i)) + sig2);
		s2     = sqrt(sd2_linf*square(f2(i)) + sig2);
		rho(i) = sd2_linf*f1(i)*f2(i)/(s1*s2);
		
		r1     = l1(i)-l1_hat(i);
		r2     = l2(i)-l2_hat(i);
		
		t1     = square(r1)/(s1*s1);
		t2     = 2.*rho(i) * (r1*r2) / (s1*s2);
		t3     = square(r2)/(s2*s2);
		t9     = 1.-square(rho(i));
		nll(i) = log(2.*M_PI)+log(s1)+log(s2)
		         +0.5*log(t9)
		         +(t1-t2+t3)/(2.*t9);
	}
	
	f = sum(nll) + 0.5*norm2(A);
	
REPORT_SECTION
	REPORT(f);
	REPORT(mu_linf);
	REPORT(sd_linf);
	REPORT(k);
	REPORT(sig);
	REPORT(muA);
	REPORT(sigA);
	REPORT(rho);
	REPORT(nll);
	REPORT(ii)
	REPORT(age);
	REPORT(dt);
	REPORT(linf);
	REPORT(l1);
	REPORT(l1_hat);
	REPORT(l2);
	REPORT(l2_hat);
	dvector l1_res = value((l1 - elem_prod(linf,f1))/sig);
	dvector l2_res = value((l2 - elem_prod(linf,f2))/sig);
	REPORT(l1_res);
	REPORT(l2_res);

TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(1600000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(1600000);
 

GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	//#include <statsLib.h>
	
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;

