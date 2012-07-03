//  ******************************************************************
//  Bayesian Hierarchical Growth Model (BHGM)
//  This model is based on the Zhang et. al. (2009) model
//  for Northern Abalone.  It also has an option for comparing
//  to the method of Fabens, as described in Hillary (2011) with
//  the student-t distribution for the probability density function.
//
//  Created by Steven James Dean Martell on 2011-11-03.
//  Copyright (c) 2011. All rights reserved.
//  Comments:
//  ******************************************************************


DATA_SECTION
	int method;
	LOC_CALCS
		int on,opt;
		method=1;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-method",opt))>-1)
		{
			method=atoi(ad_comm::argv[on+1]);
		}
		cout<<method<<endl;
		if(method !=1 && method !=2)
		{
			cout<<"Available methods are \n \t Fabens = 1 \n \t Zhang = 2\n";
			exit(1);
		}
	END_CALCS
	
	int nobs;
	!! nobs = 100;  //total number of records is 8315
	init_matrix gdata(1,nobs,1,3);
	vector ltag(1,nobs);
	vector dl(1,nobs);
	vector dt(1,nobs);
	
	!! ltag = column(gdata,1);
	!! dl   = column(gdata,2);
	!! dt   = column(gdata,3);	
	//!!cout<<gdata<<endl;
	int phz;

PARAMETER_SECTION
	//Hyper-parameters
	init_number log_linf;
	init_number log_k;
	init_number log_phi;
	
	!! log_linf = log(mean(ltag+dl));
	!! log_k    = log(0.15);
	!! log_phi  = log(50);
	
	//Random-effects for Zhang method=2
	!! if(method==2) phz=2; else phz=-1;
	init_number log_sig_linf(phz);
	random_effects_vector delta(1,nobs,phz);
	
	
	objective_function_value f;
	
	number linf;
	number k;
	number phi;
	number sig_linf;
	vector dl_hat(1,nobs);
	vector epsilon(1,nobs);
	sdreport_number sd_dl1;
	
PROCEDURE_SECTION
	
	if(method==1) fabens_method();
	if(method==2) zhang_method();
	sd_dl1 = dl_hat(1);
	
FUNCTION zhang_method
	linf = mfexp(log_linf);
	k    = mfexp(log_k);
	phi  = mfexp(log_phi);
	sig_linf = mfexp(log_sig_linf);
	
	dvar_vector linf_i=linf*exp(delta * sig_linf);
	dl_hat = elem_prod(linf_i-ltag,1.-exp(-k*dt));
	epsilon = dl - dl_hat;
	dvar_vector sig = sqrt(31.8 + square(phi*ltag));
	f = dnorm(epsilon,sig);
	f+= dnorm(delta,1);
	
FUNCTION fabens_method
	linf = mfexp(log_linf);
	k    = mfexp(log_k);
	phi  = mfexp(log_phi);
	
	dl_hat = elem_prod(linf-ltag,1.-exp(-k*dt));
	epsilon = dl - dl_hat;
	dvar_vector df = 1./(1./100. + ltag/phi);
	//cout<<epsilon(1,1)<<" "<<df(1,1)<<" "<<dstudent_t(epsilon(1,1),df(1,1))<<endl;
	//f = dstudent_t(epsilon,df);
	
	//sigm^2 from Coggins 2006 was 31.8mm
	dvar_vector sig = sqrt(31.8 + square(phi*ltag));
	f = dnorm(epsilon,sig);
	//f= norm2(epsilon);
	
FUNCTION dvariable dstudent_t( const dvar_vector& residual, const dvar_vector& df)
	{
		double pi =  3.141593;
		dvar_vector t1 = 0.5*(df+1);
		dvar_vector t2 = gammln(t1);
		dvar_vector t3 = 0.5*log(df*pi)+gammln(0.5*df);
		dvar_vector t4 = elem_prod(t1,log(1+elem_div(square(residual),df)));
		dvariable pdf = sum(t3+t4-t2);
		return( pdf );
	}


FUNCTION dvariable dnorm(const dvar_vector& residual, const double std)
	{
		double pi=3.141593;
		long n=size_count(residual);
		dvariable SS=norm2(residual);
		return(n*(0.5*log(2.*pi)+log(std))+0.5*SS/(std*std));
		//return(n*(0.5*log(2.*pi)-log(std))+0.5*SS*(std*std));
	}



FUNCTION dvariable dnorm( const dvar_vector& residual, const dvar_vector std )
	{
		RETURN_ARRAYS_INCREMENT();
		int n=size_count(residual);
		double pi =  3.141593;
		dvar_vector var=elem_prod(std,std);
		dvar_vector SS=elem_prod(residual,residual);
		dvariable tmp=0.5*n*log(2.*pi)+sum(log(std))+sum(elem_div(SS,2.*var));
		RETURN_ARRAYS_DECREMENT();
		return( tmp );
	}


REPORT_SECTION
	REPORT(linf);
	REPORT(k);
	REPORT(phi);
	REPORT(epsilon);
	REPORT(ltag);
	REPORT(dl);
	REPORT(dt);
	REPORT(dl_hat);

TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
	gradient_structure::set_MAX_NVAR_OFFSET(250504);
 

GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	//#include <stats.cxx>
	//#include <martool.cxx>
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

