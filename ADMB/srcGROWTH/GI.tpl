//  ******************************************************************
//  GI  A growth increment model for use with mark-recapture data.
//
//  Created by Martell on 2012-06-27.
//  Copyright (c) 2012. All rights reserved.
//  Comments:
//  ******************************************************************


DATA_SECTION
	init_adstring str_dataFile;
	!!ad_comm::change_datafile_name(str_dataFile);

	init_int nobs;
	init_ivector inobs(1,nobs);
	init_3darray data(1,nobs,1,inobs,1,6);
	init_int eof;
	//!! cout<<data(1)<<endl;
	!! if(eof != 999){cout<<"Error reading data\n"; exit(1);}
	
	imatrix iyr(1,nobs,1,inobs);
	imatrix loc(1,nobs,1,inobs);
	matrix l1(1,nobs,1,inobs);
	matrix l2(1,nobs,1,inobs);
	matrix dt(1,nobs,1,inobs);
	matrix dl(1,nobs,1,inobs);
	
	
	LOC_CALCS
		int i;
		for(i=1;i<=nobs;i++)
		{
			iyr(i) = ivector(column(data(i),1));
			loc(i) = ivector(column(data(i),2));
			l1(i)  = column(data(i),3);
			l2(i)  = column(data(i),4);
			dt(i)  = column(data(i),5)/365.25;
			dl(i)  = column(data(i),6);
		}
		//cout<<iyr<<endl;
	END_CALCS
	int nbins;
	!! nbins = 46;
	vector lbins(1,nbins);
	!! lbins.fill_seqadd(50,10);
	!! cout<<lbins<<endl;
	3darray tmat(1,nobs,1,nbins-1,1,nbins-1);
	3darray T(1,nobs,1,nbins-1,1,nbins-1);
	

PARAMETER_SECTION
	init_number mu_log_linf;
	init_number mu_log_k;
	init_number mu_log_phi;
	init_number sig_log_linf(3);
	init_number sig_log_k(3);
	init_number sig_log_phi(3);
	
	init_bounded_dev_vector log_linf(1,nobs,-5,5,2);
	init_bounded_dev_vector log_k(1,nobs,-5,5,2);
	init_bounded_dev_vector log_phi(1,nobs,-5,5,2);
	
	objective_function_value f;

	vector linf(1,nobs);
	vector k(1,nobs);
	vector phi(1,nobs);
	vector fvec(1,nobs);
	matrix dl_hat(1,nobs,1,inobs);
	matrix epsilon(1,nobs,1,inobs);
	matrix sig(1,nobs,1,inobs);
	matrix df(1,nobs,1,inobs);
	sdreport_number sd_mu_linf;

INITIALIZATION_SECTION
	mu_log_linf   5.0;
	mu_log_k     -2.0;
	mu_log_phi   -2.5;
	sig_log_linf -2.3;
	sig_log_k    -2.3;
	sig_log_phi  -2.3;

PROCEDURE_SECTION

	vbginc();
	calc_objective_function();
	if(mceval_phase())
	{
		calcLebesgueMeasure();
	} 
	
FUNCTION calcLebesgueMeasure
	/*
	This function is based on Rich Hillary's algorithm for calculating
	the relative overlap of each length partition based on parametric
	uncertainty and measurement error.  In this case, I have used a normal
	density function for measurement error.
	
	
	*/
	
	int h,i,j;
	static int nf = 0;
	double tau    = 1.0;         //Time interval
	double i1,i2,sig,eps;
	double j1,j2,l1i,l2i;
	double mu,nu;
	
	dvector tmp(1,2);
	if(nf == 0) 
	{
		tmat.initialize();
		T.initialize();
		cout<<"Initialized Arrays"<<endl;
	}
	nf++;
	
	for(h=1;h<=nobs;h++)
	{
		double ilinf = value(linf(h));
		double ik    = value(k(h));
		for(i=1;i < nbins;i++)
		{
			i1  = lbins(i);
			i2  = lbins(i+1);
			sig = value(phi(h))*0.5*(i1+i2);
			eps = randn(nf+i)*sig;
			
			for(j=1;j < nbins;j++)
			{
				j1    = lbins(j);
				j2    = lbins(j+1);
				
				l1i   = i1 + ginc(ilinf,ik,tau,i1) + eps;
				l2i   = i2 + ginc(ilinf,ik,tau,i2) + eps;
			
				/* Calculate Lebesque measure of intersection */
				if(l1i>j1 && l2i<j2)
				{
					tmat(h,i,j) += 1.0;
				}
				else
				{
					tmp(1)        = max(l1i,j1);
					tmp(2)        = min(l2i,j2);
					
					if(tmp(1)<tmp(2)) mu = tmp(2)-tmp(1); else mu = 0.0;
						
					nu            = l2i-l1i;
					tmat(h,i,j) += mu/nu;
				}
			}
		}
		T(h) = tmat(h)/nf;
	}
	if(nf==1000) cout<<"Transition Matrix\n"<<T(1)<<endl;
	
FUNCTION double ginc(const double& linf, const double& k, const double& tau,const double& l1)
	return( (linf-l1)*(1.-exp(-k*tau)) );

FUNCTION vbginc
	linf    = mfexp(mu_log_linf + log_linf);
	k       = mfexp(mu_log_k + log_k);
	phi     = mfexp(mu_log_phi + log_phi);
	sd_mu_linf = mean(linf);
	
	int i;
	for(i=1;i<=nobs;i++)
	{
		dl_hat(i)  = elem_div(elem_prod((linf(i)-l1(i)),(1.-exp(-k(i)*dt(i)))),dt(i));
		epsilon(i) = dl(i) - dl_hat(i);
		sig(i)     = sqrt( 88.6 + square(phi(i)*l1(i)) );
		df(i)      = 1./(1./88.6 + l1(i)/phi(i));
	}

FUNCTION dvariable deprecated_dstudent_t( const dvar_vector& residual, const dvar_vector& df)
	{
		double pi =  3.141593;
		dvar_vector t1 = 0.5*(df+1);
		dvar_vector t2 = gammln(t1);
		dvar_vector t3 = 0.5*log(df*pi)+gammln(0.5*df);
		dvar_vector t4 = elem_prod(t1,log(1+elem_div(square(residual),df)));
		dvariable pdf = sum(t3+t4-t2);
		return( pdf );
	}


FUNCTION calc_objective_function
	int i;
	fvec.initialize();
	for(i=1;i<=nobs;i++)
	{
		//fvec(i) = dstudent_t(epsilon(i),df(i));
		fvec(i) = dnorm(epsilon(i),sig(i));
	}
	f  = sum(fvec);
	f += dnorm(log_k,mfexp(sig_log_k));
	f += dnorm(log_linf,mfexp(sig_log_linf));
	f += dnorm(log_phi,mfexp(sig_log_phi));
	f += dgamma(1/mfexp(sig_log_k),1.01,1.01);
	f += dgamma(1/mfexp(sig_log_linf),1.01,1.01);
	f += dgamma(1/mfexp(sig_log_phi),1.01,1.01);

REPORT_SECTION
	REPORT(nobs);
	REPORT(linf);
	REPORT(k);
	REPORT(phi);
	REPORT(l1);
	REPORT(dl);
	REPORT(dl_hat);
	REPORT(dt);
	REPORT(epsilon);
	REPORT(sig);
	REPORT(df);
	dmatrix nu = value(elem_div(df,df-2.));
	REPORT(nu);
	REPORT(loc);
	REPORT(iyr);

TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 

GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	#include <statsLib.h>
	
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
FINAL_SECTION
	ofstream ofs("TransitionMatrix.txt");
	int i;
	ofs<<"lbins\n"<<lbins<<endl;
	ofs<<"T"<<endl;
	for(i=1;i<=nobs;i++)
	{
		ofs<<"#"<<iyr(i,1)<<endl<<T(i)<<endl;
	}
	

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

