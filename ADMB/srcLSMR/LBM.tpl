//  ******************************************************************
//  LBM
//  A test template for a length based model (LBM) using the 
//  incomplete gamma function (cumd_gamma(x,a)/gamma(x))
//  
//  Created by Steven James Dean Martell on 2011-05-31.
//  Copyright (c) 2011. All rights reserved.
//  Comments:  
//  ******************************************************************


DATA_SECTION
	init_int syr;
	init_int nyr;
	init_int lb;
	init_int ub;
	init_int bw;
	
	ivector iyr(syr,nyr);
	!! iyr.fill_seqadd(syr,1);
	
	//set up length intervals lx and midpoints xm.
	int nbins;
	!! nbins = (int)(ub-lb)/bw;
	vector xl(1,nbins+1);
	vector xm(1,nbins);
	!! xl.fill_seqadd(lb,bw);
	!! xm = xl(1,nbins)+0.5*bw;
	
	//Growth parameters
	init_number linf;
	init_number k;
	init_number beta;
	
PARAMETER_SECTION
	init_number dum;
	objective_function_value f;

	//matrix P(1,nbins,1,nbins);	//length transition matrix
	
PRELIMINARY_CALCS_SECTION
  dvar_matrix P = getLengthTransitionMatrix(linf,k,beta,xl);

PROCEDURE_SECTION
	
	
FUNCTION dvar_matrix getLengthTransitionMatrix(const double& linf, const double& k, const double& beta, dvector& xl)
	/*
	This routine returns the length transition matrix P(ij),
	which is defined as the probability of an individual
	in length interval i growing into length interval j.
	This function assumes that the increase in length is 
	based on the gamma distribution Γ(x|α_l,β), where x 
	represents the change in length Δl=(L∞-l)(1-exp(-k)).
	The mean change in length is given by α_lβ and the variance
	by βΔl.  
	
	Args: 	linf = asymptotic length
			k = brody growth coefficient
			beta = 1/coefficient of variation
			xl = vector of length intervals (lower bound and upper bounds)
	
	Assumptions:  Fish greater than L∞ are assumed to remain in the
	same size class.
	
	With the new version of ADMB 10.2, there is an issue with the 
	cumd_gamma function calculating large values of the growth
	increment (dx).  It would crash in cases where cumd_gamma(dx,alpha)
	where dx > 250 mm. Simpliest solution to this problem is to switch
	scales from mm to cm.
	*/
	
	// Change in length (this is the x in the gamma distribution)
	int i,j,n;
	n = size_count(xl)-1;
	dvector dl = (linf-xm)*(1.-exp(-k));
	dvector alpha = dl/beta;
	
	// Get midpoints (xm) of length intervals
	dvector xm = xl(1,n)+0.5*first_difference(xl);
	
	dvar_matrix P(1,n,1,n);
	P.initialize();
	
	for(i=1;i<=n;i++)
	{
		for(j=1;j<=i;j++)
		{
			double dx = xl(i)-xl(j);
			cout<<i<<"\t"<<j<<"\t"<<dx<<"\t"<<alpha(i)<<endl;
			
			double dwx = 0.5*(xl(i+1)-xl(i));
			P(i,j) = cumd_gamma(dx+dwx,alpha(i))-cumd_gamma(dx-dwx,alpha(i));
			
			cout<<P(i,j)<<endl;
		}
	}
	cout<<endl<<P<<endl;
	return(P);

	
REPORT_SECTION


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

