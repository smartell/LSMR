// ************************************************************************* //
// LSMR: Length Structured Mark-Recapture model.                             //
//                                                                           //
// Created by Steven James Dean Martell on 2011-10-26.                       //
// Copyright (c) 2011. All rights reserved.                                  //
//                                                                           //
// PROJECT LOG:  DATE       HOURS  COMMENTS                                  //
//               Sep 23/11	  4    Simulation model in R			         //
//               Sep 26/11    8    Set up HBC project & Rcode                //
//               Sep 27/11	  8    Simulation model in R                     //
//               Oct 27/11    8    Simulation model in R & ADMB code         //
//               Oct 31/11    4    LSMR.tpl code  (PAID $11,000)             //
//               Apr 26/12    4    LSMR code development.                    //
//                                                                           //
//                                                                           //
// ************************************************************************* //


DATA_SECTION
	
	
	// Read input data from data file //
	init_int syr;   //first year
	init_int nyr;   //last year
	init_number dt; //time-step
	
	// Array dimensions //
	init_int nrow;  //number of rows in observed data
	init_int ncol;  //number of cols in observed data
	init_int nbin;	//number of length intervals
	int nx;			//number of length bins (nbin-1)
	int nr;			//number of rows in predicted arrays
	!! nr = int((nyr-syr+1)/dt);
	!! nx = nbin-1;
	init_vector xbin(1,nbin);  //length-intervals (not mid-points)
	vector xmid(1,nx);
	!! xmid = xbin(1,nbin-1)+0.5*first_difference(xbin);
	
	init_matrix i_C(1,nrow,1,ncol); //Catch-at-length
	init_matrix i_M(1,nrow,1,ncol); //New Marks-at-length
	init_matrix i_R(1,nrow,1,ncol); //Recaptures-at-length
	
	matrix C(1,nrow,1,nx);
	matrix M(1,nrow,1,nx);
	matrix R(1,nrow,1,nx);
	LOC_CALCS
		for(int i=1;i<=nrow;i++)
		{
			C(i)(1,nx) = i_C(i)(3,ncol).shift(1);
			M(i)(1,nx) = i_M(i)(3,ncol).shift(1);
			R(i)(1,nx) = i_R(i)(3,ncol).shift(1);
		}
	END_CALCS
	
	init_int eof;
	!! if(eof!=999){cout<<"Error reading data\n"; exit(1);}
	!! cout<< "____ End of Data file! ____\n"<<endl;

	
PARAMETER_SECTION
	init_number log_Ninit;
	init_number log_rbar;
	init_number log_fbar;
	init_number log_m_linf(4);
	!! log_m_linf = log(0.08);
	!! log_Ninit = 7;
	!! log_rbar = 9;
	!! log_fbar = log(0.1);
	//Selectivity parameters
	init_number log_lx(2);
	init_number log_gx(2);
	!! log_lx = log(10);
	!! log_gx = log(2.5);
	
	
	init_bounded_dev_vector init_log_Ninit_devs(1,nx,-15,15,1);
	init_bounded_dev_vector log_rec_devs(syr+1,nyr,-15,15,3);
	init_bounded_dev_vector log_fi_devs(1,nrow,-15.0,15.0,2);
	
	
	
	objective_function_value f;
	
	//Variables for growth.
	number linf;
	number k;
	number beta;
	number dl_cv;
	number fpen;
	
	number m_linf;
	number lx;				//length at 50% selectivity
	number gx;				//std in length at 50% selectivity
	vector mx(1,nx);		//Mortality rate at length xmid
	vector sx(1,nx);		//Selectivity at length xmid
	vector log_rt(syr+1,nyr);
	vector fi(1,nrow);		//capture probability in period i.
	matrix N(1,nr+1,1,nx);	//Numbers(time step, length bins)
	matrix T(1,nr+1,1,nx);  //Marks-at-large (time step, length bins)
	
	// predicted observations
	matrix Chat(1,nrow,1,nx);
	matrix Mhat(1,nrow,1,nx);	//Newly marked animals
	matrix Rhat(1,nrow,1,nx);
	
RUNTIME_SECTION
    maximum_function_evaluations 500,1500,2500,25000,25000
    convergence_criteria 0.01,1.e-4,1.e-5,1.e-5


PROCEDURE_SECTION
	initialize_model();
	calc_survival_at_length();
	calc_selectivity_at_length();
	calc_numbers_at_length();
	calc_catch_at_length();
	calc_objective_function();
	
	
FUNCTION initialize_model
  {
	fpen = 0.;
	
	//Set up the initial states.
	N.initialize();
	T.initialize();
	N(1) = mfexp(log_Ninit+init_log_Ninit_devs);
	
	fi  = mfexp(log_fbar + log_fi_devs);
	
	log_rt = log_rbar + log_rec_devs;
	m_linf=mfexp(log_m_linf);
	lx = mfexp(log_lx);
	gx = mfexp(log_gx);
	
	linf=40;
	k=0.18;
	beta = 0.75;
	dl_cv = 0.75;
		
	
  }

FUNCTION calc_survival_at_length
  {
	/*
	This function calculates the length-specific survival rate
	based on the Lorenzen function, where survival increases
	with increasing length.
	
	mortality rate at length x  mx=m.linf*linf/xbin
	note that m_linf is an annual rate.
	*/
	mx = m_linf*linf/xmid;
	
  }

FUNCTION calc_selectivity_at_length
  {
	/*
	This function calculates the length-specific selectivity.
	The parametric option is a logistic curve (plogis(x,lx,gx))
	*/
	sx.initialize();
	
	sx = 1./(1+mfexp(-(xmid-lx)/gx));
	//cout<<sx<<endl;
  }

FUNCTION calc_numbers_at_length
  {
	/*	This function updates the numbers at length
		over time as a function of the numbers at length
		times the survival rate.
	  	N_{t+dt} = elem_prod(N_{t},exp(-mx)) * LTM
	
	  	These are the total number of fish (tagged + untagged)
		at large in the population.  N is used in the observation
		models to determine number of captures and recaptures.
	*/
	int i,j,t,im;
	
	dvar_matrix LTM = dv_LTM(xbin,linf,k,beta);
	//dvar_matrix LTM = dLTM(xbin,linf,k*dt,dl_cv);
	dvariable a=1/(0.2*0.2);  //a=1/cv^2
	dvariable b=10/a;		  //b=mu/a
	dvar_vector p1 = dgamma(xmid,a,b);
	p1/=sum(p1);
	
	i=0; j=1;
	for(t=syr;t<=nyr;t++)
	{
		for(im=1;im<=int(1/dt);im++)
		{
			i++;
			//recruitment if in period 1
			if(im==1 && t > syr)
			{
				N(i)+= mfexp(log_rt(t))*p1;
			}
			
			N(i+1)=elem_prod(N(i),mfexp(-mx*dt))*LTM;
			
			//Set up accounting for predicted marks at large
			//By adding observed marks to Mhat.
			//Mhat(i+1)= elem_prod(Mhat(i),exp(-mx*dt))*LTM + M(i)
			dvector newMarks(1,nx);
			newMarks.initialize();
			if(t==i_M(j,1) && im==i_M(j,2)) 
			{
				newMarks = M(j++);
				if(j>nrow) j=nrow;
			}
			
			T(i+1) = elem_prod(T(i),mfexp(-mx*dt)) * LTM + newMarks;
			
			//cout<<t<<"\t"<<i<<"\t"<<sum(N(i))<<endl;
		} //end of im
	} // end of t
  }

FUNCTION calc_catch_at_length
  {
	/*
		This function calculates the total catch at length
		at each time step in the model. It also calculates
		the number of new marks released for fish greater 
		than 150mm (15 cm), and the number of recaptures.
		
		i = index for time
		j = index for length
		f_i = capture probabilty in year i.
		C_ij = f_i*sx_j*N_ij
		N_ij = total number of (marked & unmarked) fish at large
		T_ij = number of marked fish at large.
		markRate = T_ij/N_ij
		
	*/
	int i,j,t,im,ir;
	
	for(i=1;i<=nrow;i++)
	{
		t  = i_C(i,1);		//year
		im = i_C(i,2);		//im
		ir = int((t-syr)/dt+im);
		//fi(i) = sum(C(i)) / (N(ir)*sx);
		//Chat(i) = fi(i)*elem_prod(N(ir),sx);
		dvar_vector zj = mx*dt;
		Chat(i) = elem_prod(elem_div(elem_prod(fi(i)*sx,1.-exp(-zj)),zj),N(ir));
		//cout<<t<<"\t"<<im<<"\t"<<ir<<"\tfi = "<<fi(i)<<endl;
		
		//Might also think about trying:
		Mhat(i)(12,nx) = fi(i)*elem_prod(N(ir)(12,nx)-T(ir)(12,nx),sx(12,nx));
		Rhat(i) = fi(i)*elem_prod(T(ir),sx);
		
		
		//New Marks based on proportion of Chat that is unmarked (1-T/N)
		//dvar_vector markRate = 1. - vposfun(1. - elem_div(T(ir),N(ir)),0.01,fpen);
		//cout<<min(N(ir))<<endl;
		//Mhat(i)(13,nx) = elem_prod( Chat(i)(13,nx), 1.-markRate(13,nx) );
		
		
		//Recaptures based on proportion of Chat that is marked (T/N)
		//Rhat(i) = elem_prod( Chat(i), markRate );
		
	}
  }

FUNCTION dvar_vector vposfun(const dvar_vector &x,const double eps, dvariable& fpen)
  {
	int i,i1,i2;
	i1=x.indexmin();
	i2=x.indexmax();
	dvar_vector tmp(i1,i2);
	for(i=i1;i<=i2;i++)
	{
		tmp(i) = posfun(x(i),eps,fpen);
	}
	return(tmp);
  }


FUNCTION calc_objective_function;
  {
	int i;
	dvar_vector fvec(1,3);
	dvar_vector pvec(1,3);
	
	fvec.initialize();
	pvec.initialize();
	//Priors (pvec)
	if(!last_phase())
	{
		pvec(1) = norm2(init_log_Ninit_devs);
		pvec(2) = norm2(log_rec_devs);
		pvec(3) = dnorm(log_fbar,log(0.1),0.05);
	}
	else
	{
		pvec(1) = dnorm(init_log_Ninit_devs,0,2.5);
		pvec(2) = dnorm(log_rec_devs,0,2.5);
		pvec(3) = dnorm(log_fbar,log(0.1),2.5);
	}
	
	//Likelihoods (fvec);
	fvec(1) = 10.*norm2(sqrt(C)-sqrt(Chat));
	//fvec(2) = 10.*norm2(M-Mhat);
	//fvec(3) = 10.*norm2(R-Rhat);
	
	/* Likelihood for new marks */
	for(i=1;i<=nrow;i++)
	{
		fvec(2) += dmultinom(M(i),Mhat(i)+EPS);
		fvec(3) += dmultinom(R(i),Rhat(i)+EPS);
	}
	
	
	if(fpen > 0 ) cout<<"Fpen = "<<fpen<<endl;
	f = sum(fvec) + sum(pvec);//	 + 1.e6*fpen;
  }



FUNCTION dvar_vector dgamma(dvector& x, dvariable& a, dvariable& b)
  {
	//returns the gamma density with a & b as parameters
	RETURN_ARRAYS_INCREMENT();
	dvariable t1 = 1./(pow(b,a)*exp(gammln(a)));
	dvar_vector t2 = (a-1.)*log(x)-x/b;
	RETURN_ARRAYS_DECREMENT();
	return(t1*exp(t2));
  }

FUNCTION dvar_matrix dLTM(dvector& x, const dvariable &linf, const dvariable &k, const dvariable &cv)
  {
	/*
		This function returns a length transition matrix (P_ij) with
		the probability of an individual of length (x_i) growing into
		the length interval (x_j). Length transition is based on the 
		gamma distribution dgamma(x,a,b_l), where a=1/cv^2 and b_l
		is given by b_l = (linf-x_j)*(1-exp(-k*dt))
		
		Args:
			(x)  	- discrete length interval boundaries.
			(linf)	- asymptotic length
			(k)		- Brody growth coefficient over time interval dt.
			(cv) 	- coefficient of variation in the size increment.
			
			Use the following for the mean growth increment to ensure
			the function (dl) remains positive. Beyone Linf, growth 
			increment decays exponentially.
			dl	<- log(exp((linf-xp)*(1-exp(-k)))+1)  to ensure dl is positive
			
	*/
	
		RETURN_ARRAYS_INCREMENT();
		int i,j,n;
		n = size_count(x);
		dvector xm(1,n-1);  //mid points of the length intervals
		xm = x(1,n-1)+0.5*first_difference(x);
		dvariable pen;
		dvariable alpha = 1./square(cv);
		dvar_vector beta(1,n-1);
		dvar_vector dl(1,n-1);
		dvar_matrix P(1,n-1,1,n-1);
		cout<<"Ok"<<endl;
		for(j=1;j<=n-1;j++)
		{
			//dvariable t1 = posfun(linf-xm(j),1./n,pen);
			//beta(j) = t1 * (1.-mfexp(-k));
			beta(j) = log(mfexp( (linf-xm(j))*(1.-mfexp(-k)) )+1.0);
		}
		dl = alpha*beta;
		
		for(i=1;i<=n-1;i++)
		{
			dvar_vector t2(i,n);
			for(j=i;j<=n;j++)
			{
				double dx = x(j)-x(i);
				t2(j) = cumd_gamma(dx,dl(i));
			}
			P(i)(i,n-1) = first_difference(t2);
			P(i) /= sum(P(i));
		}
		
		//cout<<P<<endl;
		
		RETURN_ARRAYS_DECREMENT();
		return(P);
  }


FUNCTION dvar_matrix dv_LTM(dvector& x, const dvariable &linf, const dvariable &k, const dvariable &beta)
  {
	/*This function computes the length transition matrix.*/
	/*
	- cumd_gamma(x,a) is the same as Igamma(a,x) in R.
	- If length interval > linf, then assume now further growth.
	- Note the use of posfun to ensure differentiable.
	*/
	RETURN_ARRAYS_INCREMENT();
	int i,j;
	double dx;
	int n = size_count(x);     //number of length intervals
	dvar_vector alpha(1,n-1);
	dvar_matrix P(1,n-1,1,n-1);
	P.initialize();
	
	//Growth increment
	dvariable pen;
	dvar_vector xm(1,n-1);
	dvar_vector dl(1,n-1);
	xm = x(1,n-1)+0.5*first_difference(x);
	for(j=1;j<=n-1;j++)
	{
		//dvariable t2=posfun(linf-xm(j),1./n,pen);
		//dl(j) = (t2)*(1.-exp(-k*dt));
		dl(j) = log(mfexp( (linf-xm(j))*(1.-mfexp(-k*dt)) )+1.0);
	}
	alpha = dl/beta;
	
	
	for(i=1;i<=n-1;i++)
	{
		dvar_vector t1(i,n);
		for(j=i;j<=n;j++)
		{
			dx = x(j)-x(i);
			//cout<<j<<" "<<x(i)<<" "<<x(j)<<" "<<cumd_gamma(dx,alpha(i))<<endl;
			t1(j) = cumd_gamma(dx,alpha(i));
		}
		P(i)(i,n-1) = first_difference(t1);
		P(i) /= sum(P(i));
	}
	RETURN_ARRAYS_DECREMENT();
	//cout<<P<<endl;
	return(P);
  }
	
REPORT_SECTION
	REPORT(f);
	REPORT(log_Ninit);
	REPORT(log_rbar);
	REPORT(log_fbar);
	REPORT(log_m_linf);
	REPORT(mx);
	
	REPORT(log_rt);
	REPORT(fi);
	
	REPORT(N);
	REPORT(T);
	REPORT(i_C);
	REPORT(i_M);
	REPORT(i_R);
	REPORT(Chat);
	REPORT(Mhat);
	REPORT(Rhat);


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
	#include <contrib.h>
	//#include <stats.cxx>
	//#include <martool.cxx>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
	dvariable dpois(const dvector& k, const dvar_vector& lambda)
	{
		RETURN_ARRAYS_INCREMENT();
		int i;
		int n = size_count(k);
		dvariable nll=0;
		for(i = 1; i <= n; i++)
		{
			nll -= k(i)*log(lambda(i))+lambda(i)+gammln(k(i)+1.);
		}
		RETURN_ARRAYS_DECREMENT();
		return nll;
	}
	
	
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

