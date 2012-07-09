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
//               Jul 08/12    8    LSMR simulation routine.                  //
//                                                                           //
// ************************************************************************* //


DATA_SECTION
	init_adstring data_file;
	init_adstring control_file;
	
	int SimFlag;
	int rseed;
	LOC_CALCS
		SimFlag = 0;
		rseed   = 1;
		int on,opt;
		//the following line checks for the "-SimFlag" command line option
		//if it exists the if statement retreives the random number seed
		//that is required for the simulation model
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			SimFlag = 1;
			rseed   = atoi(ad_comm::argv[on+1]);
			if(SimFlag) cout<<"In Simulation Mode\n";
		}
		
	END_CALCS
	
	!! ad_comm::change_datafile_name(data_file);
	
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
	
	
	init_int eof;
	!! if(eof!=999){cout<<"Error reading data\n"; exit(1);}
	!! cout<< "____ End of Data file! ____\n"<<endl;
	
	
	!! ad_comm::change_datafile_name(control_file);
	init_int npar
	init_matrix theta_control(1,npar,1,7);
	matrix trans_theta_control(1,7,1,npar);
	vector theta_ival(1,npar);
	vector theta_lbnd(1,npar);
	vector theta_ubnd(1,npar);
	ivector theta_phz(1,npar);
	ivector theta_prior(1,npar);
	LOC_CALCS
		trans_theta_control = trans(theta_control);
		theta_ival = trans_theta_control(1);
		theta_lbnd = trans_theta_control(2);
		theta_ubnd = trans_theta_control(3);
		theta_phz  = ivector(trans_theta_control(4));
		theta_prior = ivector(trans_theta_control(5));
	END_CALCS
	
	init_int nflags;
	init_vector flag(1,nflags);
	// 1) minimum size of fish for tagging.
	
	
	
	// Pre-processing some of the data.
	ivector j_min(1,nrow);
	matrix C(1,nrow,1,nx);
	matrix M(1,nrow,1,nx);
	matrix R(1,nrow,1,nx);
	LOC_CALCS
		for(int i=1;i<=nrow;i++)
		{
			C(i)(1,nx) = i_C(i)(3,ncol).shift(1);
			M(i)(1,nx) = i_M(i)(3,ncol).shift(1);
			R(i)(1,nx) = i_R(i)(3,ncol).shift(1);
			
			//set up index for minimum size for tagging fish.
			for(int j=1;j<=nx;j++)
			{
				if( M(i,j)>0 || xbin(j) > flag(1) )
				{
					j_min(i)=j;
					break;
				}
			}
		}
	END_CALCS
	
	
PARAMETER_SECTION
	init_bounded_number_vector theta(1,npar,theta_lbnd,theta_ubnd,theta_phz);
	number log_ddot_r;
	number log_bar_r;
	number log_bar_f;
	number m_infty;

	//Variables for growth.
	number l_infty;
	number vbk;
	number beta;
	
	//Variables for size distribution of new recruits
	number mu_r;
	number cv_r
	
	//Selectivity parameters
	init_number log_lx(2);
	init_number log_gx(2);
	
	init_bounded_dev_vector ddot_r_devs(1,nx,-15,15,1);
	init_bounded_dev_vector bar_r_devs(syr+1,nyr,-15,15,3);
	init_bounded_dev_vector bar_f_devs(1,nrow,-15.0,15.0,2);
	
	
	
	objective_function_value f;
	

	number dl_cv;
	number fpen;
	
	number m_linf;
	number lx;				//length at 50% selectivity
	number gx;				//std in length at 50% selectivity
	vector mx(1,nx);		//Mortality rate at length xmid
	vector sx(1,nx);		//Selectivity at length xmid
	vector rx(1,nx);		//size pdf for new recruits
	
	vector log_rt(syr+1,nyr);
	vector fi(1,nrow);		//capture probability in period i.
	matrix N(1,nr,1,nx);	//Numbers(time step, length bins)
	matrix T(1,nr,1,nx);	//Marks-at-large (time step, length bins)
	matrix P(1,nx,1,nx);	// Size-Transition Matrix
	
	// predicted observations
	matrix Chat(1,nrow,1,nx);
	matrix Mhat(1,nrow,1,nx);	//Newly marked animals
	matrix Rhat(1,nrow,1,nx);
	
INITIALIZATION_SECTION
	theta theta_ival;
	

RUNTIME_SECTION
    maximum_function_evaluations 500,1500,2500,25000,25000
    convergence_criteria 0.01,1.e-4,1.e-5,1.e-5

PRELIMINARY_CALCS_SECTION
	if(SimFlag)
	{
		cout<<"******************************"<<endl;
		cout<<"**                          **"<<endl;
		cout<<"** Running Simulation Model **"<<endl;
		cout<<"** Random seed = "<<rseed<<"        **"<<endl;
		cout<<"**                          **"<<endl;
		cout<<"******************************"<<endl;
		runSimulationModel(rseed);
		
	}


PROCEDURE_SECTION
	initParameters();
	calcSurvivalAtLength();
	calcSizeTransitionMatrix();
	initializeModel();
	calcNumbersAtLength();
	calcSelectivityAtLength();
	calcNewMarksAndRecaptures();
		
	//calc_catch_at_length();
	calc_objective_function();
	
FUNCTION void runSimulationModel(const int& seed)
  {
	int i,j;
	random_number_generator rng(seed);
	
	dvector tmp_ddot_r_devs = value(ddot_r_devs);
	dvector tmp_bar_r_devs  = value(bar_r_devs);
	tmp_ddot_r_devs.fill_randn(rng);
	tmp_bar_r_devs.fill_randn(rng);
	
	ddot_r_devs = 0.6*tmp_ddot_r_devs;
	bar_r_devs  = 0.6*tmp_bar_r_devs;
	
	/* Capture probabilities */
	//Et <- h1+(h2-h1)*cos(0.5+(sample.yrs-t)/(T-t)*pr*pi)^2
	double h1 = -1.;
	double h2 =  1.;
	for(i=1;i<=nrow;i++)
	{
		bar_f_devs(i) = h1 + (h2-h1)*square(sin((i-1.)/(nrow-1.)*1.5*PI));
	}
	
	initParameters();
	calcSurvivalAtLength();
	calcSizeTransitionMatrix();
	initializeModel();
	calcNumbersAtLength();
	calcSelectivityAtLength();
	calcNewMarksAndRecaptures();
	C=value(Chat);
	M=value(Mhat);
	R=value(Rhat);
	
	for(i=1;i<=nrow;i++)
	{
		for(j=1;j<=nx;j++)
		{
			C(i,j) = randnegbinomial(C(i,j),8.0,rng);
			if(M(i,j)>0)
			{
				M(i,j) = randnegbinomial(M(i,j),8.0,rng);
			}
			if(R(i,j)>0)
			{
				R(i,j) = randnegbinomial(R(i,j),8.0,rng);
			}
			
		}
		i_C(i)(3,ncol) = C(i)(1,nx).shift(3);
		i_M(i)(3,ncol) = M(i)(1,nx).shift(3);
		i_R(i)(3,ncol) = R(i)(1,nx).shift(3);
	}
  }

FUNCTION initParameters
  {
	/* Leading parameters */
	log_ddot_r = theta(1);
	log_bar_r  = theta(2);
	log_bar_f  = theta(3);
	m_infty    = theta(4);
	l_infty    = theta(5);
	vbk        = theta(6);
	beta       = theta(7);
	mu_r       = theta(8);
	cv_r       = theta(9);
	
	/* Selex parameters */
	lx         = mfexp(log_lx);
	gx         = mfexp(log_gx);
	
	
  }

FUNCTION calcSizeTransitionMatrix
	/*
	This function calls the necessary routines to compute the Size Transition Matrix (P)
	*/
	
	P = calcLTM(xbin,l_infty,vbk,beta);
	

FUNCTION initializeModel
  {
	int i;
	//Set up the initial states.
	N.initialize();
	T.initialize();
	
	
	/* Initialize numbers at length at the first time step. */
	dvariable a=1/(cv_r*cv_r);  //a=1/cv^2
	dvariable b=mu_r/a;		  //b=mu/a
	rx  = dgamma(xmid,a,b);
	rx /= sum(rx);
	
	for(i=1;i<=nx;i++)
	{
		dvar_vector init_r = mfexp(log_ddot_r + ddot_r_devs(i)) * rx;
		N(1) = elem_prod(N(1),mfexp(-mx)) * P + init_r;
	}
	
	
	/* Annual recruitment */
	log_rt = log_bar_r + bar_r_devs;
	
	
	/* Sampling effort at each time step. */
	fi  = mfexp(log_bar_f + bar_f_devs);		
  }

FUNCTION calcSurvivalAtLength
  {
	/*
	This function calculates the length-specific survival rate
	based on the Lorenzen function, where survival increases
	with increasing length.
	
	mortality rate at length x  mx=m.linf*linf/xbin
	note that m_linf is an annual rate.
	*/
	mx = m_infty*l_infty/xmid;
  }

FUNCTION calcSelectivityAtLength
  {
	/*
	This function calculates the length-specific selectivity.
	The parametric option is a logistic curve (plogis(x,lx,gx))
	*/
	sx.initialize();
	
	sx = 1./(1+mfexp(-(xmid-lx)/gx));
	//cout<<sx<<endl;
  }

FUNCTION calcNumbersAtLength
  {
	/*	This function updates the numbers at length
		over time as a function of the numbers at length
		times the survival rate.
	  	N_{t+dt} = elem_prod(N_{t},exp(-mx*dt)) * P + rt*rx
	
	  	These are the total number of fish (tagged + untagged)
		at large in the population.  N is used in the observation
		models to determine number of captures and recaptures is
		based on T.
		
		i = index for nrow
		j = index for size class
		t = index for year
		im= index for month/season/quarter (fraction of a year)
	*/
	int i,j,t,im;
	dvector newMarks(1,nx);	
	
	i=0; j=1;
	for(t = syr;t < nyr;t++)
	{
		for(im=1;im<=int(1/dt);im++)
		{
			i++;
			/* Recruitment if in period 1 */
			if(im==1 && t > syr)
			{
				N(i)+= mfexp(log_rt(t))*rx;
			}
			
			N(i+1)=elem_prod(N(i),mfexp(-mx*dt))*P;
			
			/*
			SJDM Might change this below.  As it currently stands the model
			is conditioned on the numbers of marks released, and these same
			data appear in the likelihood. Rather, the accounting of the Tagged
			and Untagged populations should be done in a separate routine where
			the predicted number of new marks released is based on the number of 
			untagged fish (N-T) captured at each sampling event.
			
			Set up accounting for predicted marks at large
			by adding observed marks to Mhat.
			T(i+1)= elem_prod(T(i),exp(-mx*dt))*P + M(i)
			
			
			

			newMarks.initialize();
			if(t==i_M(j,1) && im==i_M(j,2)) 
			{
				newMarks = M(j++);
				if(j>nrow) j=nrow;
			}
			
			T(i+1) = elem_prod(T(i),mfexp(-mx*dt)) * P + newMarks;
			*/
			//cout<<t<<"\t"<<i<<"\t"<<sum(N(i))<<"\t"<<sum(T(i))<<endl;
		} //end of im
	} // end of t
	//cout<<setprecision(3)<<N<<endl;
  }

FUNCTION calcNewMarksAndRecaptures
  {
	/*
	The accounting of the Tagged and Untagged populations where
	the predicted number of new marks released is based on the number of 
	untagged fish (N-T) captured at each sampling event. The precited 
	number of recaptures is based on the number of marks at large (T)
	
	j_min = index for minimum size of fish tagged (based on M data and flag(1))
	
	*/
	int i;
	
	Chat.initialize();
	Mhat.initialize();
	Rhat.initialize();
	T.initialize();
	dvar_vector zx(1,nx);
	dvar_vector ox(1,nx);
	for(i=1;i<=nrow;i++)
	{
		zx      = mx*dt;
		ox      = 1.0-mfexp(-zx);
		Chat(i) = elem_div(elem_prod(elem_prod(fi(i)*sx,ox),N(i)),zx);
		
		Rhat(i) = elem_div(elem_prod(elem_prod(fi(i)*sx,ox),T(i)),zx);
		
		Mhat(i)(j_min(i),nx) = Chat(i)(j_min(i),nx) - Rhat(i)(j_min(i),nx);
		
		if( i<nrow )
		{
			T(i+1) = elem_prod(T(i),mfexp(-zx)) * P + Mhat(i);
		}
	}
	
	//cout<<Mhat(2)<<endl;
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
		pvec(1) = dnorm(ddot_r_devs,0,0.4);
		pvec(2) = dnorm(bar_r_devs,0,0.4);
		pvec(3) = dnorm(log_bar_f,log(0.1108032),0.05);
	}
	else
	{
		pvec(1) = dnorm(first_difference(ddot_r_devs),0,0.4);
		pvec(2) = dnorm(bar_r_devs,0,2.5);
		pvec(3) = dnorm(log_bar_f,log(0.1108032),2.5);
	}
	
	//Likelihoods (fvec);
	fvec(1) = 10.*norm2(sqrt(C)-sqrt(Chat));
	fvec(2) = 10.*norm2(M-Mhat);
	fvec(3) = 10.*norm2(R-Rhat);
	//cout<<fvec<<endl;	
	
	/*
	PRIORS for estimated model parameters from the control file.
	*/
	dvar_vector priors(1,npar);
	dvariable ptmp; priors.initialize();
	for(i=1;i<=npar;i++)
	{
		if(active(theta(i)))
		{
			switch(theta_prior(i))
			{
			case 1:		//normal
				ptmp=dnorm(theta(i),theta_control(i,6),theta_control(i,7));
				break;
				
			case 2:		//lognormal CHANGED RF found an error in dlnorm prior. rev 116
				ptmp=dlnorm(theta(i),theta_control(i,6),theta_control(i,7));
				break;
				
			case 3:		//beta distribution (0-1 scale)
				double lb,ub;
				lb=theta_lbnd(i);
				ub=theta_ubnd(i);
				ptmp=dbeta((theta(i)-lb)/(ub-lb),theta_control(i,6),theta_control(i,7));
				break;
				
			case 4:		//gamma distribution
				ptmp=dgamma(theta(i),theta_control(i,6),theta_control(i,7));
				break;
				
			default:	//uniform density
				ptmp=-log(1./(theta_control(i,3)-theta_control(i,2)));
				break;
			}
			priors(i)=ptmp;	
		}
	}
	
	
	if(fpen > 0 ) cout<<"Fpen = "<<fpen<<endl;
	f = sum(fvec) + sum(pvec) + sum(priors);//	 + 1.e6*fpen;
  }



FUNCTION dvar_vector dgamma(const dvector& x, const dvariable& a, const dvariable& b)
  {
	//returns the gamma density with a & b as parameters
	RETURN_ARRAYS_INCREMENT();
	dvariable t1 = 1./(pow(b,a)*exp(gammln(a)));
	dvar_vector t2 = (a-1.)*log(x)-x/b;
	RETURN_ARRAYS_DECREMENT();
	return(t1*exp(t2));
  }

FUNCTION dvariable dgamma(const prevariable& x, const double& a, const double& b)
  {
	//returns the gamma density with a & b as parameters
	RETURN_ARRAYS_INCREMENT();
	dvariable t1 = 1./(pow(b,a)*exp(gammln(a)));
	dvariable t2 = (a-1.)*log(x)-x/b;
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


FUNCTION dvar_matrix calcLTM(dvector& x, const dvariable &linf, const dvariable &k, const dvariable &beta)
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
	int i,j,im;
	REPORT(f          );
	REPORT(log_ddot_r );
	REPORT(log_bar_r  );
	REPORT(log_bar_f  );
	REPORT(m_infty    );
	REPORT(l_infty    );
	REPORT(vbk        );
	REPORT(beta       );
	REPORT(mu_r       );
	REPORT(cv_r       );
	
	ivector yr(syr,nyr);
	yr.fill_seqadd(syr,1);
	REPORT(yr);
	REPORT(xmid);
	REPORT(mx);
	REPORT(sx);
	
	dvector Nt(syr,nyr);
	dvector Rt(syr,nyr);
	Nt.initialize();
	Rt.initialize();
	
	j = 1;
	for(i=syr;i<=nyr;i++)
	{
		if(i==syr) Rt(i) = value(exp(log_ddot_r+ddot_r_devs(nx)));
		else       Rt(i) = value(exp(log_rt(i)));
		for(im=1;im<=int(1/dt);im++)
		{
			Nt(i) += sum(value(N(j++)));
		}
	}
	REPORT(Nt);
	REPORT(Rt);
	
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
	#include <statsLib.h>
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

