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
//               Jul 10/12    8    LSMR simulation routine & likelihoods     //
//               Jul 11/12    8    R-scripts for HBC2011.dat                 //
//                                                                           //
//                                                                           //
// GIT REPOSITORY LOG:                                                       //
//  -BRANCHES: master     -> the finished product.                           //
//             ASMR       -> update using old ASMR model with presentation   //
//             developer  -> where code is developed to merge into master    //
//             |->newdata -> major code change to accomodate GILL & HOOP gear//
//                                                                           //
//                                                                           //
// TO DO LIST:                                                               //
//  -Model numbers-at-length on an annual time step, then in the observation //
//   submodel, grow and survive numbers-at-length upto the time step samples //
//   were collected.  This will require transition matrix for dt and annual  //
//   transition matrix.                                                      //
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
	!! cout<<" TOP OF DATA_SECTION "<<endl;
	
	// Read input data from data file //
	init_int syr;   	// first year
	init_int nyr;   	// last year
	init_number dt; 	// time-step
	init_int ngear; 	// number of gears
	init_int nbin;		// number of length intervals
	init_vector xbin(1,nbin);  //length-intervals (not mid-points)
	
	int nx;			//number of length bins (nbin-1)
	int nr;			//number of rows in predicted arrays
	!! nr = int((nyr-syr+1)/dt);
	!! nx = nbin;
	vector xmid(1,nbin);
	!! xmid = xbin + 0.5*(xbin(2)-xbin(1));
		
	// Array dimensions //
	init_imatrix dim_array(1,ngear,1,2);
	ivector irow(1,ngear);
	ivector ncol(1,ngear);
	ivector jcol(1,ngear);
	LOC_CALCS
		irow = column(dim_array,1);
		ncol = column(dim_array,2);
		jcol = ncol - 1;
	END_CALCS
	// Read in effort data (number of sets) //
	init_matrix effort(1,ngear,1,irow);
	vector mean_effort(1,ngear);
	LOC_CALCS
		/* Calculate mean effort for each gear, ignore 0*/
		int i,k,n;
		for(k=1;k<=ngear;k++)
		{
			n=0;
			for(i=1;i<=irow(k);i++)
			{
				if(effort(k,i)>0)
				{
					mean_effort(k) += effort(k,i);
					n++;
				}
			}
			mean_effort(k) /= n;
		}
	END_CALCS
	
	// Read in Capture, Mark, and Recapture arrays //
	init_3darray i_C(1,ngear,1,irow,1,ncol);
	init_3darray i_M(1,ngear,1,irow,1,ncol);
	init_3darray i_R(1,ngear,1,irow,1,ncol);	

	3darray C(1,ngear,1,irow,1,jcol);
	3darray M(1,ngear,1,irow,1,jcol);
	3darray R(1,ngear,1,irow,1,jcol);	

	
	init_int eof;
	!! if(eof!=999){cout<<"Error reading data\n eof = "<<eof<<endl; exit(1);}
	!! cout<<" - END OF READING DATA"<<endl;
	
	
	int fi_count;
	// colsums of Catch-at-length //
	matrix ct(1,ngear,1,irow); 
	
	LOC_CALCS
		/* number of capture probability deviates */
		fi_count=0;
		for(k=1;k<=ngear;k++)
		{
			for(i=1;i<=irow(k);i++)
			{
				if( effort(k,i) > 0 ) fi_count++;
				C(k)(i) = i_C(k)(i)(2,ncol(k)).shift(1);
				M(k)(i) = i_M(k)(i)(2,ncol(k)).shift(1);
				R(k)(i) = i_R(k)(i)(2,ncol(k)).shift(1);
				ct(k,i) = sum( C(k)(i) );
			}
		}
	END_CALCS
	
	
	// OPEN CONTROL FILE //
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
	// 1) verbose 
	// 2) minimum size of fish for tagging.
	// 3) std of catch residuals in 1st phase
	// 4) std of catch residuals in last phase
	
	
	// index for minimum size of tagged fish //
	imatrix min_tag_j(1,ngear,1,irow);
	LOC_CALCS
		int j;
		for(k=1;k<=ngear;k++)
		{
			for(i=1;i<=irow(k);i++)
			{
				for(j=1;j<=jcol(k);j++)
				{
					if( M(k)(i,j)>0 || xbin(j)>flag(2) )
					{
						min_tag_j(k,i) = j;
						break;
					}
				}	
			}
		}
	END_CALCS
	
	
	// Variables for Simulated Data
	vector true_Nt(syr,nyr);
	vector true_Rt(syr,nyr);
	vector true_Tt(syr,nyr);
	matrix true_fi(1,ngear,1,irow);
	
	!! cout<< " END OF DATA_SECTION \n"<<endl;
	
	
PARAMETER_SECTION
	init_bounded_number_vector theta(1,npar,theta_lbnd,theta_ubnd,theta_phz);
	number log_ddot_r;
	number log_bar_r;
	number m_infty;

	//Variables for growth.
	number l_infty;
	number vbk;
	number beta;
	
	//Variables for size distribution of new recruits
	number mu_r;
	number cv_r
	
	//Mean fishing mortality rates
	init_vector log_bar_f(1,ngear,2);
	
	//Overdispersion
	init_vector log_tau(1,ngear,5);
	
	//Selectivity parameters
	init_vector log_lx(1,ngear,3);
	init_vector log_gx(1,ngear,3);
	
	init_bounded_dev_vector ddot_r_devs(1,nx,-15,15,2);
	init_bounded_dev_vector bar_r_devs(syr+1,nyr,-15,15,2);
	
	
	//TODO Fix this so there is a dev vector for each gear, otherwise biased estimates of log_bar_f
	
	init_bounded_dev_vector bar_f_devs(1,fi_count,-5.0,5.0,2);
	
	
	
	objective_function_value f;
	

	
	number m_linf;
	number fpen;
	vector tau(1,ngear);	// over-dispersion parameters >1.0
	vector lx(1,ngear);		// length at 50% selectivity
	vector gx(1,ngear);		// std in length at 50% selectivity
	vector qk(1,ngear);		// catchability of gear k
	vector mx(1,nx);		// Mortality rate at length xmid
	vector rx(1,nx);		// size pdf for new recruits
	
	
	vector log_rt(syr+1,nyr);
	matrix fi(1,ngear,1,irow);// capture probability in period i.
	matrix sx(1,ngear,1,jcol);	// Selectivity at length xmid
	matrix N(syr,nyr,1,nx);		// Numbers(time step, length bins)
	matrix T(syr,nyr,1,nx);		// Marks-at-large (time step, length bins)
	matrix A(1,nx,1,nx);		// Size-transitin matrix (annual step)
	//matrix P(1,nx,1,nx);		// Size-Transition Matrix for step dt
	
	// Predicted observations //
	matrix hat_ct(1,ngear,1,irow);			// Predicted total catch
	matrix delta(1,ngear,1,irow);			// residuals in total catch
	3darray Chat(1,ngear,1,irow,1,jcol);	// Predicted catch-at-length
	3darray Mhat(1,ngear,1,irow,1,jcol);	// Predicted new marks-at-length
	3darray Rhat(1,ngear,1,irow,1,jcol);	// Predicted recaptures-at-length
	
	
	
INITIALIZATION_SECTION
	theta     theta_ival;
	log_lx    2.3;
	log_gx    0.2;
	log_bar_f -2.3;
	log_tau   1.1;

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
	if( flag(1) ) cout<<"\n TOP OF PROCEDURE_SECTION "<<endl;
	fpen.initialize();
	initParameters();  
	calcSurvivalAtLength(); 
	calcSizeTransitionMatrix(); 
	initializeModel();
	calcCaptureProbability();
	calcNumbersAtLength();
	calcSelectivityAtLength();
	calcObservations();
	calc_objective_function();
	if( flag(1) ) cout<<"\n END OF PROCEDURE_SECTION "<<endl;

//	
FUNCTION void runSimulationModel(const int& seed)
  {
	int i,j,k,im;
	random_number_generator rng(seed);
	
	dvector tmp_ddot_r_devs = value(ddot_r_devs);
	dvector tmp_bar_r_devs  = value(bar_r_devs);
	dvector tmp_bar_f_devs  = value(bar_f_devs);
	
	tmp_ddot_r_devs.fill_randn(rng);
	tmp_bar_r_devs.fill_randn(rng);
	tmp_bar_f_devs.fill_randn(rng);
	
	ddot_r_devs = 0.6*tmp_ddot_r_devs;
	bar_r_devs  = 0.6*tmp_bar_r_devs;
	bar_f_devs  = 0.2*tmp_bar_f_devs;
	/* Capture probabilities */
	
	
	initParameters();
	calcSurvivalAtLength();
	calcSizeTransitionMatrix();
	initializeModel();
	calcCaptureProbability();
	calcNumbersAtLength();
	calcSelectivityAtLength();
	calcObservations();
	
	
	/* Overwrite observations and draw from multinomial distribution */
	C=value(Chat);
	M=value(Mhat);
	R=value(Rhat);
	
	for(k=1;k<=ngear;k++)
	{
		for(i=1;i<=irow(k);i++)
		{
			//ct(k,i) = sum(C(k)(i));
			//C(k)(i)  = rmultinom(rng,int(ct(k,i)),C(k)(i));
			//M(k)(i)  = rmultinom(rng,int(sum(M(k)(i))),M(k)(i));
			//R(k)(i)  = rmultinom(rng,int(sum(R(k)(i))),R(k)(i));
			if(effort(k,i)>0)
			{
				for(j=1;j<=nx;j++)
				{
					C(k)(i)(j) = randnegbinomial(1.e-5+C(k)(i)(j),2.0,rng);
					M(k)(i)(j) = randnegbinomial(1.e-5+M(k)(i)(j),2.0,rng);
					R(k)(i)(j) = randnegbinomial(1.e-5+R(k)(i)(j),2.0,rng);
				}
			}
			
			i_C(k)(i)(2,ncol(k)) = C(k)(i)(1,nx).shift(2);
			i_M(k)(i)(2,ncol(k)) = M(k)(i)(1,nx).shift(2);
			i_R(k)(i)(2,ncol(k)) = R(k)(i)(1,nx).shift(2);
			
		}		
	}
	
	/*Save true data*/
	true_Nt.initialize();
	true_Tt.initialize();
	true_Rt.initialize();
	true_fi = value(fi);
	true_Nt = value(rowsum(N));
	true_Tt = value(rowsum(T));
	for(i=syr;i<=nyr;i++)
	{
		if(i==syr) true_Rt(i) = value(mfexp(log_ddot_r+ddot_r_devs(nx)));
		else       true_Rt(i) = value(mfexp(log_rt(i)));
	}
	
  }
//
FUNCTION initParameters
  {
	int k;
	/* Leading parameters */
	log_ddot_r = theta(1);
	log_bar_r  = theta(2);	
	m_infty    = theta(3);
	l_infty    = theta(4);
	vbk        = theta(5);
	beta       = theta(6);
	mu_r       = theta(7);
	cv_r       = theta(8);
	//for(k=1;k<=ngear;k++)
	//{
	//	log_bar_f(k)  = theta(8+k);
	//}
	
	
	/* Selex parameters */
	lx         = mfexp(log_lx);
	gx         = mfexp(log_gx);
	
	/* Catchability */
	qk         = elem_div(mfexp(log_bar_f),mean_effort);
	
	
  }
//
FUNCTION calcSizeTransitionMatrix
  {
	/*
	This function calls the necessary routines to compute the Size Transition Matrix (P)
	*/
	A = calcLTM(xmid,l_infty,vbk,beta);
	//P = calcLTM(xmid,l_infty,vbk,beta/dt);
  }
//
FUNCTION initializeModel
  {
	int i,ii,k,ik;
	//Set up the initial states.
	N.initialize();
	T.initialize();
	
	
	/* Initialize numbers at length at the first time step. */
	dvar_vector init_r(1,nx);
	dvar_matrix phi_X(1,nx,1,nx);
	dvariable a=1/(cv_r*cv_r);  //a=1/cv^2
	dvariable b=mu_r/a;		  //b=mu/a
	rx  = dgamma(xmid,a,b);
	rx /= sum(rx);
	
	/* Calculate per-recruit survivorship at length */
	phi_X(1) = rx;
	ii = 0;
	for(i=2;i<=nx;i++)
	{
		phi_X(i) = elem_prod(phi_X(i-1),mfexp(-mx)) * A;
		if( i==nx )
		{
			phi_X(i) = phi_X(i) + elem_div(phi_X(i),1.-mfexp(-mx));
		}
	}
	
	/* Initial numbers at length */
	init_r = mfexp(log_ddot_r + ddot_r_devs);
	N(syr)   = init_r * phi_X;
	
	/* Annual recruitment */
	log_rt = log_bar_r + bar_r_devs;
  }
//
FUNCTION calcCaptureProbability		
  {
	/* Capture probability at each time step. */
	int i,k,ik;
	ik = 1;
	fi.initialize();
	for(k=1;k<=ngear;k++)
	{
		for(i=1;i<=irow(k);i++)
		{
			if( effort(k,i)>0 )
			{
				//fi(k,i)  = mfexp(log_bar_f(k) + bar_f_devs(ik++));
				fi(k,i) = qk(k)*effort(k,i)*mfexp(bar_f_devs(ik++));
			}
		}
	}
	
  }
//
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
//
FUNCTION calcSelectivityAtLength
  {
	/*
	This function calculates the length-specific selectivity.
	The parametric option is a logistic curve (plogis(x,lx,gx))
	*/
	int k;
	sx.initialize();
	for(k=1;k<=ngear;k++)
	{
		//sx(k) = 1./(1+mfexp(-(xmid-lx(k))/gx(k)));
		sx(k) = plogis(xmid,lx(k),gx(k));
	}
	
  }
//
FUNCTION calcNumbersAtLength
  {
	/*	This function updates the numbers at length
		at the start of each year as a function of the 
		numbers at length times the survival rate * size transition
		and add new recuits-at-length.
		
	  	N_{t+1} = elem_prod(N_{t},mfexp(-mx*dt)) * P + rt*rx
	
	  	These are the total number of fish (tagged + untagged)
		at large in the population.  N is used in the observation
		models to determine number of captures and recaptures is
		based on T.
		
		t = index for year
	*/
	int t;
	dvariable rt;
	for(t=syr;t<nyr;t++)
	{
		rt     = mfexp(log_rt(t+1));
		N(t+1) = elem_prod(N(t),mfexp(-mx)) * A + rt*rx;
	}
	if( flag(1)==2 ) cout<<"Nt\n"<<rowsum(N)<<endl;
  }
//
FUNCTION calcObservations
  {
	/*
		Calculate the predicted total catch-at-length (Chat)
		Calculate the predicted total recaptures-at-length (Rhat)
		Calculate the predicted total new markes released-at-length (Mhat)
		Calculate and update the number of marks at large.
		
		t   = index for year
		its = index for time step
		k   = index for gear
	*/
	
	int t,its,k,i,lb;
	dvar_vector Ntmp(1,nx);
	dvar_vector Ttmp(1,nx);
	dvar_vector Utmp(1,nx);
	dvar_vector Mtmp(1,nx);
	dvar_vector zx(1,nx);
	dvar_vector ox(1,nx);
	dvar_vector ux(1,nx);
	
	i=0;
	zx = mx*dt;
	ox = elem_div(1.0-mfexp(-zx),zx);
	
	for(t=syr;t<=nyr;t++)
	{
		Ntmp = N(t);
		Ttmp = T(t);
		Utmp = posfun(Ntmp - Ttmp,0.01,fpen);
		Mtmp.initialize();
		i++;
		for(k=1;k<=ngear;k++)
		{
			if( effort(k,i)>0 )
			{
				lb           = min_tag_j(k,i);
				ux           = 1.0 - mfexp(-fi(k,i)*sx(k));
				Chat(k)(i)   = elem_prod(ux,elem_prod(Ntmp,ox));
				Mhat(k)(i)   = elem_prod(ux,elem_prod(Utmp,ox));
				Rhat(k)(i)   = elem_prod(ux,elem_prod(Ttmp,ox));
				Mtmp(lb,nx) += Mhat(k)(i)(lb,nx);
			}
		}
		
		/* Survive and grow tags-at-large and add new tags */
		if( t < nyr )
		{
			T(t+1) = elem_prod(T(t),mfexp(-mx*dt))*A + Mtmp;
		}
	}
	
	if( flag(1)==2 ) cout<<"Tt\n"<<rowsum(T)<<endl;
  }
//
FUNCTION calc_objective_function;
  {
	int i,j,k;
	/* PENALTIES TO ENSURE REGULAR SOLUTION */
	dvar_vector pvec(1,3);
	pvec.initialize();
	if(!last_phase())
	{
		pvec(1) = dnorm(first_difference(ddot_r_devs),0,0.4);
		pvec(2) = dnorm(bar_r_devs,0,0.4);
		pvec(3) = dnorm(log_bar_f,log(0.1108032),0.05);
		//pvec(3) = 1.e5 * square(log_bar_f - log(0.110));
	}
	else
	{
		pvec(1) = dnorm(first_difference(ddot_r_devs),0,0.4);
		pvec(2) = dnorm(bar_r_devs,0,2.5);
		pvec(3) = dnorm(log_bar_f,log(0.1108032),2.5);
	}
	if( flag(1) ) cout<<"Average fi = "<<mfexp(log_bar_f)<<endl;
	/* LIKELIHOODS */
	/*
		fvec(1) = likelihood of the total catch-at-length.
		fvec(2) = likelihood of the total marks-at-length.
		fvec(3) = likelihood of the total recap-at-length.
	*/
	
	dvar_vector fvec(1,4);
	fvec.initialize();
	double tiny = 1.e-10;
	tau = mfexp(log_tau)+1.01;
	
	for(k=1;k<=ngear;k++)
	{
		for(i=1;i<=irow(k);i++)
		{
			if( effort(k,i)>0 )
			{
				for(j=1;j<=nx;j++)
				{
					fvec(1) -= log_negbinomial_density(C(k)(i)(j),Chat(k)(i)(j)+tiny,tau(k));
					fvec(2) -= log_negbinomial_density(M(k)(i)(j),Mhat(k)(i)(j)+tiny,tau(k));
					fvec(3) -= log_negbinomial_density(R(k)(i)(j),Rhat(k)(i)(j)+tiny,tau(k));
				}
			}
		}
	}
	//if(flag(5))
	//{
	//	for(k=1;k<=ngear;k++)
	//	{
	//		fvec(2) += 10.* norm2(sqrt(C(k))-sqrt(Chat(k)));
	//		fvec(3) += 10.* norm2(sqrt(M(k))-sqrt(Mhat(k)));
	//		fvec(4) += 10.* norm2(R(k)-Rhat(k));
	//	}
	//}
	/*else
		{
			for(k=1;k<=ngear;k++)
			{
				for(i=1;i<=irow(k);i++)
				{
					if(effort(k,i)>0)
					{
						tau = 400;
						if( tau >0 )
						{
							fvec(2) += dmultifan( C(k)(i),Chat(k)(i),tau );
						}
				
						if( tau>0 )
						{
							fvec(3) += dmultifan( M(k)(i),Mhat(k)(i),tau );
						}
				
						if( tau>0 )
						{
							fvec(4) += dmultifan( R(k)(i),Rhat(k)(i),tau );
						}
					}
				}
			}
		}*/
	
	/*
	PRIORS for estimated model parameters from the control file.
	*/
	
	dvar_vector priors(1,npar);
	priors.initialize();
	dvariable ptmp; 
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
	
	if( flag(1) ) cout<<"Fvec\t"<<setprecision(4)<<fvec<<endl;
	if(fpen > 0) cout<<"Fpen = "<<fpen<<endl;
	f = sum(fvec) + sum(pvec) + sum(priors) + 1.e5*fpen;
  }
//
FUNCTION dvar_vector dgamma(const dvector& x, const dvariable& a, const dvariable& b)
  {
	//returns the gamma density with a & b as parameters
	RETURN_ARRAYS_INCREMENT();
	dvariable t1 = 1./(pow(b,a)*mfexp(gammln(a)));
	dvar_vector t2 = (a-1.)*log(x)-x/b;
	RETURN_ARRAYS_DECREMENT();
	return(t1*mfexp(t2));
  }
//
FUNCTION dvariable dgamma(const prevariable& x, const double& a, const double& b)
  {
	//returns the gamma density with a & b as parameters
	RETURN_ARRAYS_INCREMENT();
	dvariable t1 = 1./(pow(b,a)*mfexp(gammln(a)));
	dvariable t2 = (a-1.)*log(x)-x/b;
	RETURN_ARRAYS_DECREMENT();
	return(t1*mfexp(t2));
  }

FUNCTION dvar_vector posfun(const dvar_vector& x, const double& eps, dvariable& pen)
  {
	int i;
	dvar_vector xp(x.indexmin(),x.indexmax());
	for(i=x.indexmin();i<=x.indexmax();i++)
	{
		if(x(i)>=eps)
		{
			xp(i) = x(i);
		}
		else
		{
			pen += 0.01*square(x(i)-eps);
			xp(i) = eps/(2.-x(i)/eps);
		}
	}
	return(xp);
  }

FUNCTION ivector match(const ivector& x, const ivector& table)
  {
	//returns a vector of positions of first matches of x in table.
	int i,j;
	ivector pos(x.indexmin(),x.indexmax());
	for(i=x.indexmin();i<=x.indexmax();i++)
	{
		for(j=table.indexmin();j<=table.indexmax();j++)
		{
			if(x(i) == table(j) )
			{
				pos(i) = j;
				break;
			}
		}
	}
	return(pos);
  }
//
FUNCTION dvar_matrix calcLTM(dvector& x, const dvariable &linf, const dvariable &k, const dvariable &beta)
  {
	/*This function computes the length transition matrix.*/
	/*
	- x is the mid points of the length interval vector.
	
	- cumd_gamma(x,a) is the same as Igamma(a,x) in R.
	- If length interval > linf, then assume now further growth.
	- Note the use of posfun to ensure differentiable.
	*/
	RETURN_ARRAYS_INCREMENT();
	int i,j;
	double dx;
	double bw = 0.5*(x(2)-x(1));
	int n = size_count(x);     //number of length intervals
	dvar_vector alpha(1,n);
	dvar_matrix P(1,n,1,n);
	P.initialize();
	
	//Growth increment
	dvar_vector dl(1,n);
	for(j=1;j<=n;j++)
	{
		dl(j) = log(mfexp( (linf-x(j))*(1.-mfexp(-k)) )+1.0);
	}
	alpha = dl/beta;
	
	/* Size transition probability */
	for(i=1;i<=n;i++)
	{
		dvar_vector t1(i,n);
		for(j=i;j<=n;j++)
		{
			dx = x(j)-x(i);
			P(i)(j) = cumd_gamma(dx+bw,alpha(i)) - cumd_gamma(dx-bw,alpha(i));
			//cout<<j<<" "<<x(i)<<" "<<x(j)<<" "<<P(i)(j)<<endl;
		}
		//P(i)(i,n-1) = first_difference(t1);
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
	REPORT(m_infty    );
	REPORT(l_infty    );
	REPORT(vbk        );
	REPORT(beta       );
	REPORT(mu_r       );
	REPORT(cv_r       );
	REPORT(log_bar_f  );
	REPORT(tau        );
	
	ivector yr(syr,nyr);
	yr.fill_seqadd(syr,1);
	REPORT(yr);
	REPORT(ngear);
	REPORT(irow);
	REPORT(jcol);
	REPORT(xmid);
	REPORT(mx);
	REPORT(sx);
	
	REPORT(delta);
	
	dvector Nt = value(rowsum(N));
	dvector Tt = value(rowsum(T));
	dvector Rt(syr,nyr);
	for(i=syr;i<=nyr;i++)
	{
		if(i==syr) Rt(i) = value(mfexp(log_ddot_r+ddot_r_devs(nx)));
		else       Rt(i) = value(mfexp(log_rt(i)));
	}
	REPORT(Nt);
	REPORT(Tt);
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
	
	if(SimFlag)
	{
		REPORT(true_Nt);
		REPORT(true_Rt);
		REPORT(true_Tt);
		REPORT(true_fi);
	}

//
TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 
//
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
	
	dvariable dmultifan(const dvector& o,const dvar_vector& p,const double& s)
	{
		/*
		o is the observed numbers at length
		p is the predicted numbers at length
		s is the minimum sample size
		*/
		RETURN_ARRAYS_INCREMENT();		
		int lb     = o.indexmin();
		int nb     = o.indexmax();
		int I      = (nb-lb)+1;
		double n   = sum(o);
		if(min(n,s)<=0)
		{
			RETURN_ARRAYS_DECREMENT();
			return(0);
		} 
		double tau = 1./min(n,s);

		dvariable nll;
		dvector O     = o/sum(o);
		dvar_vector P = p/sum(p);
		dvar_vector epsilon(lb,nb);
		epsilon       = elem_prod(1.-P,P);

		dvariable T1,T2,T3;
		T1 = -0.5 * sum(log( 2.*M_PI*(epsilon+0.1/I) ));
		T2 = -0.5 * I * log(tau);
		T3 = sum( log(exp(-1.0 * elem_div(square(O-P),2.0*tau*(epsilon+0.1/I)) )+0.01) );
		
		nll = -1.0*(T1 + T2 + T3);
		
		//cout<<"dmultifan like = "<<nll<<"\tn = "<<1./tau<<"\tsum(P) = "<<sum(O)<<endl;
		//if(n==112) cout<<p<<endl;
		RETURN_ARRAYS_DECREMENT();
		return nll;
	}
	
	dvariable dnbinom(const dvector& x, const dvar_vector& mu, const prevariable& k)
	{
		//the observed counts are in x
		//mu is the predicted mean
		//k is the overdispersion parameter
		if (value(k)<0.0)
		{
			cerr<<"k is <=0.0 in dnbinom()";
			return(0.0);
		}
		RETURN_ARRAYS_INCREMENT();
		int i,imin,imax;
		imin=x.indexmin();
		imax=x.indexmax();
		dvariable loglike = 0.;

		for(i = imin; i<=imax; i++)
		{
			cout<<"mu "<<mu(i)<<endl;
			loglike += gammln(k+x(i))-gammln(k)-gammln(x(i)+1)+k*log(k)-k*log(mu(i)+k)+x(i)*log(mu(i))-x(i)*log(mu(i)+k);
		}
		RETURN_ARRAYS_DECREMENT();
		return(-loglike);
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

