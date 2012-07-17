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
//  -Change indexing of data such that input observations have the same      //
//   indexing as the model dimensions.                                       //
//                                                                           //
//                                                                           //
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
		jcol = ncol - 2;
	END_CALCS
	// Read in effort data (number of sets) //
	init_matrix effort(1,ngear,1,irow);
	
	// Read in Capture, Mark, and Recapture arrays //
	init_3darray i_C(1,ngear,1,irow,1,ncol);
	init_3darray i_M(1,ngear,1,irow,1,ncol);
	init_3darray i_R(1,ngear,1,irow,1,ncol);	

	3darray C(1,ngear,1,irow,1,jcol);
	3darray M(1,ngear,1,irow,1,jcol);
	3darray R(1,ngear,1,irow,1,jcol);	

	
	init_int eof;
	!! if(eof!=999){cout<<"Error reading data\n eof = "<<eof<<endl; exit(1);}
	
	
	
	int fi_count;
	// colsums of Catch-at-length //
	matrix ct(1,ngear,1,irow); 
	
	LOC_CALCS
		/* number of capture probability deviates */
		int i,k;
		fi_count=0;
		for(k=1;k<=ngear;k++)
		{
			for(i=1;i<=irow(k);i++)
			{
				if( effort(k,i) > 0 ) fi_count++;
				C(k)(i) = i_C(k)(i)(3,ncol(k)).shift(1);
				M(k)(i) = i_M(k)(i)(3,ncol(k)).shift(1);
				R(k)(i) = i_R(k)(i)(3,ncol(k)).shift(1);
				ct(k,i) = sum( C(k)(i) );
			}
		}
	END_CALCS
	
	
	
	
	
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
	// 2) std of catch residuals in 1st phase
	// 3) std of catch residuals in last phase
	
	
	// Pre-processing some of the data.
	imatrix j_min(1,ngear,1,irow);
	//3darray C(1,ngear,1,irow_C,1,jcol_MR);
	//3darray M(1,ngear,1,irow_MR,1,jcol_MR);
	//3darray R(1,ngear,1,irow_MR,1,jcol_MR);
	//LOC_CALCS
	//	int i,j,k;
	//	
	//	for(k=1;k<=ngear;k++)
	//	{
	//		for(i=1;i<=irow_C(k);i++)
	//		{
	//			C(k)(i) = i_C(k)(i)(2,jcol_C(k)).shift(1);
	//			ct(k,i) = sum(C(k)(i));
	//		}
	//		
	//		for(i=1;i<=irow_MR(k);i++)
	//		{
	//			M(k)(i) = i_M(k)(i)(2,jcol_C(k)).shift(1);
	//			R(k)(i) = i_R(k)(i)(2,jcol_C(k)).shift(1);
	//			//set up index for minimum size for tagging fish.
	//			for(j=1;j<=jcol_MR(k);j++)
	//			{
	//				if( M(k)(i,j)>0 || xbin(j) > flag(1) )
	//				{
	//					j_min(k)(i)=j;
	//					break;
	//				}
	//			}
	//		}
	//	}
	//END_CALCS
	
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
	init_vector log_lx(1,ngear,2);
	init_vector log_gx(1,ngear,2);
	
	init_bounded_dev_vector ddot_r_devs(1,nx,-15,15,2);
	init_bounded_dev_vector bar_r_devs(syr+1,nyr,-15,15,3);
	init_bounded_dev_vector bar_f_devs(1,fi_count,-5.0,5.0,2);
	
	
	
	objective_function_value f;
	

	number dl_cv;
	number fpen;
	
	number m_linf;
	vector lx(1,ngear);		// length at 50% selectivity
	vector gx(1,ngear);		// std in length at 50% selectivity
	vector mx(1,nx);		// Mortality rate at length xmid
	vector rx(1,nx);		// size pdf for new recruits
	
	vector log_rt(syr+1,nyr);
	matrix fi(1,ngear,1,irow);// capture probability in period i.
	matrix sx(1,ngear,1,jcol);	// Selectivity at length xmid
	matrix N(1,nr,1,nx);		// Numbers(time step, length bins)
	matrix T(1,nr,1,nx);		// Marks-at-large (time step, length bins)
	matrix P(1,nx,1,nx);		// Size-Transition Matrix
	
	// Predicted observations //
	matrix hat_ct(1,ngear,1,irow);			// Predicted total catch
	3darray Chat(1,ngear,1,irow,1,jcol);	// Predicted catch-at-length
	3darray Mhat(1,ngear,1,irow,1,jcol);	// Predicted new marks-at-length
	3darray Rhat(1,ngear,1,irow,1,jcol);	// Predicted recaptures-at-length
	
	
	
INITIALIZATION_SECTION
	theta theta_ival;
	log_lx 2.3;

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
		//runSimulationModel(rseed);
		
	}


PROCEDURE_SECTION
	cout<<"\n TOP OF PROCEDURE_SECTION "<<endl;
	initParameters();  
	calcSurvivalAtLength(); 
	calcSizeTransitionMatrix(); 
	initializeModel();cout<<"OK HERE"<<endl;
	calcNumbersAtLength();
	calcSelectivityAtLength();
	calcNewMarksAndRecaptures();
	
	cout<<"\n END OF PROCEDURE_SECTION "<<endl;
//		
//	//calc_catch_at_length();
//	calc_objective_function();
//	
//FUNCTION void runSimulationModel(const int& seed)
//  {
//	int i,j,im;
//	random_number_generator rng(seed);
//	
//	dvector tmp_ddot_r_devs = value(ddot_r_devs);
//	dvector tmp_bar_r_devs  = value(bar_r_devs);
//	tmp_ddot_r_devs.fill_randn(rng);
//	tmp_bar_r_devs.fill_randn(rng);
//	
//	ddot_r_devs = 0.6*tmp_ddot_r_devs;
//	bar_r_devs  = 0.6*tmp_bar_r_devs;
//	
//	/* Capture probabilities */
//	//Et <- h1+(h2-h1)*cos(0.5+(sample.yrs-t)/(T-t)*pr*pi)^2
//	double h1 = -1.;
//	double h2 =  1.;
//	for(i=1;i<=nrow;i++)
//	{
//		bar_f_devs(i) = h1 + (h2-h1)*square(cos((i-1.)/(nrow-1.)*0.85*PI));
//	}
//	
//	initParameters();
//	calcSurvivalAtLength();
//	calcSizeTransitionMatrix();
//	initializeModel();
//	calcNumbersAtLength();
//	calcSelectivityAtLength();
//	calcNewMarksAndRecaptures();
//	
//	/* Overwrite observations and draw from multinomial distribution */
//	C=value(Chat);
//	M=value(Mhat);
//	R=value(Rhat);
//	
//	for(i=1;i<=nrow;i++)
//	{
//		ct(i) = sum(C(i));
//		C(i)  = rmultinom(rng,int(ct(i)),C(i));
//		M(i)  = rmultinom(rng,int(sum(M(i))),M(i));
//		R(i)  = rmultinom(rng,int(sum(R(i))),R(i));
//		
//		i_C(i)(3,ncol) = C(i)(1,nx).shift(3);
//		i_M(i)(3,ncol) = M(i)(1,nx).shift(3);
//		i_R(i)(3,ncol) = R(i)(1,nx).shift(3);
//	}
//	
//	/*Save true data*/
//	true_Nt.initialize();
//	true_Rt.initialize();
//	true_fi = value(fi);
//	j = 1;
//	for(i=syr;i<=nyr;i++)
//	{
//		if(i==syr) true_Rt(i) = value(exp(log_ddot_r+ddot_r_devs(nx)));
//		else       true_Rt(i) = value(exp(log_rt(i)));
//		for(im=1;im<=int(1/dt);im++)
//		{
//			true_Nt(i) += sum(value(N(j++)));
//		}
//	}
//	
//  }
//
//
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
//
FUNCTION calcSizeTransitionMatrix
  {
	/*
	This function calls the necessary routines to compute the Size Transition Matrix (P)
	*/
	
	P = calcLTM(xmid,l_infty,vbk,beta);
  }
//
FUNCTION initializeModel
  {
	int i,k,ik;
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
	ik = 1;
	for(k=1;k<=ngear;k++)
	{
		for(i=1;i<=irow(k);i++)
		{
			if( effort(k,i)>0 )
			{
				fi(k,i)  = mfexp(log_bar_f + bar_f_devs(ik++));
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
	
	//cout<<sx<<endl;
  }
//
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
//
FUNCTION calcNewMarksAndRecaptures
  {
	/*
	The accounting of the Tagged and Untagged populations where
	the predicted number of new marks released is based on the number of 
	untagged fish (N-T) captured at each sampling event. The precited 
	number of recaptures is based on the number of marks at large (T)
	
	j_min = index for minimum size of fish tagged (based on M data and flag(1))
	
	PSUEDOCODE:
		-1) loop over years and within each time step determine if a trip occured
		-2) get corresponding row index for C, M, and R for gear k
		-3) get corresponding col index for C, M, and R for gear k
	
	*/
	
	
	cout<<"calcNewMarksAndRecaptures"<<endl;
	//
	//int i,j,k,t,im,ik,jj;
	//ivector iyr(1,ngear);
	//
	///* Match index for size bins in the model with those in the data. */
	//imatrix ix(1,ngear,1,jcol_MR);
	//for(k=1;k<=ngear;k++)
	//{
	//	ix(k) = match(x_MR(k),xbin);
	//}
	//
	//dvariable ftmp;
	//Chat.initialize();
	//Mhat.initialize();
	//Rhat.initialize();
	//T.initialize();
	//dvar_vector zx(1,nx);
	//dvar_vector ox(1,nx);
	//dvar_vector newMarks(1,nx);
	//
	//i   = 0;
	//ik  = 0;
	//iyr = 1;
	//zx  = mx*dt;
	//ox  = 1.0-mfexp(-zx);
	//for(t=syr;t<=nyr;t++)
	//{
	//	for(im=1;im<=int(1/dt);im++)
	//	{
	//		i++;	// index for row of N and T
	//		newMarks = 0;
	//		for(k=1;k<=ngear;k++)
	//		{
	//			if( i_C(k)(iyr(k))(1)==t && iyr(k)<=irow_C(k) )
	//			{
	//				ftmp = fi(k,iyr(k));
	//				
	//				dvar_vector t1 = elem_prod( ftmp*sx(k),ox(ix(k)) );
	//				dvar_vector t2 = elem_div( N(i)(ix(k)),zx(ix(k)) );
	//				dvar_vector t3 = elem_div( T(i)(ix(k)),zx(ix(k)) );
	//				
	//				Chat(k)(iyr(k)) = elem_prod(t1,t2);
	//				Rhat(k)(iyr(k)) = elem_prod(t1,t3);
	//				
	//				//jx              = j_min(k,iyr(k));
	//				//Mhat(k)(iyr(k))(jx,nx) = 
	//				
	//				//T(i) = 1.0;//T(i)(ix(k)) + Mhat(k)(iyr(k));
	//				
	//				//cout<<T(i)<<endl;
	//				if(iyr(k) < irow_C(k)) iyr(k) ++;
	//			}
	//			
	//		}
	//		if( i<nr )
	//		{
	//			T(i+1) = elem_prod(T(i),mfexp(-zx)) * P + newMarks;
	//		}
	//	}
	//}
	//cout<<Chat(1)<<endl;
	/*for(i=1;i<=nrow;i++)
		{
			zx      = mx*dt;
			ox      = 1.0-mfexp(-zx);
			Chat(i) = elem_div(elem_prod(elem_prod(fi(i)*sx,ox),N(i)),zx);
			//cout<<fi<<endl;
			hat_ct(i) = sum(Chat(i));
			
			Rhat(i) = elem_div(elem_prod(elem_prod(fi(i)*sx,ox),T(i)),zx);
			
			Mhat(i)(j_min(i),nx) = Chat(i)(j_min(i),nx) - Rhat(i)(j_min(i),nx);
			
			if( i<nrow )
			{
				T(i+1) = elem_prod(T(i),mfexp(-zx)) * P + Mhat(i);
			}
		}*/
	
	//cout<<Chat(2)<<endl;
  }

//
//FUNCTION calc_objective_function;
//  {
//	int i;
//	dvar_vector fvec(1,4);
//	dvar_vector pvec(1,3);
//	
//	fvec.initialize();
//	pvec.initialize();
//	//Priors (pvec)
//	if(!last_phase())
//	{
//		pvec(1) = dnorm(ddot_r_devs,0,0.4);
//		pvec(2) = dnorm(bar_r_devs,0,0.4);
//		pvec(3) = dnorm(log_bar_f,log(0.1108032),0.05);
//	}
//	else
//	{
//		pvec(1) = dnorm(first_difference(ddot_r_devs),0,0.4);
//		pvec(2) = dnorm(bar_r_devs,0,2.5);
//		pvec(3) = dnorm(log_bar_f,log(0.1108032),2.5);
//	}
//	
//	//Likelihoods (fvec);
//	//fvec(1) = 10.*norm2(sqrt(C)-sqrt(Chat));
//	//fvec(2) = 10.*norm2(M-Mhat);
//	//fvec(3) = 10.*norm2(R-Rhat);
//	//cout<<fvec<<endl;	
//	/*dvariable mu =3.0;
//		dvariable size = 8.0;
//		cout<<dnbinom(0,mu,size)<<endl;
//		exit(1);*/
//	double tau;
//	double itau;
//	for(i=1;i<=nrow;i++)
//	{
//		tau      = sum(C(i));
//		fvec(1) += multifan(C(i),Chat(i),tau);
//		
//		itau      = sum(M(i));
//		//cout<<endl<<Mhat(i)(j_min(i),nx)<<endl<<endl;
//		fvec(2) += dmultifan(M(i)(j_min(i),nx),Mhat(i)(j_min(i),nx),itau);
//		
//		if(i>1){
//		itau      = sum(R(i));
//		fvec(3) += dmultifan(R(i)(j_min(i),nx),Rhat(i)(j_min(i),nx),itau);
//		}
//		/*Trying to figure out why I get nan's in fvec*/
//	}
//	dvar_vector delta = log(ct)-log(hat_ct);
//	if(!last_phase())
//	{
//		fvec(4)  = dnorm(delta,0,flag(2));
//	}
//	else
//	{
//		fvec(4)  = dnorm(delta,0,flag(3));
//	}
//	
//	/*
//	PRIORS for estimated model parameters from the control file.
//	*/
//	dvar_vector priors(1,npar);
//	dvariable ptmp; priors.initialize();
//	for(i=1;i<=npar;i++)
//	{
//		if(active(theta(i)))
//		{
//			switch(theta_prior(i))
//			{
//			case 1:		//normal
//				ptmp=dnorm(theta(i),theta_control(i,6),theta_control(i,7));
//				break;
//				
//			case 2:		//lognormal CHANGED RF found an error in dlnorm prior. rev 116
//				ptmp=dlnorm(theta(i),theta_control(i,6),theta_control(i,7));
//				break;
//				
//			case 3:		//beta distribution (0-1 scale)
//				double lb,ub;
//				lb=theta_lbnd(i);
//				ub=theta_ubnd(i);
//				ptmp=dbeta((theta(i)-lb)/(ub-lb),theta_control(i,6),theta_control(i,7));
//				break;
//				
//			case 4:		//gamma distribution
//				ptmp=dgamma(theta(i),theta_control(i,6),theta_control(i,7));
//				break;
//				
//			default:	//uniform density
//				ptmp=-log(1./(theta_control(i,3)-theta_control(i,2)));
//				break;
//			}
//			priors(i)=ptmp;	
//		}
//	}
//	
//	
//	if(fpen > 0 ) cout<<"Fpen = "<<fpen<<endl;
//	f = sum(fvec) + sum(pvec) + sum(priors);//	 + 1.e6*fpen;
//  }
//
//
//
FUNCTION dvar_vector dgamma(const dvector& x, const dvariable& a, const dvariable& b)
  {
	//returns the gamma density with a & b as parameters
	RETURN_ARRAYS_INCREMENT();
	dvariable t1 = 1./(pow(b,a)*exp(gammln(a)));
	dvar_vector t2 = (a-1.)*log(x)-x/b;
	RETURN_ARRAYS_DECREMENT();
	return(t1*exp(t2));
  }
//
FUNCTION dvariable dgamma(const prevariable& x, const double& a, const double& b)
  {
	//returns the gamma density with a & b as parameters
	RETURN_ARRAYS_INCREMENT();
	dvariable t1 = 1./(pow(b,a)*exp(gammln(a)));
	dvariable t2 = (a-1.)*log(x)-x/b;
	RETURN_ARRAYS_DECREMENT();
	return(t1*exp(t2));
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
//FUNCTION dvar_matrix dLTM(dvector& x, const dvariable &linf, const dvariable &k, const dvariable &cv)
//  {
//	/*
//		This function returns a length transition matrix (P_ij) with
//		the probability of an individual of length (x_i) growing into
//		the length interval (x_j). Length transition is based on the 
//		gamma distribution dgamma(x,a,b_l), where a=1/cv^2 and b_l
//		is given by b_l = (linf-x_j)*(1-exp(-k*dt))
//		
//		Args:
//			(x)  	- discrete length interval boundaries.
//			(linf)	- asymptotic length
//			(k)		- Brody growth coefficient over time interval dt.
//			(cv) 	- coefficient of variation in the size increment.
//			
//			Use the following for the mean growth increment to ensure
//			the function (dl) remains positive. Beyone Linf, growth 
//			increment decays exponentially.
//			dl	<- log(exp((linf-xp)*(1-exp(-k)))+1)  to ensure dl is positive
//			
//	*/
//	
//		RETURN_ARRAYS_INCREMENT();
//		int i,j,n;
//		n = size_count(x);
//		dvector xm(1,n-1);  //mid points of the length intervals
//		xm = x(1,n-1)+0.5*first_difference(x);
//		dvariable pen;
//		dvariable alpha = 1./square(cv);
//		dvar_vector beta(1,n-1);
//		dvar_vector dl(1,n-1);
//		dvar_matrix P(1,n-1,1,n-1);
//		cout<<"Ok"<<endl;
//		for(j=1;j<=n-1;j++)
//		{
//			//dvariable t1 = posfun(linf-xm(j),1./n,pen);
//			//beta(j) = t1 * (1.-mfexp(-k));
//			beta(j) = log(mfexp( (linf-xm(j))*(1.-mfexp(-k)) )+1.0);
//		}
//		dl = alpha*beta;
//		
//		for(i=1;i<=n-1;i++)
//		{
//			dvar_vector t2(i,n);
//			for(j=i;j<=n;j++)
//			{
//				double dx = x(j)-x(i);
//				t2(j) = cumd_gamma(dx,dl(i));
//			}
//			P(i)(i,n-1) = first_difference(t2);
//			P(i) /= sum(P(i));
//		}
//		
//		//cout<<P<<endl;
//		
//		RETURN_ARRAYS_DECREMENT();
//		return(P);
//  }
//
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
		dl(j) = log(mfexp( (linf-x(j))*(1.-mfexp(-k*dt)) )+1.0);
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
	
//REPORT_SECTION
//	int i,j,im;
//	REPORT(f          );
//	REPORT(log_ddot_r );
//	REPORT(log_bar_r  );
//	REPORT(log_bar_f  );
//	REPORT(m_infty    );
//	REPORT(l_infty    );
//	REPORT(vbk        );
//	REPORT(beta       );
//	REPORT(mu_r       );
//	REPORT(cv_r       );
//	
//	ivector yr(syr,nyr);
//	yr.fill_seqadd(syr,1);
//	REPORT(yr);
//	REPORT(xmid);
//	REPORT(mx);
//	REPORT(sx);
//	
//	dvector Nt(syr,nyr);
//	dvector Rt(syr,nyr);
//	Nt.initialize();
//	Rt.initialize();
//	
//	j = 1;
//	for(i=syr;i<=nyr;i++)
//	{
//		if(i==syr) Rt(i) = value(exp(log_ddot_r+ddot_r_devs(nx)));
//		else       Rt(i) = value(exp(log_rt(i)));
//		for(im=1;im<=int(1/dt);im++)
//		{
//			Nt(i) += sum(value(N(j++)));
//		}
//	}
//	REPORT(Nt);
//	REPORT(Rt);
//	
//	REPORT(log_rt);
//	REPORT(fi);
//	
//	REPORT(N);
//	REPORT(T);
//	REPORT(i_C);
//	REPORT(i_M);
//	REPORT(i_R);
//	REPORT(Chat);
//	REPORT(Mhat);
//	REPORT(Rhat);
//	
//	if(SimFlag)
//	{
//		REPORT(true_Nt);
//		REPORT(true_Rt);
//		REPORT(true_Tt);
//		REPORT(true_fi);
//	}
//
//
//TOP_OF_MAIN_SECTION
//	time(&start);
//	arrmblsize = 50000000;
//	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
//	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
//	gradient_structure::set_MAX_NVAR_OFFSET(5000);
//	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
// 
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
		RETURN_ARRAYS_INCREMENT();
		dvariable like;
		dvariable tau;
		int lb=o.indexmin();
		int nb=o.indexmax();
		int I = (nb-lb)+1;

		dvar_vector epsilon(lb,nb);
		dvar_vector O=o/sum(o);
		dvar_vector P=p/sum(p);

		tau=1./s;
		epsilon=elem_prod(1.-P,P);


		like=0.5*sum(log(2.*M_PI*(epsilon+0.1/I)))+I*log(sqrt(tau));
		like+= -1.*sum(log(mfexp(-1.*elem_div(square(O-P),2.*tau*(epsilon+0.1/I)))+0.01));
		
		RETURN_ARRAYS_DECREMENT();
		return like;
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
	
	
//FINAL_SECTION
//	time(&finish);
//	elapsed_time=difftime(finish,start);
//	hour=long(elapsed_time)/3600;
//	minute=long(elapsed_time)%3600/60;
//	second=(long(elapsed_time)%3600)%60;
//	cout<<endl<<endl<<"*******************************************"<<endl;
//	cout<<"--Start time: "<<ctime(&start)<<endl;
//	cout<<"--Finish time: "<<ctime(&finish)<<endl;
//	cout<<"--Runtime: ";
//	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
//	cout<<"*******************************************"<<endl;

