//******************************************************
//	Programmer: Steve Martell
//	Project Name:
//	Date:
//	Version:
//	Comments:  For the sturgeon case, need to integrate the SRA with ASMR
//	TO DO LIST FOR HBC PROJECT WORKSHOP:
//		1) Add aging error to simulation model (**watch array bounds in agetag**)
//		2) Extend simulation for 10 year projections with input Rt values (DONE)
//		3) Write code to produce simulated data file. (DONE)
//		4) Need to compare with ASMR3t (lews code) and repeat results
//		5) 
//		6) 
//		**The plus group code in here has not been tested yet. (Ignored for this project)
//******************************************************/
DATA_SECTION
	//SM Code for adding a command line option to the program.
	int sim;					//flag for running the simulator
	int rseed;					//random number seed for generating fake data
	int retyr;					//number of retrospective years
	number ptcap;				//capture probability for simulating future data
	int plus;					//flag for using +group calculations (termage)
	int scenario;				//which data scenario to use
	LOC_CALCS
		int on,opt;
		sim=0;
		plus=0;
		retyr=0;
		rseed = 123;		//default random number seed.
		ptcap=0;
		if((on = option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			sim=1;
			rseed=atoi(ad_comm::argv[on+1]);
			scenario=atoi(ad_comm::argv[on+2]);
			ptcap=atof(ad_comm::argv[on+3]);
			//retyr=-6;
			ofstream ofs("seed.txt");
			ofs<<rseed<<endl;
		}
		if((on = option_match(ad_comm::argc,ad_comm::argv,"-retro",opt))>-1)
		{
			retyr=atoi(ad_comm::argv[on+1]);
			sim=0;		//turn off simulation model if doing retrospective analysis
		}
		if((on = option_match(ad_comm::argc,ad_comm::argv,"-scenario",opt))>-1)
		{
			scenario=atoi(ad_comm::argv[on+1]);		//scenario # for future data
			ptcap=atof(ad_comm::argv[on+2]);
		}
	END_CALCS
	
	init_adstring data_file;			//file name of data
	init_adstring scenario_file;		//file for alternative data scenario
	init_int bycohort;					//switch to allow by tagcohort consideration of recaptures
	init_int model;						//switch for model type ASMR1, ASMR2, ASMR3
	init_int likelihood;				//switch for poisson or negative binomial likelihood.
	!! cout<<"__________________________________\n";
	!! cout<<"\t-Model type = "<<model<<endl;
	!! cout<<"\t-By cohort = "<<bycohort<<endl;
	!! cout<<"\t-Likelihood: ";
	!! if(likelihood==1) cout<<"Poisson"<<endl; else cout<<"Negative binomial"<<endl;
	!! cout<<"__________________________________\n";
		
	//_______________________________________________
	!! ad_comm::change_datafile_name(data_file);
	//_______________________________________________
	//_______________________________________________
	!! cout<<"Data file name: "<<data_file<<endl;
	init_int syr;						//start year of data
	init_int nnyr;						//end year of data
	int nyr;							//end of simulation year (for retrospective analysis)
	!! nyr = nnyr-retyr;
	!! if(retyr) cout<<"Retrospective year = "<<nyr<<endl;
	init_int sage;						//youngest age
	init_int nage;						//oldest age
	vector age(sage,nage);
	!! age.fill_seqadd(sage,1);
	!! cout<<"Age vector: "<<age(sage,sage+3)<<" ..."<<age(nage)<<endl;
	ivector byr(syr,nyr-2); //brood years
	!! byr.fill_seqadd(syr-2,1);

	init_number vonbk;					//von b k parameter
	init_number linf;					//Linfinity parameter
	init_vector lenage(sage,nage);			//length at age
	init_int termage;					//last age to estimate abundance in terminal year (ASMR 2 and 3)
	!! cout<<"Terminal age = "<<termage<<endl;
	
	init_int nva;						//number of vulnerability blocks
	init_ivector epics(1,nva);			//ending year of each vulnerability epoch

	init_int m_phz;						//solution phaze for m-deviations
	!! cout<<"Phase for estimating m deviations = "<<m_phz<<endl;
	
	!! cout<<scenario_file<<endl;
	!! if(!(scenario_file == "NULL")) ad_comm::change_datafile_name(scenario_file);
	
	init_matrix i_mta(syr,nnyr,sage,nage);	//Data, number of untagged fish captured by age and year
	init_matrix i_rta(syr,nnyr,sage,nage);	//Data, number of recaptured fish by age and year
	//Data, number of tagged fish recaptured by age and year and tag cohort
	init_3darray i_rcta(syr,nnyr,syr,nnyr,sage,nage);	
	!! cout<<i_rcta(nnyr-1)(nnyr)<<endl;
	
	matrix mta(syr,nyr,sage,nage);
	matrix rta(syr,nyr,sage,nage);
	3darray rcta(syr,nyr,syr,nyr,sage,nage);
	LOC_CALCS
		int i,n_yr;
		mta.initialize();
		cout<<syr<<" "<<nnyr<<endl;
		//some accouting here to deal with retrospective or projection.
		if(retyr<=0)n_yr = nnyr; else n_yr=nyr;
		for(i=syr;i<=n_yr;i++)
		{
			mta(i)=i_mta(i);
			rta(i)=i_rta(i);
			for(int j=syr;j<=n_yr;j++)
				rcta(i)(j)=i_rcta(i)(j);
		}
		epics(nva)=nyr;
		cout<<"Ok here?\n"<<endl;
	END_CALCS
	
		
	matrix epsilon(syr,nyr,sage,nage);
	matrix delta(syr,nyr,sage,nage);
	3darray Cdelta(syr,nyr,syr,nyr,sage,nage);
	init_int eof;
	!! if(eof!=999){cout<<"Error reading data file\n"; exit(1);}
	!! cout<<"\n___---*** END OF DATA FILE ***---___\n\n \t\t"<<eof<<endl;
	
	int tau_phz;
	int UT_phz;
	int pterm_phz;
	int ahat_phz;
	int nf;
	
	//_______________________________________________
	!! ad_comm::change_datafile_name("Simdata.dat");
	//_______________________________________________
	//_______________________________________________
	//input data for simulation model.
	init_vector sim_pt(syr,nnyr);
	init_vector sim_rt(syr-10,nnyr);
	
	//_______________________________________________
	!! ad_comm::change_datafile_name("Selectivity.dat");
	//_______________________________________________
	//_______________________________________________
	//input selectivities for simulation model.
	init_matrix sel_coffs(sage,nage,1,7);  //selectivity coefficients from 7 data scenarios
	int nsim;
	!! nsim=10;	//number of simulation years into the future
	matrix o_mta(syr,nyr+nsim,sage,nage);
	matrix o_rta(syr,nyr+nsim,sage,nage);
	3darray o_rcta(syr,nyr+nsim,syr,nyr+nsim,sage,nage);
	
	
	//Classification matrix for aging errors
	matrix Q(sage,nage,sage,nage);
	LOC_CALCS
		dvector ia(sage,nage);
		dvector std(sage,nage);
		//std = 0.1*age;
		std = 0.2+(15-0.2)*plogis(age,15.,3.);
		ia.fill_seqadd(0,1);
		//cout<<ia<<endl;
		double z1,z2;
		for(int i = sage;i<=nage;i++){
			for(int j=sage; j<=nage;j++)
			{
				z1=(ia(i)+0.5 - ia(j))/std(i);
				z2=(ia(i)-0.5 - ia(j))/std(i);
				Q(i,j) = cumd_norm(z1)-cumd_norm(z2);
			}
			Q(i)=Q(i)/sum(Q(i));
		}
		//Q=inv(Q);
		Q=identity_matrix(sage,nage);
		//cout<<setprecision(2)<<Q<<endl;
		//exit(1);
		
	END_CALCS 
	
PARAMETER_SECTION
	init_bounded_number m(0,1.0,-2);
	!! if(model==1) pterm_phz=1; else pterm_phz=-1;
	init_bounded_number pterm(0,1,pterm_phz);
	!! if(likelihood==2) tau_phz=2; else tau_phz=-1;
	init_bounded_number tau(0,500,tau_phz);
	!! if(model==3) ahat_phz=-1; else ahat_phz=1;
	init_bounded_vector ahat(1,nva,0,nage,ahat_phz);
	init_bounded_vector tau_ahat(1,nva,0,50,ahat_phz);
	init_bounded_vector g2(1,nva,0,1,-2);
	
	!! if(model==2||model==3) UT_phz=1; else UT_phz=-1;
	//init_bounded_vector log_UT(sage+1,termage,-10,10,UT_phz);
	init_vector log_UT(sage+1,termage,UT_phz);
	init_bounded_dev_vector m_dev(syr,nyr,-3.0,3.0,m_phz);
	
	!! m= 0.13;
	!! tau=25.0;
	!! ahat = 3.0;
	!! tau_ahat=2.0;
	!! log_UT=5.0;
	!! g2=0.01;
	
	objective_function_value f;
	sdreport_number sd_nt4_syr;
   	sdreport_number sd_nt4_nyr;
	sdreport_number sd_rt_syr;
   	sdreport_number sd_rt_ppyr;   
	
	vector sa(sage,nage);
	vector pt(syr,nyr);					//capture probabilities
	
	
	matrix va(syr,nyr,sage,nage);
	matrix U(syr,nyr,sage,nage);	//unmarked fish at large
	matrix M(syr,nyr,sage,nage);	//marked fish at large
	matrix Mta(syr,nyr,sage,nage);	//new marks deployed
	matrix Rta(syr,nyr,sage,nage);	//recaptures
	matrix pta(syr,nyr,sage,nage);	//age-specific capture probabilities
	
	3darray Mc(syr,nyr,syr,nyr,sage,nage);		//Unmarked cohorts vulnerable to capture
	3darray Rcta(syr,nyr,syr,nyr,sage,nage);	//Predicted recaptures by cohort.

	//random_effects_vector u(syr,nyr,2);

PRELIMINARY_CALCS_SECTION
   nf=0;
	//call the data simulator from here.
	//Note, you cannot use tabs in this section
    if(sim){
        cout<<"______________________________________________________\n"<<endl;
        cout<<"    **Implementing Simulation--Estimation trial** "<<endl;
        cout<<"______________________________________________________"<<endl;
        cout<<"\tRandom Seed No.:\t"<< rseed<<endl;
        cout<<"______________________________________________________\n"<<endl;
        
        //Call the simulator  (This appears to be working now.)
        //Still need to add aging error to simulation model
        //simulation_model();
        //exit(1);
    }


PROCEDURE_SECTION
	//   *************************************  \\
	//				MAIN ROUTINES				\\
	//  **************************************  \\
	
	get_survival();
	get_vulnerabilities();
	unmarked_and_marked_fish();	
	get_capture_probabilities();
	calc_negloglike();
	dvar_vector nt4=rowsum(trans(trans(U).sub(4,nage)))+rowsum(trans(trans(M).sub(4,nage)));
	dvar_vector rt=get_drt();
	sd_nt4_syr = log(nt4(syr));	
	sd_nt4_nyr = log(nt4(nyr));
	sd_rt_syr = log(rt(syr));	
	sd_rt_ppyr = log(rt(max(byr)));
	
	if(mceval_phase()) write_mcmc_report();
	//  **************************************  \\
	//  **************************************  \\
	
	
	
	
FUNCTION get_survival
	//to do: add mdevs
	//survival rate  ma=m*linf/lenage
	dvector t1 = linf/lenage;
	sa = mfexp(-m*t1);
	//cout<<"sa\n"<<sa<<endl;
	//exit(1);
	// May 2, 2010. SM compared this routine with Coggins ASMR3t.tpl
	// it gives the exact same results.
	
FUNCTION get_vulnerabilities	
	int i,j;
	//vulnerabilities
	j=1;
	for(i=syr;i<=nyr;i++)
	{
		if(i>epics(j)) j++;
		va(i)=plogis(age,ahat(j),tau_ahat(j));
		if(active(g2)) va(i)=eplogis(age,ahat(j),tau_ahat(j),g2(j));
	}
	
FUNCTION unmarked_and_marked_fish
	int i,j,k;
	double tiny = 1.e-30;
	
	//terminal abundances and VPA reconstruction
	U.initialize();
	U=tiny;
	U(nyr,sage)=999.;
	
	//Add catches to terminal age column
	dvector d1 = column(mta,nage);
	for(i=syr; i<=nyr;i++)U(i,nage)+=d1(i);

	if(model==1) U(nyr) = elem_div((mta(nyr)+tiny),(pterm*va(nyr)));
	
	if(model==2 || model==3){ //ASMR 2
		U(nyr)(sage+1,termage)+=exp(log_UT);
		for(j=termage+1; j<=nage; j++)
		{
			U(nyr,j) = U(nyr,j-1)*sa(j-1);
		}
	}
	for(i=nyr-1; i>=syr; i--)
	{
		for(j=nage-1; j>=sage; j--)
		{
			U(i,j) = (U(i+1,j+1)/sa(j)) + mta(i,j);	
		}
	}
	
	
	//marked fish
	if(!bycohort)
	{
		//marked fish
		M.initialize();
		for(i=syr+1;i<=nyr;i++)
		{
			M(i)(sage+1,nage) =++elem_prod(M(i-1)(sage,nage-1)+mta(i-1)(sage,nage-1),sa(sage,nage-1));
		}
	}
	
	if(bycohort)
	{
		Mc.initialize();
		for(k=syr;k<=nyr;k++){ //loop over cohorts
			Mc(k)=tiny;
			for(i=k+1;i<=nyr;i++){ //loop over cohort years
				for(j=sage+1;j<=nage;j++) //loop over ages
				{
					if(i==k+1){Mc(k,i,j)= mta(i-1,j-1)*sa(j-1);}
					else {Mc(k,i,j)=Mc(k,i-1,j-1)*sa(j-1);}
				}
			}
		}
		M.initialize();
		M=tiny;
		for(k=syr;k<=nyr;k++)
			for(i=syr;i<=nyr;i++)
				for(j=sage;j<=nage;j++)
					M(i,j) += Mc(k,i,j);
					
		//cout<<M<<endl;
		//exit(1);
	}


FUNCTION get_capture_probabilities
	int i,j;
	dvariable t1,t2;
	Rta.initialize();
	Mta.initialize();
	Rcta.initialize(); 
	dmatrix iQ = trans(trans(Q.sub(sage+1,nage)).sub(sage+1,nage));
	if(model != 3)
	{
		for(i=syr; i<=nyr;i++)
		{
			t1=sum(mta(i)+rta(i));
			t2=sum(elem_prod(va(i),U(i)+M(i)));
			pt(i)=t1/t2;
			Mta(i) = elem_prod(U(i),va(i)*pt(i));
			Rta(i) = elem_prod(M(i),va(i)*pt(i));
			
			if(bycohort)
			{
				for(j=syr;j<=nyr;j++)
				{
					Rcta(i)(j)=elem_prod(Mc(i)(j),va(i)*pt(i));
				}
			}
		}
	}

	if(model==3)  //ASMR3
	{
		dvar_vector tmp(sage,nage);
		for(i=syr; i<=nyr; i++)
		{
			tmp = elem_div(mta(i)+rta(i),U(i)+M(i));  //MLE of capture probability
			//cout<<mta(i)+rta(i)<<endl;
			Mta(i) = elem_prod(U(i),tmp);
			Rta(i) = elem_prod(M(i),tmp);
			
			t1=sum(mta(i)+rta(i));
			t2=sum(U(i)+M(i));
			pt(i)=t1/t2;
			pta(i)=tmp;
			if(bycohort)
			{
				for(j=syr;j<=nyr;j++)
				{
					tmp = elem_div(mta(j)+rta(j),U(j)+M(j));  //MLE of capture probability
					Rcta(i)(j)=elem_prod(Mc(i)(j),tmp+1.e-30);
				}
			}
		}
		//cout<<setprecision(2)<<endl<<Rcta(syr)<<endl;
		//exit(1);
	}
	//cout<<endl<<Mta<<endl<<endl;
FUNCTION calc_negloglike
	/*
		There are two options here, 1: poisson, and 2: negative binomial
		
	*/
	int i,j;
	dvar_vector lvec(1,4);
	lvec.initialize();
	epsilon.initialize();
	delta.initialize();
	Cdelta.initialize();
	
	
	//Likelihood of the data given a poisson distribution
	if(likelihood==1)
	{
		for(i=syr; i<=nyr;i++)
		{
			lvec(1)+=dpois_residual(mta(i)(sage,nage),Mta(i)(sage,nage),epsilon(i));
			
			if(!bycohort)
				lvec(2)+=dpois_residual(rta(i)(sage+1,nage),Rta(i)(sage+1,nage),delta(i));
				
			if(bycohort)
			{
				for(j=syr;j<=nyr;j++)
				{
					lvec(2)+=dpois_residual(rcta(i)(j)(sage+1,nage),Rcta(i)(j)(sage+1,nage),Cdelta(i)(j));
				}
			}
		}
	}
	
	//Likelihood of the data given a negative binomial distribution.
	if(likelihood==2)
	{
		for(i=syr; i<=nyr;i++)
		{
			lvec(1)+=dnbinom(mta(i)(sage,nage),Mta(i)(sage,nage),tau,epsilon(i));
			
			if(!bycohort)
				lvec(2)+=dnbinom(rta(i)(sage+1,nage),Rta(i)(sage+1,nage),tau,delta(i));
				
			if(bycohort)
			{
				for(j=syr;j<=nyr;j++)
				{
					lvec(2)+=dnbinom(rcta(i)(j)(sage+1,nage),Rcta(i)(j)(sage+1,nage),tau,Cdelta(i)(j));
				}
			}
		}
	} 
	
	//if(active(m)) lvec(3) = dlnorm(m,log(0.13-0.0005),0.1); 	//prior for adult m.
	//lvec(4) = 0.1*square(log_UT(sage)-log(mta(nyr,sage)/pt(nyr)));	//use termial capture prob for age2 in terminal year
	//cout<<"Likelihood\t"<<lvec<<endl;
	
	f=sum(lvec);
	nf++;
	//exit(1);
	
TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);

GLOBALS_SECTION
	#include <admodel.h>
	#include <time.h>
	#include <ASMR.cxx>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
	
	
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
	
	if(!sim) write_trial();
	if(sim)simulate_future_data(); 

FUNCTION write_trial
	//write seed, pt, nyr, rt(2009,2015), 
	int seed;
	ifstream ifs("seed.txt");
	ifs>>seed;
	ofstream ofs("Trials.txt",ios::app);
	adstring tt="\t";
	dvector rt=get_rt();
	ofs<<scenario<<tt<<seed<<tt<<ptcap<<tt<<rt(2009,2015)<<endl;


REPORT_SECTION
	cout<<"Report Section: Nyr is = "<<nyr<<endl;
	report<<"M\n"<<m<<endl;
	report<<"age\n"<<age<<endl;
	dvector yr(syr,nyr);
	yr.fill_seqadd(syr,1);
	report<<"yr\n"<<yr<<endl;

	report<<"va\n"<<va<<endl;
	report<<"sa\n"<<sa<<endl;
	report<<"pt\n"<<pt<<endl;
	report<<"pta\n"<<pta<<endl;
	report<<"epsilon\n"<<epsilon<<endl;
	report<<"delta\n"<<delta<<endl;
	report<<"mta\n"<<mta<<endl;
	report<<"rta\n"<<rta<<endl;
	report<<"nt2\n"<<rowsum(trans(trans(U).sub(2,nage)))+rowsum(trans(trans(M).sub(2,nage)))<<endl;
	report<<"nt3\n"<<rowsum(trans(trans(U).sub(3,nage)))+rowsum(trans(trans(M).sub(3,nage)))<<endl;
	report<<"nt4\n"<<rowsum(trans(trans(U).sub(4,nage)))+rowsum(trans(trans(M).sub(4,nage)))<<endl;
	REPORT(byr);
	report<<"rt\n"<<get_rt()<<endl<<endl;
	report<<"Q\n"<<Q<<endl;
	
	if(model==3) report<<"cdelta\n"<<Cdelta<<endl;

FUNCTION dvar_vector get_drt()
	//back-calculated-recruitment
	int pyr = min(byr);
	int ppyr = max(byr);
	dvar_vector rt(pyr,ppyr); rt.initialize();
	rt(syr,ppyr)=(column(U,sage)(syr,ppyr));
	dvar_vector lx(sage,nage); lx(sage)=1;
	for(int i=sage;i<sage+2;i++)
	{
		lx(i+1) = lx(i)*(sa(i));
		rt(syr-i+(sage-1)) = (U(syr,i+1))/lx(i+1);
	}
	return(rt); 
	
FUNCTION dvector get_rt()
	//back-calculated-recruitment
	int pyr = min(byr);
	int ppyr = max(byr);
	dvector rt(pyr,ppyr); rt.initialize();
	rt(syr,ppyr)=value(column(U,sage)(syr,ppyr));
	dvector lx(sage,nage); lx(sage)=1;
	for(int i=sage;i<sage+2;i++)
	{
		lx(i+1) = lx(i)*value(sa(i));
		rt(syr-i+(sage-1)) = value(U(syr,i+1))/lx(i+1);
	}
	return(rt);

FUNCTION write_mcmc_report
	if(nf==1)
	{
		ofstream ofs("asmr.mcmc");
		ofs<<"m\t tau\t pterm\t"<<endl;
		ofstream of2("nt4.mcmc");
		ofstream of3("rt.mcmc");
		ofstream of4("pt.mcmc");
	}
	ofstream ofs("asmr.mcmc",ios::app);
	ofs<<m<<"\t"<<tau<<"\t"
		 <<pterm<<"\t"<<endl;

	ofstream of2("nt4.mcmc",ios::app);
	of2<<rowsum(trans(trans(U).sub(4,nage))+trans(trans(M).sub(4,nage)))<<endl;
	
	ofstream of3("rt.mcmc",ios::app);
	of3<<get_rt()<<endl;
	
	ofstream of4("pt.mcmc",ios::app);
	of4<<pt<<endl;

FUNCTION int length2age(const double& LenTag)
	//SM   Need to change this polynomial to the one that is used to convert length to age from
	//the real data in the Rcode.   DONE
	//ifelse(round(-1.05952E-08*len^4+1.25697E-05*len^3-0.005058694*len^2
	//	+0.878738367*len^1-53.55077549)>alast,alast,
	
	int age=(-1.05952E-08*pow(LenTag,4)+1.25697E-05*pow(LenTag,3)-0.005058694*pow(LenTag,2)+0.878738367*LenTag-53.55077549);
	age<sage?age=sage:NULL;
	return(age);
	//return ( -0.0000000019605 * pow(LenTag,4) + 0.0000019314 * pow(LenTag,3) - 0.00039319 * pow(LenTag,2) + 0.024165 * LenTag + 1.8451);
	//return ceil( -0.0000000019605 * pow(LenTag,4) + 0.0000019314 * pow(LenTag,3) - 0.00039319 * pow(LenTag,2) + 0.024165 * LenTag + 1.8451);
	
	
	//inverse of vonb growth, very approximate
	/*int agetag = nage;
		if(linf>LenTag) agetag = -1./vonbk * log(1.-LenTag/linf);
		//cout<<agetag<<"\t"<<vonbk<<endl;
		return (agetag);*/
	
FUNCTION void simulate_future_data()
	/*
	This function is used to concatenate future data (post 2009) to
	the full data set (Scenario 1) using an IBM model where the MLE
	values from S1 are used to condition this simulation model. Capture
	probabilities will differ according to changes in the sampling
	regimes (Scenarios 2-7) and should result in different data sets
	for a given random number seed owing to minor differences in 
	capture probabilities and reduced sampling effort.
	
	Inputs into this routine include:
		- A vector of future recruits
	*/ 
	cout<<"****\n"<<"Simulating Future Data\n"<<"****\n"<<endl;
	
	int i,j,k;			//index for year, individual, age=k
	random_number_generator rng(rseed);
	int agetag, yeartag;
	bool tag;
	double cv_len=0.1;
	double sigma=3.0;
	double pt=ptcap;
	
	dvector rt(nyr+1,nyr+nsim); //future age-2 recruits  
	dvector s_nt4(syr,nyr+nsim);
	
	/*------------------------------------------------------------------
	The following code deals with taging and recapturing new recruits
	only.  
	pcap is based on the general capture probability and the age-specic
	selectivity cofficients that come from each of the 7 scenarios.
	I'm not completely satisfied with this approach because there is a 
	lot of variability in the pta matrix in terms of age-specific 
	selectivity coefficients.
	------------------------------------------------------------------*/
	
	for(i=syr;i<=nyr;i++)
	{
		o_mta(i)=mta(i);
		o_rta(i)=rta(i);
	}
	
	rt=2000; rt(2011)=1000;
	for(i=nyr+1;i<=nyr+nsim;i++)
	{
		for(j=1;j<=rt(i);j++) //loop over individuals
		{   
			tag=false;
			agetag=0;
			double l_dev=cv_len*randn(rng);
			for(k=sage;k<=nage;k++)
			{
				int iyr = i+k-sage;		//index for cohort
				//cout<<"year\t"<<iyr<<endl;
				if(iyr>nyr+nsim) break;		//break if past nyr+nsim
				double len=lenage(k)+l_dev*lenage(k);
				double survival=exp(-value(m)*linf/len);	//length based survival rate
				//double pcap = pt(iyr)*plogis(len,100,6);	//check selectivity from results 
				//Change pcap to be age-based and use lowess 1/4 estimates from average pta values
				double pcap=pt*sel_coffs(k,scenario);
				double xx=randu(rng);
				if(xx<=pcap)	//individual is captured or recaptured.
				{    
					if(!tag)	//first capture
					{
						tag=true;
						double lentag=len+sigma*randn(rng);	//measured length at tagging
						//agetag=k;	//add aging error here from ALK.
						//cout<<sage<<" "<<k<<endl;
						agetag=length2age(lentag);  //assignment of age based on Lew's polynomial
						yeartag=iyr;
						 
						o_mta(iyr,agetag)++;
						//mta(iyr,agetag)++;
					}
					else		//recapture
					{   
						o_rta(iyr,agetag)++;
						o_rcta(yeartag,iyr,agetag)++;
					}
				} 
				if(k>=4)s_nt4(iyr)++;
				if(randu(rng)>survival) break;		//fish dies
				if(agetag!=0) agetag++;				//fish lives another year
				
			}
		}
		//cout<<"Year\t"<<i<<endl;
	}
	//cout<<o_mta<<endl;
	//exit(1);
	
	/*The following code uses the MLE estimates of Untagged individuals U
	to generate new marks.
	*/
	int iage;
	ivector iU=value(U(nyr)); //number of unmarked animals in the terminal year
	iU(sage)=rt(nyr+1);			  //assume 2000 recruits in terminal year
	//cout<<"iU\n"<<iU<<"-----\n"<<endl;
	for(iage=sage+1;iage<=nage;iage++)
	{   
		for(j=1;j<=iU(iage);j++)
		{   
			tag=false;
			agetag=0;
			double l_dev=cv_len*randn(rng);
			for(k=iage;k<=nage;k++)
			{
             	//int iyr = nyr+k-sage+1;		//index for cohort
				int iyr = nyr+1+k-iage;
				//cout<<"year\t"<<iyr<<endl;
				
				if(iyr>nyr+nsim) break;		//break if past nyr+nsim
				double len=lenage(k)+l_dev*lenage(k);
				double survival=exp(-value(m)*linf/len);	//length based survival rate
				//Change pcap to be age-based and use lowess 1/4 estimates from average pta values
				double pcap=pt*sel_coffs(k,scenario);
				double xx=randu(rng);
				if(xx<=pcap)	//individual is captured or recaptured.
				{    
					if(!tag)	//first capture
					{
						tag=true;
						double lentag=len+sigma*randn(rng);	//measured length at tagging
						//agetag=k;	//add aging error here from ALK.
						agetag=length2age(lentag);  //assignment of age based on Lew's polynomial
						//cout<<iyr<<" "<<k<<" "<<agetag<<endl;
						yeartag=iyr; 
						o_mta(iyr,agetag)++;
						//mta(iyr,agetag)++;
					}
					else		//recapture
					{   
						o_rta(iyr,agetag)++;
						o_rcta(yeartag,iyr,agetag)++;
					}
				}
				if(k>=4)s_nt4(iyr)++;
				if(randu(rng)>survival) break;				//fish dies
				if(agetag!=0) agetag++;				//fish lives another year
				  
			}
			
		}
	}
	//cout<<o_mta<<endl;
	//exit(1);
	/*Next need to write the code to look for recaptures of previously marked cohorts
	from nyr*/
	ivector iM=value(M(nyr));
	//cout<<iM<<endl;
	for(iage=sage+1;iage<=nage;iage++)
	{
		for(j=1;j<=iM(iage);j++)
		{   
			double l_dev=cv_len*randn(rng);
			for(k=iage;k<=nage;k++)
			{
				int iyr = nyr+k-iage+1;		//index for cohort
				if(iyr>nyr+nsim) break;		//break if past nyr+nsim
				double len=lenage(k)+l_dev*lenage(k);
				double survival=exp(-value(m)*linf/len);	//length based survival rate
				
				//Change pcap to be age-based and use lowess 1/4 estimates from average pta values
				double pcap=pt*sel_coffs(k,scenario);
				double xx=randu(rng);
				if(xx<=pcap)
				{
					/*All fish are tagged here, so these are only recaptures*/ 
					//cout<<iyr<<" "<<k<<endl;
					o_rta(iyr,k)++;
					//not sure how to deal with cohort based groups, b/c I don't know the year tagged for iM
				}
				if(k>=4)s_nt4(iyr)++;
				if(randu(rng)>survival) break;				//fish dies 
			}
		}
	}
    //cout<<o_rta<<endl;

	ofstream ofs("SimData.dat");
	ofs<<syr<<endl<<nyr+nsim<<endl;
	ofs<<sage<<endl<<nage<<endl;
	ofs<<vonbk<<endl<<linf<<endl;
	ofs<<"#Length-at-age\n"<<lenage<<endl;
	ofs<<"#Term age\n"<<termage<<endl;
	ofs<<nva<<endl<<epics<<endl;
	ofs<<"#M_phz\n"<<m_phz<<endl;
	ofs<<"#mta\n"<<o_mta<<endl;
	ofs<<"#rta\n"<<o_rta<<endl;
	ofs<<"#rcta\n"<<o_rcta<<endl;
	ofs<<"#eof\n"<<999<<endl;
    ofs.close();

	//Simulation values used for comparing estimates
	ofstream ofss("SimValues.dat");
	ivector yr(syr,nyr+nsim);
	yr.fill_seqadd(syr,1);
	dvector rtt(syr,nyr+nsim);
	
	rtt(syr,nyr)=value(column(U,sage));
	rtt(nyr+1,nyr+nsim)=rt(nyr+1,nyr+nsim);
	rtt(nyr)=iU(sage);
	s_nt4(syr,nyr)=value(rowsum(trans(trans(U).sub(4,nage))+trans(trans(M).sub(4,nage))));
	
	ofss<<"yr\n"<<yr<<endl;
	ofss<<"rt\n"<<rtt<<endl;
	ofss<<"nt4\n"<<s_nt4<<endl;
	cout<<"*** ---FINISHED WRITING SIMULATED DATA--- ***"<<endl;
	
	
	

FUNCTION void simulation_model()
	random_number_generator rng(rseed);
	long ntag=sum(mta);
	long recaps=sum(rta);
	mta.initialize();
	rta.initialize();
	rcta.initialize();
	
	int i,j,k; 						//index for year, individual,age;
	int ssyr=syr-(nage-sage+1); 	//birth year of oldest cohort in syr;
	int agetag,yeartag;
	bool tag;						//boolean values for capture and recapture
	double cv_len=0.1;
	double sigma=3.0;
	dvector rt(ssyr,nyr);
	dvector pt(ssyr,nyr);			//capture probability
	dvector nt4(ssyr,nyr);			//Age-4+ abundance.
	dmatrix tmp_mta(syr,nyr,sage,nage); tmp_mta.initialize();
	
	if(rseed % 2)rt=3500+700; else rt=3500-700;
	
	//SM Notes:  get historical recruits and capture probabilities from Scenario 1
	//Keep historical data the same and project new program from 2009 forwards.
	//I.e., make this like an MSE evaluating changes in capture probabilities
	//associated with changes in sampling design.
	
	rt(syr-10,nnyr)=sim_rt;
	//rt(2009)=3500-700;
	pt=0.2; pt(ssyr,syr-1)=0;
	pt(syr,nnyr)=sim_pt(syr,nnyr);
	
	for(i=ssyr; i<=nyr; i++)
	{
		for(j=1; j<=rt(i); j++)		//loop over individuals
		{
			tag=false;
			agetag=0;
			double l_dev=cv_len*randn(rng);
			for(k=sage; k<=nage; k++)
			{
				int iyr = i+k-sage;		//index for cohort
				//cout<<"year\t"<<iyr<<endl;
				if(iyr>nyr) break;		//break if past nyr
				double len=lenage(k)+l_dev*lenage(k);
				double survival=exp(-value(m)*linf/len);	//length based survival rate
				double pcap = pt(iyr)*plogis(len,100,6);	//check selectivity from results
				double xx=randu(rng);
				if(xx<=pcap)	//individual is captured or recaptured.
				{
					if(!tag)	//first capture
					{
						tag=true;
						double lentag=len+sigma*randn(rng);	//measured length at tagging
						agetag=k;	//add aging error here from ALK.
						//agetag=length2age(lentag);
						yeartag=iyr; 
						
						mta(iyr,agetag)++;
					}
					else		//recapture
					{
						rta(iyr,agetag)++;
						rcta(yeartag,iyr,agetag)++;
					}
				}
				if(k>=4)nt4(iyr)++;					//age-4 or older
				if(randu(rng)>survival) break;		//fish dies
				if(agetag!=0) agetag++;				//fish lives another year
			}
		}
	}
	
	cout<<"******************************************************"<<endl;
	cout<<"**\t\t\tObserved\tSimulated"<<endl;
	cout<<"** -Marks released:\t"<<ntag<<"\t\t"<<sum(mta)<<endl;
	cout<<"** -Marks recaptured:\t"<<recaps<<"\t\t"<<sum(rta)<<endl;
	cout<<"******************************************************"<<endl;
	cout<<"******************************************************"<<endl;
	
	//-- Write true abundance of age 4+ fish and age-2 recruits
	ofstream ofs("SimValues.dat");
	ivector yr(ssyr,nyr);
	yr.fill_seqadd(ssyr,1);
	ofs<<"yr\n"<<yr<<endl;
	ofs<<"rt\n"<<rt<<endl;
	ofs<<"nt4\n"<<nt4<<endl;
	
	write_simulation_data();

FUNCTION void write_simulation_data()
	//-- This routine writes simulation data
	
	ofstream ofs("SimulationData.dat");
	ofs<<syr<<endl<<nyr<<endl;
	ofs<<sage<<endl<<nage<<endl;
	ofs<<vonbk<<endl<<linf<<endl;
	ofs<<"#Length-at-age\n"<<lenage<<endl;
	ofs<<"#Term age\n"<<termage<<endl;
	ofs<<nva<<endl<<epics<<endl;
	ofs<<"#M_phz\n"<<m_phz<<endl;
	ofs<<"#mta\n"<<mta<<endl;
	ofs<<"#rta\n"<<rta<<endl;
	ofs<<"#rcta\n"<<rcta<<endl;
	ofs<<"#eof\n"<<999<<endl;

	cout<<"*** ---FINISHED WRITING SIMULATED DATA--- ***"<<endl;



