//*************************************
//**Model Description 				 **
//**Author  						 **
//**Date     						 **
//*************************************

DATA_SECTION
//This is the section where data is read in 
	init_number cntrl; 
	init_int yrs;
	init_int ages;
	init_number math;
	init_number mats;
	init_number vblinf;
	init_number vbk;
	init_number vbto;
	init_number lwa;
	init_number lwb;
	init_number iselh;
	init_number isels;
	init_number iRo;
	init_number ih;
	init_number pM;
	init_number psdM; 
	init_vector ocpue(1,yrs);
	init_vector olandings(1,yrs);
	init_matrix ocat(1,yrs,1,ages);
	init_int eof;
	matrix opat(1,yrs,1,ages);
	vector len(1,ages);
	vector wgt(1,ages);
	vector mat(1,ages);
	vector age(1,ages);
	vector fec(1,ages);
//Thes two lines set up a counter for MCMC Iterations
	int iter;
	!!iter=0;
//This local calulation heck that the data file was read in full
	LOC_CALCS
		if(eof != 999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}
		for(int i=1; i<=yrs; i++)
		{
			opat(i)=ocat(i)/sum(ocat(i));
		}
		opat=opat+1e-6;  //Moved from the procedure section so yu were changing the data.
		olandings = olandings/1e6;
		age.fill_seqadd(1.0,1.0);
		len=vblinf*(1.0-mfexp(-vbk*(age-vbto)));
		wgt=lwa*pow(len,lwb);
		mat=1.0/(1.0+exp(-mats*(age-math)));
		fec=0.5*elem_prod(wgt,mat);
	END_CALCS
	
PARAMETER_SECTION
//This is where all leading and derived parameters are declared	
	!!int init_dev_phz=2;
	!!if(cntrl) init_dev_phz=-1;
	init_bounded_number log_Ro(-40,40,4);
	init_bounded_number h(0.2,1.0,4);
	init_bounded_number log_avg_rec(-40,40,1);
	init_bounded_vector log_rec_dev(1,yrs,-15.,15.,2);
	init_number log_M;
	init_number log_selh(4);
	init_number log_sels(4);
	
	init_bounded_vector init_log_rec_dev(2,ages,-15.,15.,init_dev_phz);
	
	init_bounded_number log_Fbar(-30.,3.0,1);
	init_bounded_vector log_F_dev(1,yrs,-30.,3.0,2);
	init_bounded_number rho(0.001,0.999,3);//proportion of the observation error
	//Total precision is 1/var, so you need to increase the upper bound for this parameter.
	//Also estimate it in the last phase.
	init_bounded_number varphi(0.001,99,3);//total precision in the CPUE & Rec anomalies.
	//init_number log_rsig;
	//init_number log_zsig;
	//init_number log_csig;
	LOCAL_CALCS
		log_Ro=log(iRo);
		h=ih;
		log_avg_rec=.7*log_Ro;
		log_Fbar=log(0.3);
		log_M=log(pM);
		log_selh=log(iselh);
		log_selh=log(isels);
		//log_rsig=log(.4);
		//log_zsig=log(.4);
		//log_csig=log(.4);
	END_CALCS
	
	
	number Ro;
	vector Ft(1,yrs);
	number M;
	number selh;
	number sels;
	number rsig;
	number zsig;
	number csig;

	vector sel(1,ages);
	vector Mage(1,ages);
	vector lxo(1,ages);
	
	number ralpha;
	number rbeta;
	number phie;
	
	matrix Nt(1,yrs,1,ages);
	vector log_Rt(1,yrs);
	vector Rt_bar(2,yrs);
	vector rdelta(2,yrs);
	vector Bt(1,yrs);
	vector Zt(1,yrs);
	vector plandings(1,yrs);
	vector SSBt(1,yrs);
	matrix Fat(1,yrs,1,ages);
	matrix Zat(1,yrs,1,ages);
	matrix pcat(1,yrs,1,ages);
	matrix ppat(1,yrs,1,ages);

	vector cdev(1,yrs);
	vector zdev(1,yrs);
	matrix pdev(1,yrs,1,ages);
	number tau;
	number sig;
	number fpen;
	
	objective_function_value f;//This is necessary
	
	//sdreport_number ####;//This is needed if you want stdev estimated for a parameter or vector 
	//sdreport_vector ####;
	//likeprof_number ####;//This is needed if you want to get the likelihood profile for a parameter

PROCEDURE_SECTION
	//*** Main function calls  ***//
		initial_calculations();
		recursive_calculations();
		recruitment_calculations();
		observation_model();
		calc_objective_function();
		if(mceval_phase()) mcmc_stuff();		
	//****************************//
	
FUNCTION initial_calculations
	{
		fpen=0;
		Ro=mfexp(log_Ro);
		Ft=mfexp(log_Fbar+log_F_dev);
		M=mfexp(log_M);
		selh=mfexp(log_selh);
		sels=mfexp(log_sels);
		//rsig=mfexp(log_rsig);
		//csig=mfexp(log_csig);
		//zsig=mfexp(log_zsig);
		Mage=ages*M/sum(3.69*pow(wgt*1000.,-0.305))*3.69*pow(wgt*1000,-0.305);
		lxo(1)=1;
		for(int i=1; i<=ages; i++)
		{
			sel(i)=1.0/(1.0+exp(-sels*(i-selh)));
			if(i>1) lxo(i)=lxo(i-1)*mfexp(-Mage(i-1));
		}
		lxo(ages)/=(1-mfexp(-Mage(ages)));
		phie=sum(elem_prod(lxo,fec));
		ralpha=4.0*Ro*h/(5.0*h-1);
		rbeta=(phie*Ro*(1-h))/(5.0*h-1);
		Fat=outer_prod(Ft,sel);
		for(int i=1; i<=yrs; i++) Zat(i)=Fat(i)+Mage;
		// initialize recruitment and numbers
		if(cntrl){
			log_Rt(1)=log_Ro;
			Nt(1)=Ro*lxo;
		}else{
			cout<<"here"<<endl;
			ad_exit(1);
			Nt(1)(2,ages)=elem_prod(log_avg_rec+init_log_rec_dev,lxo(2,ages));
			Nt(ages)/=(1.-mfexp(-Mage(ages)));
		}
		tau=(1.-rho)/varphi;
		log_Rt=log_avg_rec+log_rec_dev - 0.5* tau*tau;
		Nt.colfill(1,mfexp(log_Rt));
		Bt(1)=sum(elem_prod(Nt(1),wgt));
		SSBt(1)=sum(elem_prod(Nt(1),fec));
		pcat(1)=elem_prod(elem_prod(elem_div(Fat(1),Zat(1)),elem_prod(Nt(1),(1-mfexp(-Zat(1))))),wgt);
		plandings(1)=sum(pcat(1));
		ppat(1)=pcat(1)/sum(pcat(1));
		//cout<<ralpha*SSBt(1)/(1+rbeta*SSBt(1))<<endl;
		//ad_exit(1);
	}	
FUNCTION recursive_calculations
	{
		for(int i=2; i<=yrs; i++)
		{
			Nt(i)(2,ages)=++elem_prod(Nt(i-1)(1,ages-1),mfexp(-Zat(i-1)(1,ages-1)));
			Nt(i)(ages)+=Nt(i-1,ages)*mfexp(-Zat(i-1,ages));
			Rt_bar(i)=ralpha*SSBt(i-1)/(rbeta+SSBt(i-1));
			Bt(i)=sum(elem_prod(Nt(i),wgt));
			SSBt(i)=sum(elem_prod(Nt(i),fec));
			pcat(i)=elem_prod(elem_div(Fat(i),Zat(i)),elem_prod(Nt(i),(1-mfexp(-Zat(i)))));
			ppat(i)=pcat(i)/sum(pcat(i));
			plandings(i)=sum(elem_prod(pcat(i),wgt));
		}
		Zt=log(ocpue)-log(Bt);
		ppat=ppat+1e-6;
		//opat=opat+1e-6;  //Note that every iteration you keep changing the data, by adding 1e-6.
	}
FUNCTION recruitment_calculations
	{
		
		rdelta=log(Rt_bar)-log_Rt(2,yrs);//+0.5*tau*tau;
	}
		
FUNCTION observation_model
	{
		sig=rho/varphi;
		cdev=log(olandings)-log(plandings);
		pdev=log(opat)-log(ppat);
		zdev=Zt-mean(Zt);
	}

FUNCTION calc_objective_function
	{
		dvar_vector likevec(1,8);
		likevec.initialize();
	    dvar_matrix tau2(1,yrs,1,ages);
	    dvariable s_tau2;
		//likelihoods
		//cpue residuals
		likevec(1) = dnorm(zdev,sig);
		// process residuals
		likevec(2) = dnorm(rdelta,tau);
		//catch residuals
		//likevec(3) = 0.5*(yrs-1)*log(sum(pow(cdev,2)));
		//likevec(3) = 0.5*(yrs-1)*log(norm2(cdev));
		likevec(3) = dnorm(cdev,0.01);
		//proprtions at age from Maunder 2011
		tau2=(1.0/(yrs*ages)*sum(elem_prod(ppat,pow(pdev,2))))/opat; 
		s_tau2 = (1.0)/(yrs*ages)*sum(elem_prod(ppat,square(log(opat)-log(ppat)))); 
		//cout<<1.0/(yrs*ages)*sum(elem_prod(opat,pow(pdev,2)))/opat<<endl;
		for(int i=1; i<=yrs; i++)
		{
			//if(min(tau2(i))<=0) exit(1);
			//likevec(4)+=dnorm(pdev(i),sqrt(tau2(i)));
			likevec(4)+=dnorm(pdev(i),sqrt(s_tau2/ppat(i)));
			//cout<<i<<endl;
			//ad_exit(1);
		}
		//ad_exit(1);
		
		if(last_phase())
		{
			likevec(5) = dnorm(log_F_dev,2.0);
			likevec(5)+=dnorm(log_Fbar,log(0.3),2.0);
			likevec(6) = dnorm(log_rec_dev,2.0);
			likevec(7) = dnorm(init_log_rec_dev,2.0);
		}
		else
		{
			likevec(5) =100.*norm2(log_F_dev);
			likevec(5)+=dnorm(log_Fbar,log(0.3),0.01);
			likevec(6)=100.*norm2(log_rec_dev);
			likevec(7)=100.*norm2(init_log_rec_dev);
		}
		likevec(8)=log(M)+log(psdM)+square(log(M)-log(pM))/(2*psdM*psdM);
		
		cout<<likevec<<endl;
		f=sum(likevec);
		f+= dbeta(rho,1.01,1.01);
		f+= dgamma(varphi,1.01,1.01);
		f+= dbeta((h-0.2)/0.8,1.01,1.01);
	}

REPORT_SECTION
	REPORT(Ro);
	REPORT(h);
	REPORT(Ft);
	REPORT(M);
	REPORT(selh);
	REPORT(sels);
	REPORT(sel)
	REPORT(rho);
	REPORT(varphi);
	REPORT(tau);
	REPORT(sig);
	REPORT(cdev);
	REPORT(zdev);
	REPORT(rdelta);
	
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
	#define REPORT(object) \
	report << #object "\n" \
	<< object << endl;

	#include <admodel.h>
	#include <time.h>
	#include <stats.cxx>

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

FUNCTION mcmc_stuff

		adstring fl1=adstring("LMCMC.rep");

	 	if(iter==0){
			//ofstream ofs1(fl1);
			//ofs1<<"HEADER STUFF\t"<<endl;
		}
		iter++;
		ofstream ofs1(fl1,ios::app);	
		//ofs1<<"dump mcmc results"<<endl;
