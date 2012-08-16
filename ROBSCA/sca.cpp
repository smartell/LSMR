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
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <sca.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  cntrl.allocate("cntrl");
  yrs.allocate("yrs");
  ages.allocate("ages");
  math.allocate("math");
  mats.allocate("mats");
  vblinf.allocate("vblinf");
  vbk.allocate("vbk");
  vbto.allocate("vbto");
  lwa.allocate("lwa");
  lwb.allocate("lwb");
  iselh.allocate("iselh");
  isels.allocate("isels");
  iRo.allocate("iRo");
  ih.allocate("ih");
  pM.allocate("pM");
  psdM.allocate("psdM");
  ocpue.allocate(1,yrs,"ocpue");
  olandings.allocate(1,yrs,"olandings");
  ocat.allocate(1,yrs,1,ages,"ocat");
  eof.allocate("eof");
  opat.allocate(1,yrs,1,ages);
  len.allocate(1,ages);
  wgt.allocate(1,ages);
  mat.allocate(1,ages);
  age.allocate(1,ages);
  fec.allocate(1,ages);
iter=0;
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
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
int init_dev_phz=2;
if(cntrl) init_dev_phz=-1;
  log_Ro.allocate(-40,40,4,"log_Ro");
  h.allocate(0.2,1.0,4,"h");
  log_avg_rec.allocate(-40,40,1,"log_avg_rec");
  log_rec_dev.allocate(1,yrs,-15.,15.,2,"log_rec_dev");
  log_M.allocate("log_M");
  log_selh.allocate(4,"log_selh");
  log_sels.allocate(4,"log_sels");
  init_log_rec_dev.allocate(2,ages,-15.,15.,init_dev_phz,"init_log_rec_dev");
  log_Fbar.allocate(-30.,3.0,1,"log_Fbar");
  log_F_dev.allocate(1,yrs,-30.,3.0,2,"log_F_dev");
  rho.allocate(0.001,0.999,3,"rho");
  varphi.allocate(0.001,99,3,"varphi");
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
  Ro.allocate("Ro");
  #ifndef NO_AD_INITIALIZE
  Ro.initialize();
  #endif
  Ft.allocate(1,yrs,"Ft");
  #ifndef NO_AD_INITIALIZE
    Ft.initialize();
  #endif
  M.allocate("M");
  #ifndef NO_AD_INITIALIZE
  M.initialize();
  #endif
  selh.allocate("selh");
  #ifndef NO_AD_INITIALIZE
  selh.initialize();
  #endif
  sels.allocate("sels");
  #ifndef NO_AD_INITIALIZE
  sels.initialize();
  #endif
  rsig.allocate("rsig");
  #ifndef NO_AD_INITIALIZE
  rsig.initialize();
  #endif
  zsig.allocate("zsig");
  #ifndef NO_AD_INITIALIZE
  zsig.initialize();
  #endif
  csig.allocate("csig");
  #ifndef NO_AD_INITIALIZE
  csig.initialize();
  #endif
  sel.allocate(1,ages,"sel");
  #ifndef NO_AD_INITIALIZE
    sel.initialize();
  #endif
  Mage.allocate(1,ages,"Mage");
  #ifndef NO_AD_INITIALIZE
    Mage.initialize();
  #endif
  lxo.allocate(1,ages,"lxo");
  #ifndef NO_AD_INITIALIZE
    lxo.initialize();
  #endif
  ralpha.allocate("ralpha");
  #ifndef NO_AD_INITIALIZE
  ralpha.initialize();
  #endif
  rbeta.allocate("rbeta");
  #ifndef NO_AD_INITIALIZE
  rbeta.initialize();
  #endif
  phie.allocate("phie");
  #ifndef NO_AD_INITIALIZE
  phie.initialize();
  #endif
  Nt.allocate(1,yrs,1,ages,"Nt");
  #ifndef NO_AD_INITIALIZE
    Nt.initialize();
  #endif
  log_Rt.allocate(1,yrs,"log_Rt");
  #ifndef NO_AD_INITIALIZE
    log_Rt.initialize();
  #endif
  Rt_bar.allocate(2,yrs,"Rt_bar");
  #ifndef NO_AD_INITIALIZE
    Rt_bar.initialize();
  #endif
  rdelta.allocate(2,yrs,"rdelta");
  #ifndef NO_AD_INITIALIZE
    rdelta.initialize();
  #endif
  Bt.allocate(1,yrs,"Bt");
  #ifndef NO_AD_INITIALIZE
    Bt.initialize();
  #endif
  Zt.allocate(1,yrs,"Zt");
  #ifndef NO_AD_INITIALIZE
    Zt.initialize();
  #endif
  plandings.allocate(1,yrs,"plandings");
  #ifndef NO_AD_INITIALIZE
    plandings.initialize();
  #endif
  SSBt.allocate(1,yrs,"SSBt");
  #ifndef NO_AD_INITIALIZE
    SSBt.initialize();
  #endif
  Fat.allocate(1,yrs,1,ages,"Fat");
  #ifndef NO_AD_INITIALIZE
    Fat.initialize();
  #endif
  Zat.allocate(1,yrs,1,ages,"Zat");
  #ifndef NO_AD_INITIALIZE
    Zat.initialize();
  #endif
  pcat.allocate(1,yrs,1,ages,"pcat");
  #ifndef NO_AD_INITIALIZE
    pcat.initialize();
  #endif
  ppat.allocate(1,yrs,1,ages,"ppat");
  #ifndef NO_AD_INITIALIZE
    ppat.initialize();
  #endif
  cdev.allocate(1,yrs,"cdev");
  #ifndef NO_AD_INITIALIZE
    cdev.initialize();
  #endif
  zdev.allocate(1,yrs,"zdev");
  #ifndef NO_AD_INITIALIZE
    zdev.initialize();
  #endif
  pdev.allocate(1,yrs,1,ages,"pdev");
  #ifndef NO_AD_INITIALIZE
    pdev.initialize();
  #endif
  tau.allocate("tau");
  #ifndef NO_AD_INITIALIZE
  tau.initialize();
  #endif
  sig.allocate("sig");
  #ifndef NO_AD_INITIALIZE
  sig.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  f =0.0;
	//*** Main function calls  ***//
		initial_calculations();
		recursive_calculations();
		recruitment_calculations();
		observation_model();
		calc_objective_function();
		if(mceval_phase()) mcmc_stuff();		
	//****************************//
}

void model_parameters::initial_calculations(void)
{
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
}

void model_parameters::recursive_calculations(void)
{
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
}

void model_parameters::recruitment_calculations(void)
{
	{
		rdelta=log(Rt_bar)-log_Rt(2,yrs);//+0.5*tau*tau;
	}
}

void model_parameters::observation_model(void)
{
	{
		sig=rho/varphi;
		cdev=log(olandings)-log(plandings);
		pdev=log(opat)-log(ppat);
		zdev=Zt-mean(Zt);
	}
}

void model_parameters::calc_objective_function(void)
{
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
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
}

void model_parameters::final_calcs()
{
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
}

void model_parameters::mcmc_stuff(void)
{
		adstring fl1=adstring("LMCMC.rep");
	 	if(iter==0){
			//ofstream ofs1(fl1);
			//ofs1<<"HEADER STUFF\t"<<endl;
		}
		iter++;
		ofstream ofs1(fl1,ios::app);	
		//ofs1<<"dump mcmc results"<<endl;
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
