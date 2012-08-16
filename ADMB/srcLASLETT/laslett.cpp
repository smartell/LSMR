	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
	#include <admodel.h>
	#include <time.h>
	//#include <statsLib.h>
	
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
#include <admodel.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <laslett.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
 ad_comm::change_datafile_name("HBCdata.dat");
		int on,opt;
		if((on = option_match(ad_comm::argc,ad_comm::argv,"-iyr",opt))>-1)
		{
			iyr = atoi(ad_comm::argv[on+1]);
		}
		else
		{
			iyr = 1989;
		}
  nobs.allocate("nobs");
  data.allocate(1,nobs,1,6,"data");
		int i;
		i1=0;
		i2=0;
		for(i=1;i<=nobs;i++)
		{
			if(data(i,1)==iyr && i1==0) i1 = i;
			if(data(i,1)==iyr+1) break;
			i2 = i;
		}
  ii.allocate(i1,i2);
 ii.fill_seqadd(i1,1);
  l1.allocate(i1,i2);
  l2.allocate(i1,i2);
  dt.allocate(i1,i2);
		for(i=i1;i<=i2;i++)
		{
			l1(i) = data(i,3);
			l2(i) = data(i,4);
			dt(i) = data(i,5)/365.25;
		}
}

void model_parameters::initializationfunction(void)
{
  k.set_initial_value(0.2);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  mu_linf.allocate(50,550,1,"mu_linf");
  sd_linf.allocate(0,100,2,"sd_linf");
  k.allocate(0,10,1,"k");
  sig.allocate(0,150,1,"sig");
  muA.allocate(0,50,1,"muA");
  sigA.allocate(0,25,2,"sigA");
  A.allocate(i1,i2,2,"A");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  f.allocate("f");  /* ADOBJECTIVEFUNCTION */
  sig2.allocate("sig2");
  #ifndef NO_AD_INITIALIZE
  sig2.initialize();
  #endif
  sd2_linf.allocate("sd2_linf");
  age.allocate(i1,i2,"age");
  #ifndef NO_AD_INITIALIZE
    age.initialize();
  #endif
  l1_hat.allocate(i1,i2,"l1_hat");
  #ifndef NO_AD_INITIALIZE
    l1_hat.initialize();
  #endif
  l2_hat.allocate(i1,i2,"l2_hat");
  #ifndef NO_AD_INITIALIZE
    l2_hat.initialize();
  #endif
  linf.allocate(i1,i2,"linf");
  #ifndef NO_AD_INITIALIZE
    linf.initialize();
  #endif
  f1.allocate(i1,i2,"f1");
  #ifndef NO_AD_INITIALIZE
    f1.initialize();
  #endif
  f2.allocate(i1,i2,"f2");
  #ifndef NO_AD_INITIALIZE
    f2.initialize();
  #endif
  rho.allocate(i1,i2,"rho");
  #ifndef NO_AD_INITIALIZE
    rho.initialize();
  #endif
  nll.allocate(i1,i2,"nll");
  #ifndef NO_AD_INITIALIZE
    nll.initialize();
  #endif
}
void model_parameters::userfunction(void)
{
  f =0.0;
	calc_relative_age();
	calc_predicted_length();
	calc_objective_function();
	
}

void model_parameters::calc_relative_age(void)
{
	int i;
	dvariable d,tmp;
	for(i=i1;i<=i2;i++)
	{
		//d      = mfexp(-k*dt(i));
		//tmp    = (l2(i)-l1(i))/(l2(i)-d*l1(i));
		//age(i) = -1./k*log(tmp)*exp(A(i)*sigA);  //add random effect here
		//tmp      = -log(fabs(mu_linf-l1(i))/mu_linf) / k;
		//age(i)   = tmp*exp(A(i)*sigA);
		age(i) = muA*exp(A(i)*sigA);
	}
}

void model_parameters::calc_predicted_length(void)
{
	int i;
	dvariable t1,t3;
	sig2     = square(sig);
	sd2_linf = square(sd_linf);
	for(i=i1;i<=i2;i++)
	{
		f1(i)      = (1.-exp(-k*age(i)));
		f2(i)      = (1.-exp(-k*(age(i)+dt(i))));
		t1         = sd2_linf/(sig2+sd2_linf*(square(f1(i))+square(f2(i))));
		t3         = f1(i)*(l1(i)-mu_linf*f1(i))+f2(i)*(l2(i)-mu_linf*f2(i));
		linf(i)    = mu_linf+t1*t3;  
		l1_hat(i)  = mu_linf*f1(i);
		l2_hat(i)  = mu_linf*f2(i);
	}
	
}

void model_parameters::calc_objective_function(void)
{
	int i;
	dvariable s1,s2,t1,t2,t3,t9;
	dvariable r1,r2;
	for(i=i1;i<=i2;i++)
	{
		s1     = sqrt(sd2_linf*square(f1(i)) + sig2);
		s2     = sqrt(sd2_linf*square(f2(i)) + sig2);
		rho(i) = sd2_linf*f1(i)*f2(i)/(s1*s2);
		
		r1     = l1(i)-l1_hat(i);
		r2     = l2(i)-l2_hat(i);
		
		t1     = square(r1)/(s1*s1);
		t2     = 2.*rho(i) * (r1*r2) / (s1*s2);
		t3     = square(r2)/(s2*s2);
		t9     = 1.-square(rho(i));
		nll(i) = log(2.*M_PI)+log(s1)+log(s2)
		         +0.5*log(t9)
		         +(t1-t2+t3)/(2.*t9);
	}
	
	f = sum(nll) + 0.5*norm2(A);
	
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
	REPORT(f);
	REPORT(mu_linf);
	REPORT(sd_linf);
	REPORT(k);
	REPORT(sig);
	REPORT(muA);
	REPORT(sigA);
	REPORT(rho);
	REPORT(nll);
	REPORT(ii)
	REPORT(age);
	REPORT(dt);
	REPORT(linf);
	REPORT(l1);
	REPORT(l1_hat);
	REPORT(l2);
	REPORT(l2_hat);
	dvector l1_res = value((l1 - elem_prod(linf,f1))/sig);
	dvector l2_res = value((l2 - elem_prod(linf,f2))/sig);
	REPORT(l1_res);
	REPORT(l2_res);
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
  long int arrmblsize=0;

int main(int argc,char * argv[])
{
  ad_set_new_handler();
  ad_exit=&ad_boundf;
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(1600000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(1600000);
 
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    df1b2variable::noallocate=1;
    df1b2_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;

    function_minimizer::random_effects_flag=1;
    df1b2variable::noallocate=0;
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

void df1b2_parameters::user_function(void)
{
  f =0.0;
	calc_relative_age();
	calc_predicted_length();
	calc_objective_function();
	
}

void df1b2_parameters::calc_relative_age(void)
{
	int i;
	df1b2variable d,tmp;
	for(i=i1;i<=i2;i++)
	{
		//d      = mfexp(-k*dt(i));
		//tmp    = (l2(i)-l1(i))/(l2(i)-d*l1(i));
		//age(i) = -1./k*log(tmp)*exp(A(i)*sigA);  //add random effect here
		//tmp      = -log(fabs(mu_linf-l1(i))/mu_linf) / k;
		//age(i)   = tmp*exp(A(i)*sigA);
		age(i) = muA*exp(A(i)*sigA);
	}
}

void df1b2_parameters::calc_predicted_length(void)
{
	int i;
	df1b2variable t1,t3;
	sig2     = square(sig);
	sd2_linf = square(sd_linf);
	for(i=i1;i<=i2;i++)
	{
		f1(i)      = (1.-exp(-k*age(i)));
		f2(i)      = (1.-exp(-k*(age(i)+dt(i))));
		t1         = sd2_linf/(sig2+sd2_linf*(square(f1(i))+square(f2(i))));
		t3         = f1(i)*(l1(i)-mu_linf*f1(i))+f2(i)*(l2(i)-mu_linf*f2(i));
		linf(i)    = mu_linf+t1*t3;  
		l1_hat(i)  = mu_linf*f1(i);
		l2_hat(i)  = mu_linf*f2(i);
	}
	
}

void df1b2_parameters::calc_objective_function(void)
{
	int i;
	df1b2variable s1,s2,t1,t2,t3,t9;
	df1b2variable r1,r2;
	for(i=i1;i<=i2;i++)
	{
		s1     = sqrt(sd2_linf*square(f1(i)) + sig2);
		s2     = sqrt(sd2_linf*square(f2(i)) + sig2);
		rho(i) = sd2_linf*f1(i)*f2(i)/(s1*s2);
		
		r1     = l1(i)-l1_hat(i);
		r2     = l2(i)-l2_hat(i);
		
		t1     = square(r1)/(s1*s1);
		t2     = 2.*rho(i) * (r1*r2) / (s1*s2);
		t3     = square(r2)/(s2*s2);
		t9     = 1.-square(rho(i));
		nll(i) = log(2.*M_PI)+log(s1)+log(s2)
		         +0.5*log(t9)
		         +(t1-t2+t3)/(2.*t9);
	}
	
	f = sum(nll) + 0.5*norm2(A);
	
}

void df1b2_pre_parameters::setup_quadprior_calcs(void) 
{ 
df1b2_gradlist::set_no_derivatives(); 
quadratic_prior::in_qp_calculations=1; 
}  

void df1b2_pre_parameters::begin_df1b2_funnel(void) 
{ 
(*re_objective_function_value::pobjfun)=0; 
other_separable_stuff_begin(); 
f1b2gradlist->reset();  
if (!quadratic_prior::in_qp_calculations) 
{ 
df1b2_gradlist::set_yes_derivatives();  
} 
funnel_init_var::allocate_all();  
}  

void df1b2_pre_parameters::end_df1b2_funnel(void) 
{  
lapprox->do_separable_stuff(); 
other_separable_stuff_end(); 
} 

void model_parameters::begin_df1b2_funnel(void) 
{ 
if (lapprox)  
{  
{  
begin_funnel_stuff();  
}  
}  
}  

void model_parameters::end_df1b2_funnel(void) 
{  
if (lapprox)  
{  
end_df1b2_funnel_stuff();  
}  
} 
void df1b2_parameters::allocate(void) 
{
  mu_linf.allocate(50,550,1,"mu_linf");
  sd_linf.allocate(0,100,2,"sd_linf");
  k.allocate(0,10,1,"k");
  sig.allocate(0,150,1,"sig");
  muA.allocate(0,50,1,"muA");
  sigA.allocate(0,25,2,"sigA");
  A.allocate(i1,i2,2,"A");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  f.allocate("f");  /* ADOBJECTIVEFUNCTION */
  sig2.allocate("sig2");
  #ifndef NO_AD_INITIALIZE
  sig2.initialize();
  #endif
  sd2_linf.allocate("sd2_linf");
  age.allocate(i1,i2,"age");
  #ifndef NO_AD_INITIALIZE
    age.initialize();
  #endif
  l1_hat.allocate(i1,i2,"l1_hat");
  #ifndef NO_AD_INITIALIZE
    l1_hat.initialize();
  #endif
  l2_hat.allocate(i1,i2,"l2_hat");
  #ifndef NO_AD_INITIALIZE
    l2_hat.initialize();
  #endif
  linf.allocate(i1,i2,"linf");
  #ifndef NO_AD_INITIALIZE
    linf.initialize();
  #endif
  f1.allocate(i1,i2,"f1");
  #ifndef NO_AD_INITIALIZE
    f1.initialize();
  #endif
  f2.allocate(i1,i2,"f2");
  #ifndef NO_AD_INITIALIZE
    f2.initialize();
  #endif
  rho.allocate(i1,i2,"rho");
  #ifndef NO_AD_INITIALIZE
    rho.initialize();
  #endif
  nll.allocate(i1,i2,"nll");
  #ifndef NO_AD_INITIALIZE
    nll.initialize();
  #endif
}
