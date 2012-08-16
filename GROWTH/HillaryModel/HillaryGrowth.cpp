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
	
#include <admodel.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <HillaryGrowth.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
		int on,opt;
		method=1;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-method",opt))>-1)
		{
			method=atoi(ad_comm::argv[on+1]);
		}
		cout<<method<<endl;
		if(method !=1 && method !=2)
		{
			cout<<"Available methods are \n \t Fabens = 1 \n \t Zhang = 2\n";
			exit(1);
		}
 nobs = 100;  //total number of records is 8315
  gdata.allocate(1,nobs,1,3,"gdata");
  ltag.allocate(1,nobs);
  dl.allocate(1,nobs);
  dt.allocate(1,nobs);
 ltag = column(gdata,1);
 dl   = column(gdata,2);
 dt   = column(gdata,3);	
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  log_linf.allocate("log_linf");
  log_k.allocate("log_k");
  log_phi.allocate("log_phi");
 log_linf = log(mean(ltag+dl));
 log_k    = log(0.15);
 log_phi  = log(50);
 if(method==2) phz=2; else phz=-1;
  log_sig_linf.allocate(phz,"log_sig_linf");
  delta.allocate(1,nobs,phz,"delta");
  f.allocate("f");  /* ADOBJECTIVEFUNCTION */
  linf.allocate("linf");
  #ifndef NO_AD_INITIALIZE
  linf.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  phi.allocate("phi");
  #ifndef NO_AD_INITIALIZE
  phi.initialize();
  #endif
  sig_linf.allocate("sig_linf");
  #ifndef NO_AD_INITIALIZE
  sig_linf.initialize();
  #endif
  dl_hat.allocate(1,nobs,"dl_hat");
  #ifndef NO_AD_INITIALIZE
    dl_hat.initialize();
  #endif
  epsilon.allocate(1,nobs,"epsilon");
  #ifndef NO_AD_INITIALIZE
    epsilon.initialize();
  #endif
  sd_dl1.allocate("sd_dl1");
}
void model_parameters::userfunction(void)
{
	
	if(method==1) fabens_method();
	if(method==2) zhang_method();
	sd_dl1 = dl_hat(1);
	
}

void model_parameters::zhang_method(void)
{
	linf = mfexp(log_linf);
	k    = mfexp(log_k);
	phi  = mfexp(log_phi);
	sig_linf = mfexp(log_sig_linf);
	
	dvar_vector linf_i=linf*exp(delta * sig_linf);
	dl_hat = elem_prod(linf_i-ltag,1.-exp(-k*dt));
	epsilon = dl - dl_hat;
	dvar_vector sig = sqrt(31.8 + square(phi*ltag));
	f = dnorm(epsilon,sig);
	f+= dnorm(delta,1);
	
}

void model_parameters::fabens_method(void)
{
	linf = mfexp(log_linf);
	k    = mfexp(log_k);
	phi  = mfexp(log_phi);
	
	dl_hat = elem_prod(linf-ltag,1.-exp(-k*dt));
	epsilon = dl - dl_hat;
	dvar_vector df = 1./(1./100. + ltag/phi);
	//cout<<epsilon(1,1)<<" "<<df(1,1)<<" "<<dstudent_t(epsilon(1,1),df(1,1))<<endl;
	//f = dstudent_t(epsilon,df);
	
	//sigm^2 from Coggins 2006 was 31.8mm
	dvar_vector sig = sqrt(31.8 + square(phi*ltag));
	f = dnorm(epsilon,sig);
	//f= norm2(epsilon);
	
}

dvariable model_parameters::dstudent_t( const dvar_vector& residual, const dvar_vector& df)
{
	{
		double pi =  3.141593;
		dvar_vector t1 = 0.5*(df+1);
		dvar_vector t2 = gammln(t1);
		dvar_vector t3 = 0.5*log(df*pi)+gammln(0.5*df);
		dvar_vector t4 = elem_prod(t1,log(1+elem_div(square(residual),df)));
		dvariable pdf = sum(t3+t4-t2);
		return( pdf );
	}
}

dvariable model_parameters::dnorm(const dvar_vector& residual, const double std)
{
	{
		double pi=3.141593;
		long n=size_count(residual);
		dvariable SS=norm2(residual);
		return(n*(0.5*log(2.*pi)+log(std))+0.5*SS/(std*std));
		//return(n*(0.5*log(2.*pi)-log(std))+0.5*SS*(std*std));
	}
}

dvariable model_parameters::dnorm( const dvar_vector& residual, const dvar_vector std )
{
	{
		RETURN_ARRAYS_INCREMENT();
		int n=size_count(residual);
		double pi =  3.141593;
		dvar_vector var=elem_prod(std,std);
		dvar_vector SS=elem_prod(residual,residual);
		dvariable tmp=0.5*n*log(2.*pi)+sum(log(std))+sum(elem_div(SS,2.*var));
		RETURN_ARRAYS_DECREMENT();
		return( tmp );
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
	REPORT(linf);
	REPORT(k);
	REPORT(phi);
	REPORT(epsilon);
	REPORT(ltag);
	REPORT(dl);
	REPORT(dt);
	REPORT(dl_hat);
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
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
	gradient_structure::set_MAX_NVAR_OFFSET(250504);
 
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
	
	if(method==1) fabens_method();
	if(method==2) zhang_method();
	sd_dl1 = dl_hat(1);
	
}

void df1b2_parameters::zhang_method(void)
{
	linf = mfexp(log_linf);
	k    = mfexp(log_k);
	phi  = mfexp(log_phi);
	sig_linf = mfexp(log_sig_linf);
	
	df1b2vector linf_i=linf*exp(delta * sig_linf);
	dl_hat = elem_prod(linf_i-ltag,1.-exp(-k*dt));
	epsilon = dl - dl_hat;
	df1b2vector sig = sqrt(31.8 + square(phi*ltag));
	f = dnorm(epsilon,sig);
	f+= dnorm(delta,1);
	
}

void df1b2_parameters::fabens_method(void)
{
	linf = mfexp(log_linf);
	k    = mfexp(log_k);
	phi  = mfexp(log_phi);
	
	dl_hat = elem_prod(linf-ltag,1.-exp(-k*dt));
	epsilon = dl - dl_hat;
	df1b2vector df = 1./(1./100. + ltag/phi);
	//cout<<epsilon(1,1)<<" "<<df(1,1)<<" "<<dstudent_t(epsilon(1,1),df(1,1))<<endl;
	//f = dstudent_t(epsilon,df);
	
	//sigm^2 from Coggins 2006 was 31.8mm
	df1b2vector sig = sqrt(31.8 + square(phi*ltag));
	f = dnorm(epsilon,sig);
	//f= norm2(epsilon);
	
}

df1b2variable df1b2_parameters::dstudent_t( const df1b2vector& residual, const df1b2vector& df)
{
	{
		double pi =  3.141593;
		df1b2vector t1 = 0.5*(df+1);
		df1b2vector t2 = gammln(t1);
		df1b2vector t3 = 0.5*log(df*pi)+gammln(0.5*df);
		df1b2vector t4 = elem_prod(t1,log(1+elem_div(square(residual),df)));
		df1b2variable pdf = sum(t3+t4-t2);
		return( pdf );
	}
}

df1b2variable df1b2_parameters::dnorm(const df1b2vector& residual, const double std)
{
	{
		double pi=3.141593;
		long n=size_count(residual);
		df1b2variable SS=norm2(residual);
		return(n*(0.5*log(2.*pi)+log(std))+0.5*SS/(std*std));
		//return(n*(0.5*log(2.*pi)-log(std))+0.5*SS*(std*std));
	}
}

df1b2variable df1b2_parameters::dnorm( const df1b2vector& residual, const df1b2vector std )
{
	{
		RETURN_ARRAYS_INCREMENT();
		int n=size_count(residual);
		double pi =  3.141593;
		df1b2vector var=elem_prod(std,std);
		df1b2vector SS=elem_prod(residual,residual);
		df1b2variable tmp=0.5*n*log(2.*pi)+sum(log(std))+sum(elem_div(SS,2.*var));
		RETURN_ARRAYS_DECREMENT();
		return( tmp );
	}
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
  log_linf.allocate("log_linf");
  log_k.allocate("log_k");
  log_phi.allocate("log_phi");
 log_linf = log(mean(ltag+dl));
 log_k    = log(0.15);
 log_phi  = log(50);
 if(method==2) phz=2; else phz=-1;
  log_sig_linf.allocate(phz,"log_sig_linf");
  delta.allocate(1,nobs,phz,"delta");
  f.allocate("f");  /* ADOBJECTIVEFUNCTION */
  linf.allocate("linf");
  #ifndef NO_AD_INITIALIZE
  linf.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  phi.allocate("phi");
  #ifndef NO_AD_INITIALIZE
  phi.initialize();
  #endif
  sig_linf.allocate("sig_linf");
  #ifndef NO_AD_INITIALIZE
  sig_linf.initialize();
  #endif
  dl_hat.allocate(1,nobs,"dl_hat");
  #ifndef NO_AD_INITIALIZE
    dl_hat.initialize();
  #endif
  epsilon.allocate(1,nobs,"epsilon");
  #ifndef NO_AD_INITIALIZE
    epsilon.initialize();
  #endif
  sd_dl1.allocate("sd_dl1");
}
