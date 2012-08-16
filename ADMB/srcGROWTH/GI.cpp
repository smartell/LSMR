	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
	#include <admodel.h>
	#include <time.h>
	#include <statsLib.h>
	
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <GI.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  str_dataFile.allocate("str_dataFile");
ad_comm::change_datafile_name(str_dataFile);
  nobs.allocate("nobs");
  inobs.allocate(1,nobs,"inobs");
  data.allocate(1,nobs,1,inobs,1,6,"data");
  eof.allocate("eof");
 if(eof != 999){cout<<"Error reading data\n"; exit(1);}
  iyr.allocate(1,nobs,1,inobs);
  loc.allocate(1,nobs,1,inobs);
  l1.allocate(1,nobs,1,inobs);
  l2.allocate(1,nobs,1,inobs);
  dt.allocate(1,nobs,1,inobs);
  dl.allocate(1,nobs,1,inobs);
		int i;
		for(i=1;i<=nobs;i++)
		{
			iyr(i) = ivector(column(data(i),1));
			loc(i) = ivector(column(data(i),2));
			l1(i)  = column(data(i),3);
			l2(i)  = column(data(i),4);
			dt(i)  = column(data(i),5)/365.25;
			dl(i)  = column(data(i),6);
		}
		//cout<<iyr<<endl;
 nbins = 46;
  lbins.allocate(1,nbins);
 lbins.fill_seqadd(50,10);
 cout<<lbins<<endl;
  tmat.allocate(1,nobs,1,nbins-1,1,nbins-1);
  T.allocate(1,nobs,1,nbins-1,1,nbins-1);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  mu_log_linf.allocate("mu_log_linf");
  mu_log_k.allocate("mu_log_k");
  mu_log_phi.allocate("mu_log_phi");
  sig_log_linf.allocate(3,"sig_log_linf");
  sig_log_k.allocate(3,"sig_log_k");
  sig_log_phi.allocate(3,"sig_log_phi");
  log_linf.allocate(1,nobs,-5,5,2,"log_linf");
  log_k.allocate(1,nobs,-5,5,2,"log_k");
  log_phi.allocate(1,nobs,-5,5,2,"log_phi");
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  linf.allocate(1,nobs,"linf");
  #ifndef NO_AD_INITIALIZE
    linf.initialize();
  #endif
  k.allocate(1,nobs,"k");
  #ifndef NO_AD_INITIALIZE
    k.initialize();
  #endif
  phi.allocate(1,nobs,"phi");
  #ifndef NO_AD_INITIALIZE
    phi.initialize();
  #endif
  fvec.allocate(1,nobs,"fvec");
  #ifndef NO_AD_INITIALIZE
    fvec.initialize();
  #endif
  dl_hat.allocate(1,nobs,1,inobs,"dl_hat");
  #ifndef NO_AD_INITIALIZE
    dl_hat.initialize();
  #endif
  epsilon.allocate(1,nobs,1,inobs,"epsilon");
  #ifndef NO_AD_INITIALIZE
    epsilon.initialize();
  #endif
  sig.allocate(1,nobs,1,inobs,"sig");
  #ifndef NO_AD_INITIALIZE
    sig.initialize();
  #endif
  df.allocate(1,nobs,1,inobs,"df");
  #ifndef NO_AD_INITIALIZE
    df.initialize();
  #endif
  sd_mu_linf.allocate("sd_mu_linf");
}

void model_parameters::initializationfunction(void)
{
  mu_log_linf.set_initial_value(5.0);
  mu_log_k.set_initial_value(-2.0);
  mu_log_phi.set_initial_value(-2.5);
  sig_log_linf.set_initial_value(-2.3);
  sig_log_k.set_initial_value(-2.3);
  sig_log_phi.set_initial_value(-2.3);
}

void model_parameters::userfunction(void)
{
  f =0.0;
	vbginc();
	calc_objective_function();
	if(mceval_phase())
	{
		calcLebesgueMeasure();
	} 
}

void model_parameters::calcLebesgueMeasure(void)
{
	/*
	This function is based on Rich Hillary's algorithm for calculating
	the relative overlap of each length partition based on parametric
	uncertainty and measurement error.  In this case, I have used a normal
	density function for measurement error.
	*/
	int h,i,j;
	static int nf = 0;
	double tau    = 1.0;         //Time interval
	double i1,i2,sig,eps;
	double j1,j2,l1i,l2i;
	double mu,nu;
	dvector tmp(1,2);
	if(nf == 0) 
	{
		tmat.initialize();
		T.initialize();
		cout<<"Initialized Arrays"<<endl;
	}
	nf++;
	for(h=1;h<=nobs;h++)
	{
		double ilinf = value(linf(h));
		double ik    = value(k(h));
		for(i=1;i < nbins;i++)
		{
			i1  = lbins(i);
			i2  = lbins(i+1);
			sig = value(phi(h))*0.5*(i1+i2);
			eps = randn(nf+i)*sig;
			for(j=1;j < nbins;j++)
			{
				j1    = lbins(j);
				j2    = lbins(j+1);
				l1i   = i1 + ginc(ilinf,ik,tau,i1) + eps;
				l2i   = i2 + ginc(ilinf,ik,tau,i2) + eps;
				/* Calculate Lebesque measure of intersection */
				if(l1i>j1 && l2i<j2)
				{
					tmat(h,i,j) += 1.0;
				}
				else
				{
					tmp(1)        = max(l1i,j1);
					tmp(2)        = min(l2i,j2);
					if(tmp(1)<tmp(2)) mu = tmp(2)-tmp(1); else mu = 0.0;
					nu            = l2i-l1i;
					tmat(h,i,j) += mu/nu;
				}
			}
		}
		T(h) = tmat(h)/nf;
	}
	if(nf==1000) cout<<"Transition Matrix\n"<<T(1)<<endl;
}

double model_parameters::ginc(const double& linf, const double& k, const double& tau,const double& l1)
{
	return( (linf-l1)*(1.-exp(-k*tau)) );
}

void model_parameters::vbginc(void)
{
	linf    = mfexp(mu_log_linf + log_linf);
	k       = mfexp(mu_log_k + log_k);
	phi     = mfexp(mu_log_phi + log_phi);
	sd_mu_linf = mean(linf);
	int i;
	for(i=1;i<=nobs;i++)
	{
		dl_hat(i)  = elem_div(elem_prod((linf(i)-l1(i)),(1.-exp(-k(i)*dt(i)))),dt(i));
		epsilon(i) = dl(i) - dl_hat(i);
		sig(i)     = sqrt( 88.6 + square(phi(i)*l1(i)) );
		df(i)      = 1./(1./88.6 + l1(i)/phi(i));
	}
}

dvariable model_parameters::deprecated_dstudent_t( const dvar_vector& residual, const dvar_vector& df)
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

void model_parameters::calc_objective_function(void)
{
	int i;
	fvec.initialize();
	for(i=1;i<=nobs;i++)
	{
		//fvec(i) = dstudent_t(epsilon(i),df(i));
		fvec(i) = dnorm(epsilon(i),sig(i));
	}
	f  = sum(fvec);
	f += dnorm(log_k,mfexp(sig_log_k));
	f += dnorm(log_linf,mfexp(sig_log_linf));
	f += dnorm(log_phi,mfexp(sig_log_phi));
	f += dgamma(1/mfexp(sig_log_k),1.01,1.01);
	f += dgamma(1/mfexp(sig_log_linf),1.01,1.01);
	f += dgamma(1/mfexp(sig_log_phi),1.01,1.01);
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
	REPORT(nobs);
	REPORT(linf);
	REPORT(k);
	REPORT(phi);
	REPORT(l1);
	REPORT(dl);
	REPORT(dl_hat);
	REPORT(dt);
	REPORT(epsilon);
	REPORT(sig);
	REPORT(df);
	dmatrix nu = value(elem_div(df,df-2.));
	REPORT(nu);
	REPORT(loc);
	REPORT(iyr);
}

void model_parameters::final_calcs()
{
	ofstream ofs("TransitionMatrix.txt");
	int i;
	ofs<<"lbins\n"<<lbins<<endl;
	ofs<<"T"<<endl;
	for(i=1;i<=nobs;i++)
	{
		ofs<<"#"<<iyr(i,1)<<endl<<T(i)<<endl;
	}
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
