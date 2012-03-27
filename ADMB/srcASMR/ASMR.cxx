/*

*/

#include<admodel.h>

//normal distribution
dvariable dnorm(const dvariable& x, const double& mu, const double& std)
{
	double pi=3.141593;
	return 0.5*log(2.*pi)+log(std)+0.5*square(x-mu)/(std*std);
}


//plogis
double plogis(const double& x, const double& mu, const double& std)
{
	return 1./(1.+mfexp((mu-x)/std));
}

dvariable plogis(const dvariable& x, const double& mu, const dvariable& std)
{
	return 1./(1.+mfexp((mu-x)/std));
}

dvar_vector plogis(const dvector& x, const dvariable& mu, const dvariable& std)
{
	return 1./(1.+mfexp((mu-x)/std));
}

dvector plogis(const dvector& x, const double& mu, const double& std)
{
	return 1./(1.+mfexp((mu-x)/std));
}

dvar_vector plogis(const dvar_vector& x, const dvariable& mu, const dvariable& std)
{
	return 1./(1.+mfexp((mu-x)/std));
}

dvar_vector eplogis(const dvar_vector& x, const dvariable& alpha, const dvariable& beta, const dvariable& gamma)
{
	//exponential logistic based on Grant Thompson (1994) Paper, CJFAS.
	return (1./(1.-gamma))*pow((1.-gamma)/gamma,gamma)*elem_div(exp(alpha*gamma*(beta-x)),1.+exp(alpha*(beta-x)));
}

dvector eplogis(const dvector& x, const double& alpha, const double& beta, const double& gamma)
{
	//exponential logistic based on Grant Thompson (1994) Paper, CJFAS.
	return (1./(1.-gamma))*pow((1.-gamma)/gamma,gamma)*elem_div(exp(alpha*gamma*(beta-x)),1.+exp(alpha*(beta-x)));
}

//log normal distribution
dvariable dlnorm(const dvariable& x, const double& mu, const double& std)
{
	double pi=3.141593;
	return 0.5*log(2.*pi)+log(std)+log(x)+square(log(x)-mu)/(2.*std*std);
}


dvariable dpois(const double& k, const dvariable& lambda)
{
	return -k*log(lambda)+lambda + gammln(k+1.);
}

dvariable dpois(const dvector& k, const dvar_vector& lambda)
{
	/*	A modification to the poisson distribution
			where we first loop over observations k and
			get the indexes for non-zero values.
	*/
	int i,imin,imax;
	imin=k.indexmin();
	imax=k.indexmax();
	dvariable loglike=0.;
	for(i = imin; i<=imax;i++)
		if(k(i)>0) loglike += lambda(i)-k(i)*log(lambda(i));
	//return sum(lambda - elem_prod(k,log(lambda)));
	return loglike;
}

dvariable dnbinom(const dvector& x, const dvar_vector& lambda, const dvariable& tau, dvector& residual)
{
	//the observed counts are in x
	//lambda is the predicted count
	//tau is the overdispersion parameter
	RETURN_ARRAYS_INCREMENT();
	int i,imin,imax;
	double o=1e-30;
	imin=x.indexmin();
	imax=x.indexmax();
	dvariable loglike=0.;
	residual.initialize();

	for(i = imin; i<=imax; i++)
		if(x(i)>0){
			dvariable p=tau/(tau+lambda(i));
			loglike += gammln(tau+x(i)) - gammln(tau) -factln(x(i))//- gammln(x(i)+1)
								+tau*log(p) + x(i)*log(1.-p);
			residual(i) =value((x(i)-lambda(i))/sqrt(o+lambda(i)+square(lambda(i))/tau)) ;
		}
	//cout<<"OK in dnbinom"<<endl;
	RETURN_ARRAYS_DECREMENT();
	return(-1.0 *loglike);
}

dvariable dpois_residual(const dvector& k, const dvar_vector& lambda, dvector& residual)
{
	/*	A modification to the poisson distribution
			where we first loop over observations k and
			get the indexes for non-zero values.
			
			This was not done in the original ASMR model
			so I have commented it out for now, and added
			a tiny number to lambda.
	*/
	RETURN_ARRAYS_INCREMENT();
	double o=1.e-30;
	int i,imin,imax;
	imin=k.indexmin();
	imax=k.indexmax();
	dvariable loglike=0.;
	for(i = imin; i<=imax;i++)
	{
		//if(k(i)>0) {
			loglike += lambda(i)-k(i)*log(lambda(i)+o);
			//cout<<k(i)<<"\t"<<lambda(i)<<endl;
			residual(i) = value((k(i)-lambda(i))/sqrt(lambda(i)+o));
		//}
	}
	//return sum(lambda - elem_prod(k,log(lambda)));
	RETURN_ARRAYS_DECREMENT();
	return loglike;
}

dvar_matrix ALK(dvar_vector mu, dvar_vector sig, dvector x)
{
	//This function returns an Age-Length Key
	int i, j;
	dvariable z1;
	dvariable z2;
	int si,ni; si=mu.indexmin(); ni=mu.indexmax();
	int sj,nj; sj=x.indexmin(); nj=x.indexmax();
	dvar_matrix pdf(si,ni,sj,nj);
	pdf.initialize();
	double xs=0.5*(x[sj+1]-x[sj]);
	for(i=si;i<=ni;i++) //loop over ages
	{
		 for(j=sj;j<=nj;j++) //loop over length bins
		{
			z1=((x(j)-xs)-mu(i))/sig(i);
			z2=((x(j)+xs)-mu(i))/sig(i);
			pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
		//pdf(i)/=sum(pdf(i));
	}//end nage
	pdf/=sum(pdf);
	return(pdf);
}

dmatrix ALK(dvector mu, dvector sig, dvector x)
{
	//This function returns an Age-Length Key
	int i, j;
	double z1;
	double z2;
	int si,ni; si=mu.indexmin(); ni=mu.indexmax();
	int sj,nj; sj=x.indexmin(); nj=x.indexmax();
	dmatrix pdf(si,ni,sj,nj);
	pdf.initialize();
	double xs=0.5*(x[sj+1]-x[sj]);
	for(i=si;i<=ni;i++) //loop over ages
	{
		 for(j=sj;j<=nj;j++) //loop over length bins
		{
			z1=((x(j)-xs)-mu(i))/sig(i);
			z2=((x(j)+xs)-mu(i))/sig(i);
			pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
		//pdf(i)/=sum(pdf(i));
	}//end nage
	pdf/=sum(pdf);
	return(pdf);
}

