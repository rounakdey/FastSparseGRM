	#include<RcppArmadillo.h>
	//[[Rcpp::depends(RcppArmadillo)]]
	// [[Rcpp::plugins(cpp11)]]

	using namespace Rcpp;
	using namespace arma;
    	using namespace std;

    	#include <fstream>
	#include<thread>
	#include <vector>	

	char ** g_data;
	char ** g_data_samp;
	int g_nsnp=0;
	int g_nbytes=0, g_nbytes_samp=0;
	int g_nsamp=0;
	int g_extra=8, g_extra_samp=8;
	int g_nthreads=1;
	int hethetTab[256][256];
	int homoppTab[256][256];
	vector<double> g_mu;
	vector<double> g_sd;

vector<int> getRange(int nObj, int nPart, int index)
{
	int chunkdiv=nObj/nPart;
	int chunkrem=nObj % nPart;
	vector<int> IndRange(2,0);

	if(index<(nPart-chunkrem))
	{
		IndRange[0]=index*chunkdiv;
		IndRange[1]=(index+1)*chunkdiv;
	} else {
		IndRange[0]=(nPart-chunkrem)*chunkdiv+(index-nPart+chunkrem)*(chunkdiv+1);
		IndRange[1]=(nPart-chunkrem)*chunkdiv+(index-nPart+chunkrem+1)*(chunkdiv+1);
	}
	
	return IndRange;
}


void transposedata(int beginInd, int endInd)
{
	for(int i=beginInd*4;i<min(endInd*4,g_nsnp);i++)
	{
		int pos=2*(i%4);
		unsigned int mask = 255-3*pow(2,pos);
		int nonNullInBlock=8;
		for(int j=0;j<g_nbytes;j++)		
		{
			char c = g_data[i][j];
			if(j==g_nbytes-1)	nonNullInBlock=g_extra;
    			for (int k = 0; k < nonNullInBlock; k++)	
   			{
        			char cur = ((c >> k) & 3);
				g_data_samp[j*4+k/2][i/4] = (g_data_samp[j*4+k/2][i/4] & mask) | (cur << pos);
        			k++;
    			}
		}
	}

}


// [[Rcpp::export]]
void readfile(string fname, int nsamp, vector<int> RandomLags, int num_threads)
{

	Rcout<<"Reading the file: "<<fname<<'\n';
	g_nsamp=nsamp;
	g_nsnp=RandomLags.size();
	if(nsamp % 4 == 0)
	{
		g_nbytes=nsamp/4;
	} else {
		g_nbytes=nsamp/4+1;
		g_extra=2*(nsamp%4);
	}

	g_data = new char*[g_nsnp];
    	for (int i = 0; i < g_nsnp; ++i)
        	g_data[i] = new char[g_nbytes];
	
	ifstream bedfile (fname.c_str(), ios::in | ios::binary );

	bedfile.ignore(3);                   // Ignore first two magic numbers and snp/sample major determining byte of plink

	for(int i=0;i<g_nsnp;i++)
	{
		bedfile.ignore(RandomLags[i]*g_nbytes);
		bedfile.read(g_data[i],g_nbytes);
	}
	Rcout<<"Reading BED file complete \n";
	Rcout<<"Number of SNPs read: "<<g_nsnp<<'\n';
	Rcout<<"Number of samples: "<<g_nsamp<<'\n';

	int num_threads_max = thread::hardware_concurrency();
	if(num_threads_max<1) num_threads_max=1;
	if((num_threads > num_threads_max) | (num_threads==0))	num_threads=num_threads_max;
	Rcout<<"Using "<<num_threads<<" threads\n";
	g_nthreads=num_threads;

	if(g_nsnp % 4 == 0)
	{
		g_nbytes_samp=g_nsnp/4;
	} else {
		g_nbytes_samp=g_nsnp/4+1;
		g_extra_samp=2*(g_nsnp%4);
	}
	g_data_samp = new char*[g_nsamp];
    	for (int i = 0; i < g_nsamp; ++i)
        	g_data_samp[i] = new char[g_nbytes_samp];
  	
	if(g_extra_samp>0)
	{
		for(int i=0;i<g_nsamp;i++)	g_data_samp[i][g_nbytes_samp-1] &= 0;
	}

	thread th[g_nthreads];

	for(int i=0;i<g_nthreads;i++)
	{
		vector<int> IndRange=getRange(g_nbytes_samp,g_nthreads,i);
		int beginInd=IndRange[0];
		int endInd=IndRange[1];

		th[i]=thread(transposedata,beginInd,endInd);
	}
	for(int i=0;i<g_nthreads;i++)
	{
		th[i].join();
	}

	Rcout<<"Processing samples complete \n";
	g_mu.reserve(g_nsnp);
	g_sd.reserve(g_nsnp);

}


void meanth(int beginInd, int endInd, vector<int> idmean, int noscale)
{
	int nids=idmean.size(),j=0,l=0;
	for(int i=beginInd;i<endInd;i++)
	{
		double mucur=0;
		int missmu=0;
		for(int idj=0;idj<nids;idj++)
		{
			j=idmean[idj]/4;
			l=2*(idmean[idj]%4);
			char cur = ((g_data[i][j] >> l) & 3);
			if(cur==2) mucur+=1;
			if(cur==1) missmu+=1;
        		if(cur==0) mucur+=2;
    		}
		mucur/=(2*(nids-missmu));
		g_mu[i]=mucur;
		if (noscale==0)
		{
			g_sd[i]=sqrt(2*mucur*(1-mucur));
		} else {
			g_sd[i]=1;
		}
	}


}

// [[Rcpp::export]]
void setmeans(vector<int> idmean,int noscale=0)
{
	thread th[g_nthreads];

	for(int i=0;i<g_nthreads;i++)
	{
		vector<int> IndRange=getRange(g_nsnp,g_nthreads,i);
		int beginInd=IndRange[0];
		int endInd=IndRange[1];

		th[i]=thread(meanth,beginInd,endInd,idmean,noscale);
	}
	for(int i=0;i<g_nthreads;i++)
	{
		th[i].join();
	}


	Rcout<<"Allele frequencies calculated \n";

}

// [[Rcpp::export]]
void deletedata()
{
	for(int i = 0; i < g_nsnp; ++i) 
	{
   		delete [] g_data[i];
	}
	delete [] g_data;
}

// [[Rcpp::export]]
void deletetranspose()
{
	for(int i = 0; i < g_nsamp; ++i) 
	{
   		delete [] g_data_samp[i];
	}
	delete [] g_data_samp;
}

void postmth(int beginInd, int endInd, const vector<double> &sumu, double u[], arma::mat &multout, int ncol, const vector<int> &idmult) 
{
	int nids=idmult.size(), j=0, l=0;	//j=byte position, l=within byte position
	for(int i=beginInd;i<endInd;i++)
	{
		vector<double> u1(ncol,0.0), u2(ncol,0.0), missu(ncol,0.0);
		int jold=0;
		char curdata = g_data[i][jold];
		for(int idj=0;idj<nids;idj++)
		{
			j=idmult[idj]/4;
			l=2*(idmult[idj]%4);
			if(!(j==jold))	curdata = g_data[i][j]; 
			jold=j;
			char cur = ((curdata >> l) & 3);
			if((int)cur==2)
			{
				for(int k=0;k<ncol;k++)		u1[k]+=u[k*nids+idj];
			}else if((int)cur==1)
			{
				for(int k=0;k<ncol;k++)		missu[k]+=u[k*nids+idj];
			}else if((int)cur==0)
			{
				for(int k=0;k<ncol;k++)		u2[k]+=u[k*nids+idj];
			}			

		}

		for(int k=0;k<ncol;k++)
		{
			multout(i,k)=(u1[k]+2*u2[k]+2*g_mu[i]*(missu[k]-sumu[k]))/g_sd[i];
		}
	}
} 


// [[Rcpp::export]]
arma::mat postmultiply(NumericMatrix u, vector<int> idmult)
{

	int ncol=u.ncol();

	thread th[g_nthreads];
	mat multout(g_nsnp,ncol);
	vector<double> sumu(ncol,0.0);
	int nids=idmult.size();
	for(int j=0;j<ncol;j++)
	{
		for(int k=0;k<nids;k++)
		{
			sumu[j]+=u(k,j);
		}
	}
	
	for(int i=0;i<g_nthreads;i++)
	{
		vector<int> IndRange=getRange(g_nsnp,g_nthreads,i);
		int beginInd=IndRange[0];
		int endInd=IndRange[1];

		th[i]=thread(postmth,beginInd,endInd,ref(sumu),u.begin(),ref(multout),ncol,ref(idmult));
	}
	for(int i=0;i<g_nthreads;i++)
	{
		th[i].join();
	}

	return multout;
}


void premth(int beginInd, int endInd, const vector<double> &ratio, const vector<double> &muratio, const vector<double> &muratiosum, arma::mat &multout, int ncol, const vector<int> &idmult) 
{
	for(int idi=beginInd;idi<endInd;idi++)
	{
		int i=idmult[idi];
		int nonNullInBlock=8;
		vector<double> u1(ncol,0.0), u2(ncol,0.0), missu(ncol,0.0);
		for(int j=0;j<g_nbytes_samp;j++)
		{
			char c = g_data_samp[i][j];
			if(j==g_nbytes_samp-1)	nonNullInBlock=g_extra_samp;
			if((c & 255)!=255)
			{
			
    			for (int l = 0; l < nonNullInBlock; l++)
   			{
				int snpj=j*4+l/2;
        			char cur = ((c >> l) & 3); 
				if(cur==2)
				{
					for(int k=0;k<ncol;k++)
						u1[k]+=ratio[snpj*ncol+k];
				}else if(cur==1)
				{
					for(int k=0;k<ncol;k++)
						missu[k]+=muratio[snpj*ncol+k];
				}else if(cur==0)
				{
					for(int k=0;k<ncol;k++)
						u2[k]+=ratio[snpj*ncol+k];
				}
        			l++;
    			}
			}
		}
		for(int k=0;k<ncol;k++)
		{
			multout(idi,k)=u1[k]+2*u2[k]-muratiosum[k]+missu[k];
		}
	}
} 



// [[Rcpp::export]]
arma::mat premultiply(NumericMatrix u, vector<int> idmult)
{
	int ncol=u.ncol();
	thread th[g_nthreads];

	int nids=idmult.size();

	mat multout(nids,ncol);

	vector<double> muratiosum(ncol,0.0), ratio(g_nsnp*ncol,0.0), muratio(g_nsnp*ncol,0.0);
	for(int k=0;k<g_nsnp;k++)
	{
		for(int j=0;j<ncol;j++)
		{
			ratio[k*ncol+j]=u(k,j)/g_sd[k];
			muratio[k*ncol+j]=ratio[k*ncol+j]*g_mu[k]*2;
			muratiosum[j]+=muratio[k*ncol+j];
		}	
	}

	for(int i=0;i<g_nthreads;i++)
	{
		vector<int> IndRange=getRange(nids,g_nthreads,i);
		int beginInd=IndRange[0];
		int endInd=IndRange[1];

		th[i]=thread(premth,beginInd,endInd,ref(ratio),ref(muratio),ref(muratiosum),ref(multout),ncol,ref(idmult));
	}
	for(int i=0;i<g_nthreads;i++)
	{
		th[i].join();
	}

	return multout;

}

void convertGTth(int thbeginInd, int thendInd, int beginBlockInd, arma::mat &GT, arma::mat &mu, const arma::mat &beta, const arma::mat &score)
{
	mu.rows(thbeginInd,thendInd-1)=beta.rows(thbeginInd,thendInd-1)*score.t();
	mu.rows(thbeginInd,thendInd-1)=clamp(mu.rows(thbeginInd,thendInd-1),0.0001,2-0.0001);

	int j=0,l=0; 	//j=byte position, l=within byte position
	for(int i=thbeginInd; i<thendInd;i++)
	{
		int idi=beginBlockInd+i;
		int jold=0;
		char curdata = g_data[idi][jold];
		for(int idj=0;idj<g_nsamp;idj++)
		{
			j=idj/4;
			l=2*(idj%4);
			if(!(j==jold))	curdata = g_data[idi][j]; 
			jold=j;
			char cur = ((curdata >> l) & 3);
			if((int)cur==2)
			{
				GT(i,idj)=1;
			}else if((int)cur==1)
			{
				GT(i,idj)=0;
				mu(i,idj)=0;
			}else if((int)cur==0)
			{
				GT(i,idj)=2;
			}			

		}
	}
	GT.rows(thbeginInd,thendInd-1)=GT.rows(thbeginInd,thendInd-1)-mu.rows(thbeginInd,thendInd-1);
	mu.rows(thbeginInd,thendInd-1).transform( [](double val) {return sqrt(2*val-pow(val,2)); } );
}

void calcRelth(int thbeginInd, int thendInd, const arma::mat &GT,const arma::mat &mu, arma::mat &kinout,const vector<int> &id1,const vector<int> &id2)
{
	for(int i=thbeginInd;i<thendInd;i++)
	{
		mat Num1=GT.col(id1[i]).t()*GT.col(id2[i]);
		mat Den1=mu.col(id1[i]).t()*mu.col(id2[i]);
		kinout(i,0)=Num1(0);
		kinout(i,1)=Den1(0);
	}
}


// [[Rcpp::export]]
arma::mat calcRel(int beginInd, int endInd, mat beta, mat score, vector<int> id1, vector<int> id2)
{
    auto start = std::chrono::system_clock::now();
	int nkin=id1.size(), nblock=endInd-beginInd;	//j=byte position, l=within byte position
	mat GT(nblock,g_nsamp);
	thread th[g_nthreads];
	vector<double> subg_mu(g_mu.begin()+beginInd,g_mu.begin()+endInd);
	vec meanmu(subg_mu);
	beta.col(0)=2*meanmu;
	mat mu(nblock,g_nsamp);

	for(int i=0;i<g_nthreads;i++)
	{
		vector<int> IndRange=getRange(nblock,g_nthreads,i);
		int thbeginInd=IndRange[0];
		int thendInd=IndRange[1];

		th[i]=thread(convertGTth,thbeginInd,thendInd,beginInd,ref(GT),ref(mu),ref(beta),ref(score));
	}
	for(int i=0;i<g_nthreads;i++)
	{
		th[i].join();
	}

    auto end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "Calculating ISAF using current block finished at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";

    start = std::chrono::system_clock::now();

	mat kinout(nkin,2);

	for(int i=0;i<g_nthreads;i++)
	{
		vector<int> IndRange=getRange(nkin,g_nthreads,i);
		int thbeginInd=IndRange[0];
		int thendInd=IndRange[1];

		th[i]=thread(calcRelth,thbeginInd,thendInd,ref(GT),ref(mu),ref(kinout),ref(id1),ref(id2));
	}
	for(int i=0;i<g_nthreads;i++)
	{
		th[i].join();
	}

    end = std::chrono::system_clock::now();
   elapsed_seconds = end-start;
    end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "Calculating kinships using current block finished at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";

	return kinout;
}

void createTable()
{
	for(unsigned int i=0;i<256;i++)
	{
		for(unsigned int j=0;j<256;j++)
		{
			int nhethet=0, nhomopp=0;
			for(int k=0;k<8;k++)
			{
				char curi=((i >> k) & 3);
				char curj=((j >> k) & 3);
				if(curi==2 && curj==2)
				{
					nhethet+=1;
				} else if((curi==3 && curj==0) || (curi==0 && curj==3)){
					nhomopp+=1;
				}
				k++;
			}
			hethetTab[i][j]=nhethet;
			homoppTab[i][j]=nhomopp;
		}
	}
}

void calcDiv(int beginInd, int endInd, double cutoff, const vector<int> &indexInkin0, const vector<int> &nhet, vector<int> &divout)
{
	int ntot=endInd-beginInd;
	for(int i=beginInd;i<endInd;i++)
	{
		int sampi=indexInkin0[i];
		for(int sampj=0;sampj<g_nsamp;sampj++)
		{
			if(sampj!=sampi)
			{
				int nhethet=0, nhomopp=0;
				for(int k=0;k<g_nbytes_samp;k++)
				{
					char curi = g_data_samp[sampi][k], curj=g_data_samp[sampj][k];
					nhethet+=hethetTab[(unsigned int) (curi & 255)][(unsigned int) (curj & 255)];
					nhomopp+=homoppTab[(unsigned int) (curi & 255)][(unsigned int) (curj & 255)];
				}
				double div=((double)(nhethet-2*nhomopp))/((double) (nhet[sampi]+nhet[sampj]));
				if(div<cutoff)	divout[i]++;
			}
		}
	}
}

// [[Rcpp::export]]
vector<int> calculateDivergence(vector<int> idkin0, double cutoff = -0.025)
{
	int nInkin0=idkin0.size();

	vector<int> nhet(g_nsamp,0);
	vector<int> divout(nInkin0,0);
	createTable();
	for(int i=0;i<g_nsamp;i++)
	{
		for(int k=0;k<g_nsnp;k++)
		{	
			char curi = (g_data_samp[i][k/4] >> 2*(k % 4)) & 3;
			if(curi==2)	nhet[i]++;
		}
	}

	thread th[g_nthreads];

	int numParts;
	if(nInkin0>100)
	{
		numParts=20;
	} else {
		numParts=1;
	}
	for(int bigi=0;bigi<numParts;bigi++)
	{
    		auto start = std::chrono::system_clock::now();
		vector<int> IndRange=getRange(nInkin0,numParts,bigi);
		int beginInd=IndRange[0];
		int endInd=IndRange[1];
		int setlen=endInd-beginInd;
		for(int i=0;i<g_nthreads;i++)
		{
			vector<int> IndRange=getRange(setlen,g_nthreads,i);
			int thbeginInd=IndRange[0];
			int thendInd=IndRange[1];

			th[i]=thread(calcDiv,thbeginInd+beginInd,thendInd+beginInd,cutoff,ref(idkin0),ref(nhet),ref(divout));
		}
		for(int i=0;i<g_nthreads;i++)
		{
			th[i].join();
		}
    		auto end = std::chrono::system_clock::now();
   		std::chrono::duration<double> elapsed_seconds = end-start;
    		std::time_t end_time = std::chrono::system_clock::to_time_t(end);

		Rcout<<(bigi+1)*5<<"% completed in "<<elapsed_seconds.count() << "s\n";
	}
	return divout;
}
