//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: UBC ADMB May Workshop													 
//Date:	May 5-9, 2013														 
//Purpose:pm for hake data set											 
//Notes: 				 
//							 
//																 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	init_int nyrs; //read in the max age from the data file
	init_vector yt(1,nyrs); // vector for ages
	init_vector ct(1,nyrs); //vector for lengths
	init_int eof;
	int iter;
	!!iter=0;
	LOCAL_CALCS
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS

PARAMETER_SECTION
	init_bounded_number r(0.001,1);//use bounded number to prevent  instability
	init_number log_k;
	init_number log_q;
	init_number log_sig;
	number fpen;
	!!r=0.6;// set initial values can also use the .pin file for this. Need to add this later
	!!log_k=log(3.5*max(ct));
	!!log_q=log(0.001);
	!!log_sig=log(0.3);
	objective_function_value nll;//must define the objective function
	
	sdreport_number k;
	sdreport_number q;
	sdreport_number sig;
	vector Bt(1,nyrs+1);// we are not searching on this vector but do need to keep track of the derivatives of these values
	vector it(1,nyrs);// we are not searching on this vector but do need to keep track of the derivatives of these values
	vector epsilon(1,nyrs);// we are not searching on this vector but do need to keep track of the derivatives of these values

PROCEDURE_SECTION
	statedynamics();
	observationmodel();
	objectivefunction();

	if(mceval_phase()) mcmc_output();

FUNCTION statedynamics
	k=mfexp(log_k);//transfor the log parameter that are being searched on to regular space.
	Bt(1) = k;
	fpen = 0.0;
	for(int i=1; i<=nyrs; i++)
	{
		Bt(i+1)=posfun(Bt(i)+r*Bt(i)*(1.-Bt(i)/k)-ct(i),0.001,fpen);
	}

FUNCTION observationmodel
	q=mfexp(log_q);
	it=q*Bt(1,nyrs);

FUNCTION objectivefunction 
	//calculate objective function
	sig=mfexp(log_sig);
	epsilon=log(it)-log(yt);//calculate deviations from the mean
	//cout<<fpen<<endl;
	nll=dnorm(epsilon,sig);//use the statslib dnorm function to caalculate the total negative log likelihood.
	nll+=1.e5*fpen;

FUNCTION mcmc_output
	
	if(iter==0)
	{
		ofstream ofs("refpar.mcmc");
		ofs<<"fmsy\t bmsy\t msy\t b/bmsy\t f/fmsy"<<endl;
		ofstream ofs1("Bt.mcmc");
		//ofs<<"r\t k\t q\t sig\t"<<endl;
	}
	iter++;
	double fratio=value(-log(1-ct[nyrs]/Bt[nyrs])/(r/2));
	double bratio=value(Bt[nyrs]/(k/2));
	ofstream ofs("refpar.mcmc",ios::app);
	ofs<<r/2<<"\t"<<k/2<<"\t"<<r*k/4<<"\t"<<bratio<<"\t"<<fratio<<endl;
	ofstream ofs1("Bt.mcmc",ios::app);
	ofs1<<Bt<<endl;
	

REPORT_SECTION
	REPORT(k);// write to the report file using the REPORT function
	REPORT(r);// write to the report file using the REPORT function
	REPORT(q);// write to the report file using the REPORT function
	REPORT(sig);// write to the report file using the REPORT function
	REPORT(ct);// write to the report file using the REPORT function
	REPORT(yt);// write to the report file using the REPORT function
	REPORT(Bt);// write to the report file using the REPORT function
	REPORT(it);// write to the report file using the REPORT function

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
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	#include<contrib.h>//IF you have ADMB-11
	//#include<stats.cxx>//If you have ADMB-10 and make sure stats.cxx is in your working directory
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;


