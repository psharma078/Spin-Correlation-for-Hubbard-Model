#pragma once 
#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//computing spin correlation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void spinCorrelator(MPS psi, const SiteSet sites){
    //auto sites = siteInds(psi);
    auto N = length(psi);
    printfln("i j <Si^z Sj^z> <Si^+ Sj^-> <Si^- Sj^+>");

    for (int i=1;i<=N;i++){
	psi.position(i);
	auto szi = sites.op("Sz",i);
     	auto spi = sites.op("S+",i);
     	auto smi = sites.op("S-",i);

     	for (int j=i;j<=N;j++){
	    auto C = psi(i)*szi;
	    auto Cpm = psi(i)*spi;
	    auto Cmp = psi(i)*smi;
	    if (j==i) {
		C *=  prime(szi);
 	    	C *= dag(prime(prime(psi(i),"Site"),"Site"));
	    	//auto Cc = psi(i)*sites.op("Sz2",i)*dag(prime(psi(i),"Site"));
	   	Cpm *= prime(sites.op("S-",i));
	    	Cpm *= dag(prime(prime(psi(i),"Site"),"Site"));
	    	Cmp *= prime(sites.op("S+",i));
	    	Cmp *= dag(prime(prime(psi(i),"Site"),"Site"));
	    }
	    else {
	    	auto ir = commonIndex(psi(i),psi(i+1),"Link");
	    	C *= dag(prime(prime(psi(i),"Site"),ir));
	    	Cpm *= dag(prime(prime(psi(i),"Site"),ir));
	    	Cmp *= dag(prime(prime(psi(i),"Site"),ir));

	    	for (int k=i+1;k<j;k++){
		    C *= psi(k);
		    C *= dag(prime(psi(k),"Link"));
		    Cpm *= psi(k);
		    Cpm *= dag(prime(psi(k),"Link"));
		    Cmp *= psi(k);
		    Cmp *= dag(prime(psi(k),"Link"));
	    	}

		C *= psi(j)*sites.op("Sz",j);
		Cpm *= psi(j)*sites.op("S-",j);
		Cmp *= psi(j)*sites.op("S+",j);
		auto il = commonIndex(psi(j),psi(j-1),"Link");
		C *= dag(prime(prime(psi(j),"Site"),il));
		Cpm *= dag(prime(prime(psi(j),"Site"),il));
		Cmp *= dag(prime(prime(psi(j),"Site"),il));
	    }
	printfln("",i," ",j," ",elt(C)," ",elt(Cpm)," ",elt(Cmp));
     }
 }

}

