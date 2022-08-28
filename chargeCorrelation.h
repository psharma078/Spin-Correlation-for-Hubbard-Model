#pragma once 
#include "itensor/util/print_macro.h"

using namespace itensor;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//computing charge correlation <ntot_i ntot_j>
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void chargeCorrelator(MPS psi, const SiteSet sites){
    auto N = length(psi);
    printfln("m n <Ntot_m Ntot_n>");

    for (int i=1;i<=N;i++){
	psi.position(i);
	auto ni = sites.op("Ntot",i);

     	for (int j=i;j<=N;j++){
	    auto C = psi(i)*ni;
	    if (j==i) {
		C *=  prime(ni);
 	    	C *= dag(prime(prime(psi(i),"Site"),"Site"));
	    }
	    else {
	    	auto ir = commonIndex(psi(i),psi(i+1),"Link");
	    	C *= dag(prime(prime(psi(i),"Site"),ir));

	    	for (int k=i+1;k<j;k++){
		    C *= psi(k);
		    C *= dag(prime(psi(k),"Link"));
	    	}

		C *= psi(j)*sites.op("Ntot",j);
		auto il = commonIndex(psi(j),psi(j-1),"Link");
		C *= dag(prime(prime(psi(j),"Site"),il));
	    }
	printfln("",i," ",j," ",elt(C));
     }
 }

}

