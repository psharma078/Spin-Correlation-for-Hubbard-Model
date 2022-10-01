#pragma once 
#include "itensor/util/print_macro.h"

using namespace itensor;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//This file computes superconducting pairing correlation
//<c^dag_{i,sigma} c^dag_{j,sigmap} c_{k,sigmap} c_{l,sigma}>
//with i < j < k < l for a given wafefunction.  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Real pairing_corr(MPS psi, const SiteSet sites, const string& op1, int i, const string& op2, int j, const string& op3, int k, const string& op4, int l){
    //auto N = length(psi);
    auto Adag_i = op(sites,op1,i);
    auto Adag_j = op(sites,op2,j);
    auto A_k = op(sites,op3,k);
    auto A_l = op(sites,op4,l);

    ITensor Fi, Fj, Fk, Fl;
    if (op1 == "Adagup") Fi = op(sites,"F",i);
    else Fi = op(sites,"Id",i);

    if (op2 == "Adagdn") Fj = op(sites,"F",j);
    else Fj = op(sites,"Id",j);

    if (op3 == "Aup") Fk = op(sites,"F",k);
    else Fk = op(sites,"Id",k);

    if (op4 == "Adn") Fl = op(sites,"F",l);
    else Fl = op(sites,"Id",l);
    
    psi.position(i);
    auto psidag = dag(psi);
    psidag.prime();

    auto li = leftLinkIndex(psi,i);
    auto rl = rightLinkIndex(psi,l);

    auto Pijkl = prime(psi(i),li)*Fi;
    Pijkl *= prime(Adag_i)*prime(psidag(i),"Site");

    for(int p = i+1; p < j; ++p)
    {
	Pijkl *= psi(p);
	Pijkl *= op(sites,"F",p); //Jordan-Wigner string
	Pijkl *= psidag(p);
    }

    Pijkl *= psi(j);
    Pijkl *= Adag_j*prime(Fj);
    Pijkl *= prime(psidag(j),"Site");
   
    for (int q = j+1; q < k; ++q)
    {
	Pijkl *= prime(psi(q),"Site");
	Pijkl *= psidag(q);
    }
    Pijkl *= psi(k);
    Pijkl *= Fk;
    Pijkl *= prime(A_k);
    Pijkl *= prime(psidag(k),"Site");
    
    for(int m = k+1; m < l; ++m)
    {
        Pijkl *= psi(m);
        Pijkl *= op(sites,"F",m); //Jordan-Wigner string
        Pijkl *= psidag(m);
    }

    Pijkl *= prime(psi(l),rl);
    Pijkl *= A_l * prime(Fl);
    Pijkl *= prime(psidag(l),"Site");

   // print(i," ",j," ",k," ",l," ", elt(Pijkl),"\n");

  return elt(Pijkl);

}
