#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include "bondPairs.h"
#include "spinCorrelation.h"
#include "S2.h"
#include <iostream>
#include <stdlib.h>
#include <vector>

using namespace itensor;

//**** finding entry in an array ****
template<class C, typename T>
bool contains(C&& c, T t) {
    return std::find(std::begin(c), std::end(c), t) != std::end(c);
};

//**** random number generator ****
template <typename T>
std::vector<int> random_sample(int n, int k, T gen)
{
    std::vector<int> out;
    std::vector<int> site_indices;
    for (int i = 1 ; i <= n ; ++i) {
        site_indices.push_back(i);
    }
    std::sample(site_indices.begin(), site_indices.end(), std::back_inserter(out), k, gen);
    return out;
}

//selects odd random number
std::vector<int> selectRandomOddNumbers(int n, int N) {
    std::vector<int> oddNumbers(N / 2);
    std::iota(oddNumbers.begin(), oddNumbers.end(), 1);
    std::shuffle(oddNumbers.begin(), oddNumbers.end(), std::mt19937(std::random_device()()));

    oddNumbers.resize(std::min(n, static_cast<int>(oddNumbers.size())));
    std::vector<int> oddnum;
    for (int i : range(n)){
       oddnum.push_back(oddNumbers[i]*2-1);
    }

    return oddnum;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int
main(int argc, char* argv[])
{
 auto Nx = 12,
      Ny = 6;
 auto N = Nx*Ny;

 auto Nup = atof(argv[1]);
 auto Ndn = atof(argv[2]);
 
 printf("Nx = %d, Ny = %d\n", Nx, Ny);
 printf("Nup = %d, Ndn = %d\n", Nup, Ndn);

 auto t = 1;
 auto U = 20;	
 auto h = 0.2;
 printf("t = %d\n", t);
 printf("U = %d\n", U);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 auto sites = Electron(N);
 //
 //Make a Hamiltonian on the triangular lattice using AutoMPO
 //
 auto ampo = AutoMPO(sites);
 auto bonds = bondPairs(Nx, Ny); 
 printfln("Triangular XC6 lattice NN pairs");
 
 for(auto& bnd : bonds)
    {
    printfln("Bond from site %d -> %d",bnd.first,bnd.second);
    }

 for(auto& bnd : bonds)
    {
    ampo += -t,"Cdagup",bnd.first,"Cup",bnd.second;
    ampo += -t,"Cdagup",bnd.second,"Cup",bnd.first;
    ampo += -t,"Cdagdn",bnd.first,"Cdn",bnd.second;
    ampo += -t,"Cdagdn",bnd.second,"Cdn",bnd.first;
    }
 for(auto j : range1(N))
    {
    ampo += U,"Nupdn",j;
    }

 auto H = toMPO(ampo);

 for(auto j : range1(N))
    {
    ampo += -h,"S2",j;
    }
 
 auto HS2 = toMPO(ampo); 
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //initializing wavefunction MPS 

 auto state = InitState(sites);

 std::mt19937 gen(1);
 int doubly_occ = int(Nup+Ndn-N);
 std::vector<int> r = random_sample(N, doubly_occ, gen);

 std::vector<int> ar;
 for (auto j : range1(N)){
        bool exists = contains(r, j);
        if (exists) {state.set(j, "UpDn");}
        else {
        state.set(j,"Emp");
        ar.push_back(j);
        }
 }

 random_shuffle(std::begin(ar), std::end(ar));
 //random_shuffle(std::begin(ar), std::end(ar));
 
 int nupAdditional = Nup - doubly_occ;
 //int ndnAdditional = Ndn - doubly_occ;
 for (int i : range(ar.size())){
        if (i < nupAdditional)
                state.set(ar[i], "Up");
        else
                state.set(ar[i], "Dn");
 }
 auto psi0 = MPS(state);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //setting sweep parameters
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    println("Warmup sweeps\n");
    auto sweeps = Sweeps(11);
    sweeps.maxdim() = 8,16,24,50,100,250,500,1000,1000,2000,2000;
    sweeps.cutoff() = 1E-8;
    sweeps.niter() = 4;//change to 4
    sweeps.noise() = 1E-7,1E-8,0.0;
    println(sweeps);

    //warming up initial state
 auto [Energy, psi1] = dmrg(HS2,psi0,sweeps,"Quiet");    
 println("Done warming up\n");
    //
    //final sweeps
    //
    auto fsweeps = Sweeps(20);
    fsweeps.maxdim() = 2000,3000,3000,3000,5000,5000,5000,7000,7000,7000,7000,10000,10000,10000,10000,10000,14000;
    fsweeps.cutoff() = 1E-8;
    fsweeps.niter() = 4;//change to 4
    fsweeps.noise() = 1E-7,1E-8,0.0;
    println(fsweeps);
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //DMRG calculation for ground state
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 auto [energy,psi] = dmrg(H,psi1,fsweeps,{"WriteM",17000}); //"Quiet");

 printfln("Initial energy = %.5f", inner(psi0,H,psi0) );
 printfln("\nGround State Energy from DMRG= %.10f",energy);
 printfln("\nEnergy using overlap = %.10f", inner(psi,H,psi));
 printfln("\nEnergy per site = %.10f",energy/N);

 auto E2 = inner(H,psi,H,psi);
 printfln("\n<psi|H^2|psi> = %.10f",E2);
 printfln("\nvariance = %.10f", E2-energy*energy);

 printfln("\noverlap with initial state= %.10f", inner(psi0,psi) );

 println("\nTotal QN of Ground State = ",totalQN(psi));

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //Saving MPS to a file
 writeToFile("NoField_sites_file_nup="+str(Nup)+"_ndn="+str(Ndn),sites);
 writeToFile("psi_file_nup="+str(Nup)+"_ndn="+str(Ndn),psi);
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //claculating <Nup_i>, <Ndn_i>, <Nupdn_i>
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 /*std::cout<<"i <Nup_i> <Ndn_i> <Nupdn_i>"<<std::endl;
 for (int i : range1(N)){
        auto nn = sites.op("Nupdn",i);
	auto nup =sites.op("Nup",i);
	auto ndn =sites.op("Ndn",i);
        psi.position(i);
        auto Cupdn = psi(i)*nn;
	auto Cup = psi(i)*nup;
	auto Cdn = psi(i)*ndn;
	Cup *= dag(prime(psi(i),"Site"));
	Cdn *= dag(prime(psi(i),"Site"));
        Cupdn *= dag(prime(psi(i),"Site"));
	std::cout<< i <<" "<< elt(Cup)<<" "<<elt(Cdn)<<" "<<elt(Cupdn)<<std::endl;
        }
std::cout<<std::endl;*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//calculating <S_i^z>
printfln("magnetization");
for (int i : range1(N)){
        auto sz = sites.op("Sz",i);
        psi.position(i);
        //auto C = dag(prime(psi(i),"Site"))*sites.op("Sz",i)*psi(i);
        auto C = psi(i)*sz;
	C *= dag(prime(psi(i),"Site"));
        printfln("site ",i," <Sz>= ",elt(C));
}
std::cout<<std::endl;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//computing spin correlation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
spinCorrelator(psi, sites);

printfln("Initial state");
printfln("i"," nup", " ndn"," nupdn");
for (int i : range1(N)){
        auto nn = sites.op("Nupdn",i);
        auto nup =sites.op("Nup",i);
        auto ndn =sites.op("Ndn",i);
        psi0.position(i);
        auto Cupdn = psi0(i)*nn;
        auto Cup = psi0(i)*nup;
        auto Cdn = psi0(i)*ndn;
        Cup *= dag(prime(psi0(i),"Site"));
        Cdn *= dag(prime(psi0(i),"Site"));
        Cupdn *= dag(prime(psi0(i),"Site"));
        std::cout<< i <<" "<< elt(Cup)<<" "<<elt(Cdn)<<" "<<elt(Cupdn)<<std::endl;
        }

//auto s2 = sites.op("S2",2);
//PrintData(s2);

return 0;
}

