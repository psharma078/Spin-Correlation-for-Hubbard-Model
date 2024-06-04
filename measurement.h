#pragma once
#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include <vector>
using namespace itensor;

template <typename T1, typename T2, typename T3>
void collectdata(const T1& x, const T2& y, const T3& corr, const std::string& fname)
{
    std::ofstream outfile(fname, std::ios::app);
    outfile << x << "," << y << "," << corr << "\n";
    outfile.close();
}

template <typename T1, typename T2>
void collectdata(const T1& x, const T2& corr, const std::string& fname)
{
    std::ofstream outfile(fname, std::ios::app);
    outfile << x << "," << corr << "\n";
    outfile.close();
}

std::vector<std::vector<Real>> DimerCorr(MPS& psi, const SiteSet sites)
{
    auto N = length(psi);
    auto spsm = correlationMatrix(psi,sites,"S+","S-");
    auto smsp = correlationMatrix(psi,sites,"S-","S+");
    auto szsz = correlationMatrix(psi,sites,"Sz","Sz");
    std::vector<std::vector<Real>> Corr(N, std::vector<Real>(N, 0.0));
    for (int i : range(N))
    {
        for (int j : range(N))
        {
            Corr[i][j] = 0.5*(spsm[i][j]+smsp[i][j])+szsz[i][j];
        }
    }
    return Corr;
}

void holeDimerCorr(MPS psi, const SiteSet sites, std::string fname)
{
    auto N = length(psi);
    //projecting MPS with hole fixed at N/2
    auto so = N/2;
    auto hole_op = sites.op("Cup*Cdagup*Cdn*Cdagdn",so);
    auto A = hole_op*psi(so);
    A.noPrime();
    psi.set(so,A);
    psi.position(so);
    psi /= norm(psi);
    auto corr = DimerCorr(psi,sites);
    auto Sz = expectC(psi,sites,"Sz");
    for (int x : range(N))
    {
         collectdata(x+1,Sz[x].real(),"magnetization_Sz_N=100_tp=0_1H49M.csv");
    }
    for (int i : range(N))
    {
        for (int j : range(N))
        {
            if (i>j)
            {
            collectdata(i+1,j+1,corr[i][j],fname);
            }
        }
    }
}

void holeSingletCorr(MPS psi, const SiteSet sites, std::string fname)
{
    auto N = length(psi);
    auto hole_op = sites.op("Cup*Cdagup*Cdn*Cdagdn",N/2);
    auto A = hole_op*psi(N/2);
    A.noPrime();
    psi.set(N/2,A);
    psi.position(N/2);
    psi /= norm(psi);
    auto NupNdn = correlationMatrix(psi,sites,"Nup","Ndn");
    for (int i : range(N))
    {
        for (int j : range(N))
        {
            auto corr = NupNdn[i][j]+NupNdn[j][i];
            collectdata(i+1,j+1,corr,fname);
        }
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//hole hoping sum_s (Cdagup_i Cdag_i+1 + Cdagdn_i Cdn_i+1)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MPS twoHop(MPS psi, const SiteSet sites,int sh)
{
    auto hop = -1*sites.op("Adagup*F",sh)*sites.op("F",sh+1)*sites.op("Aup",sh+2);
    hop += -1*sites.op("Adagdn",sh)*sites.op("F",sh+1)*sites.op("F*Adn",sh+2);
    auto B = psi(sh)*psi(sh+1)*psi(sh+2)*hop;
    B.noPrime();
    auto l1 = leftLinkIndex(psi,sh);
    auto [U1,S1,V1] = svd(B,{l1,sites(sh)});
    auto C = S1*V1;
    auto l2 = commonIndex(U1,C);
    auto [U2,S2,V2] = svd(C,{l2,sites(sh+1)});
    psi.set(sh,U1);
    psi.set(sh+1,U2);
    psi.set(sh+2,S2*V2);
    psi.position(sh);
    psi /= norm(psi);
    return psi;
}

MPS oneHop(MPS psi, const SiteSet sites,int sh)
{
    auto hop = -1*sites.op("Adagup*F",sh)*sites.op("Aup",sh+1);
    hop += -1*sites.op("Adagdn",sh)*sites.op("F*Adn",sh+1);
    auto B = psi(sh)*psi(sh+1)*hop;
    B.noPrime();
    auto l1 = leftLinkIndex(psi,sh);
    auto [U,S,V] = svd(B,{l1,sites(sh)});
    psi.set(sh,U);
    psi.set(sh+1,S*V);
    psi.position(sh);
    psi /= norm(psi);
    return psi;
}

void holeDynamics(MPS psi, const SiteSet sites, std::string fname)
{
    auto N = length(psi);
    //projecting MPS with hole fixed at N/2
    auto hole_op = sites.op("Cup*Cdagup*Cdn*Cdagdn",N/2);
    auto A = hole_op*psi(N/2);
    A.noPrime();
    psi.set(N/2,A);
    psi.position(N/2);
    psi /= norm(psi);

    auto phi = psi;

    psi = oneHop(psi, sites, N/2);//first hop
    psi = oneHop(psi,sites,N/2+1);//second hop
    phi = twoHop(phi, sites, N/2);//direct hop
    psi.position(N/2);

    auto wf = sum(psi,phi,{"MaxDim", 500, "Cutoff", 1E-8});//while adding MPS, orthogonality center must match, other return an error "segmentation fault".
    wf.position(N/2);
    wf /= norm(wf);
    auto corr1 = DimerCorr(psi,sites);
    auto corr2 = DimerCorr(phi,sites);
    auto corr3 = DimerCorr(wf,sites);
    std::string fname1 = "holeDyn_twoOneHop.csv";
    std::string fname2 = "holeDyn_twoHop.csv";
    std::string fname3 = "holeDyn_finalCombo.csv";
    for (int i : range(N))
    {
        for (int j : range(N))
        {
            if (i>j)
            {
            collectdata(i+1,j+1,corr1[i][j],fname1);
            collectdata(i+1,j+1,corr2[i][j],fname2);
            collectdata(i+1,j+1,corr3[i][j],fname3);
            }
        }
    }
}
~           
