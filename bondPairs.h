#pragma once

#include<vector>

std::vector<std::pair<int, int> > bondPairs(const int Nx, const int Ny){
    int N = Nx*Ny;
    std::pair<int,int> p1;
    std::vector<std::pair<int,int> > bonds;
 for (int i=1;i<N;i++)
    {
        if (i%Ny==1 and i<N-Ny){
                p1 = std::make_pair(i,i+1);
                bonds.push_back(p1);
                p1 = std::make_pair(i,i+Ny);
                bonds.push_back(p1);
                p1 = std::make_pair(i,i+Ny+1);
                bonds.push_back(p1);
                p1 = std::make_pair(i,i+Ny-1);
                bonds.push_back(p1);
                p1 = std::make_pair(i,i+2*Ny-1);
                bonds.push_back(p1);
            }
        else if (i%Ny==0 and i<N){
                p1 = std::make_pair(i,i+Ny);
                bonds.push_back(p1);
            }

	    else if (i>N-Ny){
                p1 = std::make_pair(i,i+1);
                bonds.push_back(p1);
                if (i==N-Ny+1){
                        p1 = std::make_pair(i,i+Ny-1);
                        bonds.push_back(p1);
                    }
            }

        else if (i%2==1){
                p1 = std::make_pair(i,i+1);
                bonds.push_back(p1);
                p1 = std::make_pair(i,i+Ny-1);
                bonds.push_back(p1);
                p1 = std::make_pair(i,i+Ny);
                bonds.push_back(p1);
                p1 = std::make_pair(i,i+Ny+1);
                bonds.push_back(p1);
            }

        else if (i%2==0){
                p1 = std::make_pair(i,i+1);
                bonds.push_back(p1);
                p1 = std::make_pair(i,i+Ny);
                bonds.push_back(p1);
            }
    }

    return bonds;
}
