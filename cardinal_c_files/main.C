#include "mathfuncs.h"
#include <iostream>
int main()
{
    /*
    std::vector<Legendre*> leg_vec{};
    for (int i = 0; i <=10; i++)
    {
        Legendre* leg = new Legendre(i);
        leg->setExpandedCoeff(1.0);
        leg_vec.push_back(leg);
        std::cout<< leg->calcValue(.5)<<std::endl;
    }
    for (int i = 11; i < 100; i++)
    {
        Legendre* leg = new Legendre(i, leg_vec[i-1], leg_vec[i-2]);
        leg->setExpandedCoeff(1.0);
        leg_vec.push_back(leg);
        std::cout<< leg->calcValue(.5)<<std::endl;
    }
    */
    
    for (int n = 0; n < 7; n++)
    {
        std::vector<int> _m{n};
        for (int m = 0; m < n; m++)
        {
            _m.push_back(_m[m]-2);
        }
        for (int i = 0; i <= n; i++)
        {
            std::cout<<"Order"<<n<<" , "<<_m[i]<<std::endl;
            std::cout<<"    "<<calculateZernike(n,_m[i],.5,0.0,1.0)<<std::endl;
        }
        
    }
    
    return 0;
}