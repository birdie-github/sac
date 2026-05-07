#ifndef _DOMAIN_H
#define _DOMAIN_H

#include "../global.h"
#include <cmath>

class LogDomain {
  public:
    static constexpr int scale=256;
    static constexpr int dbits=12;
    static constexpr int dmin=-2047;
    static constexpr int dmax=2047;
    static constexpr int dscale=dmax-dmin+1;
    int fmin,fmax;
    LogDomain()
    {
      for (int i=0;i<PSCALE;i++)
      {
        double f=std::max(i,1)/(double)PSCALE;
        double q=std::log(f / (1.0-f))*scale;
        FwdTbl[i]=(int)std::round(q);
      };
      fmin=FwdTbl[0];
      fmax=FwdTbl[PSCALE-1];
      // 12-Bit
      for (int i=dmin;i<=dmax;i++)
      {
        double q=PSCALE/(1.0+exp(-double(i)/double(scale)));
        InvTbl[i-dmin]=(int)std::round(q);
      };
    }
    inline int Fwd(int p)
    {
       return FwdTbl[p];
    }
    inline int Inv(int x)
    {
       if (x<dmin) return 1;
       else if (x>dmax) return PSCALEm;
       else return InvTbl[x-dmin];
    }
    void Check()
    {
      int sum=0;
      for (int i=0;i<PSCALE;i++)
      {
        int p=Inv(Fwd(i));
        sum+=(p-i)*(p-i);
      }
      printf(" mse: %0.4f\n",double(sum)/double(PSCALE));
    }
  protected:
    int FwdTbl[PSCALE];
    int InvTbl[dscale];
};

inline LogDomain myDomain;

#endif
