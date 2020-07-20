
namespace weld {

   namespace bezier_functionality {

      class Bernstein {
         private: 
            int n;
            mutable std::vector<double> scratch;
         public:
            Bernstein(int _n) : n(_n), scratch(_n+1) {}
            double eval(int i, double u) const {
               double val=-1.0;
               for(int j=0;j<=n;j++)
                  scratch[j]=0.0;
               scratch[n-i]=1.0;
               double u1=1-u;
               for(int k=1;k<=n;k++)
                  for(int j=n;j>=k;j--)
                     temp[j]=u1*scratch[j]+u*scratch[j-1];
            
               val=scratch[n];
               return val;
            }
      }
   }
}

