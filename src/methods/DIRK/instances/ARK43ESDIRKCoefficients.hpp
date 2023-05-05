#ifndef ARK43ESDIRKCOEFFICIENTS_DEFINED__
#define ARK43ESDIRKCOEFFICIENTS_DEFINED__

#include "SingleRateMethodCoefficients.hpp"

using namespace arma;

class ARK43ESDIRKCoefficients: public SingleRateMethodCoefficients {
public:
	ARK43ESDIRKCoefficients() {
		name = "ARK4(3)6L[2]SA-ESDIRK";
		
		num_stages = 6;
		primary_order = 4.0;
		secondary_order = 3.0;
		
		A = mat(num_stages,num_stages,fill::zeros);
      	A(1,0) = 1.0/4.0; 
      	A(1,1) = 1.0/4.0;
      		
      	A(2,0) = 8611.0/62500.0; 
      	A(2,1) = -1743.0/31250.0; 
      	A(2,2) = 1.0/4.0;
      	
      	A(3,0) = 5012029.0/34652500.0;
      	A(3,1) = -654441.0/2922500.0; 
      	A(3,2) = 174375.0/388108.0;
      	A(3,3) = 1.0/4.0;
      	
      	A(4,0) = 15267082809.0/155376265600.0; 
      	A(4,1) = -71443401.0/120774400.0;
      	A(4,2) = 730878875.0/902184768.0;
      	A(4,3) = 2285395.0/8070912.0; 
      	A(4,4) = 1.0/4.0;
      
      	A(5,0) = 82889.0/524892.0;
      	A(5,1) = 0.0;
      	A(5,2) = 15625.0/83664.0;
      	A(5,3) = 69875.0/102672.0; 
      	A(5,4) = -2260.0/8211.0; 
      	A(5,5) = 1.0/4.0;

		b = vec(num_stages,fill::zeros);
		b(0) = 82889.0/524892.0;
      	b(1) = 0.0;
      	b(2) = 15625.0/83664.0;
      	b(3) = 69875.0/102672.0; 
      	b(4) = -2260.0/8211.0; 
      	b(5) = 1.0/4.0;

		d = vec(num_stages,fill::zeros);
		d(0) = 4586570599.0/29645900160.0; 
		d(1) = 0.0;
		d(2) = 178811875.0/945068544.0;
      	d(3) = 814220225.0/1159782912.0; 
      	d(4) = -3700637.0/11593932.0; 
      	d(5) = 61727.0/225920.0;

		c = A*vec(num_stages,fill::ones);
	}
};

#endif