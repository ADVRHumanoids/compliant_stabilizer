#include <compliant_stabilizer/FilterHClass.h>

FilterH::FilterH(){
//Default constructure. initialize all filter variables in 0

	Ncoeff=6; 					// Number of coefficients
	int i=0;
	for (i=0;i<MAX_FILTER_LENGTH;i++){ //Init filter in 0
		x[i]=0;
		y[i]=0;
		a[i]=0;
		b[i]=0;
	}
};

void FilterH::clear_filter(){
//clear the filter as constructure does
	int i=0;
	for (i=0;i<MAX_FILTER_LENGTH;i++){ //init filter in 0
		x[i]=0;
		y[i]=0;
		a[i]=0;
		b[i]=0;
	}
};

void FilterH::initialValue(double initialState){
//clear the filter as constructure does
    int i=0;
    for (i=0;i<MAX_FILTER_LENGTH;i++){ //init filter in 0
        x[i]=initialState;
        y[i]=initialState;
    }
};
void FilterH::least_squares_filter(ftype T, int N){
//recibe T as sample Time (s)
//N as the order for the filter

	ftype freq=1/T;  //Define the frequency
	clear_filter();  //prepared the filter
	
	//accoriding to the order N do
	if (N == 4) {
		Ncoeff = 4;
		b[0] = -0.3*freq;
		b[1] = -0.1*freq;
		b[2] =  0.1*freq;
		b[3] =  0.3*freq;
	}
	else if (N == 8) {
		Ncoeff = 8;
		b[0] = -0.0833*freq;
		b[1] = -0.0595*freq;
		b[2] = -0.0357*freq;
		b[3] = -0.0119*freq;
		b[4] =  0.0119*freq;
		b[5] =  0.0357*freq;
		b[6] =  0.0595*freq;
		b[7] =  0.0833*freq;
	}
	else { //Fail gracefully 
		Ncoeff = 2;
		b[0] = -freq;
		b[1] =  freq;
	}
};

void FilterH::moving_average_filter(int N){
// recive the order of the filter through the variable N
	clear_filter();
	int cnt;
	ftype C; //constant divisor
	if (N > MAX_FILTER_LENGTH)
		N = MAX_FILTER_LENGTH;
	
	clear_filter();
	Ncoeff = N;
	C = 1.0/N;
	for (cnt = 0; cnt < N; cnt ++) {
		a[cnt] = C;
	}
};


/** Build a butterworth filter. T is the sample period, 
 *  cutoff is the cutoff frequency in hertz, N is the order (1,2,3 or 4)
 */
void FilterH::butterworth(ftype T, ftype cutoff, int N){
	//recive sample time (s) T , cute off frequency (Hz) cutoff, and the order of the filter 
	//Generate the filter coefficients where a[0] is the cero order factor in the numerator
	clear_filter();
	ftype A;
	if (N>4)
		N=4;
	if (N==0)
		N=1;
	ftype C=1/tan(M_PI*cutoff*T);
	if (N==1){
		A=1/(1+C);
		a[0]=A;
		a[1]=A;
		b[0]=1;
		b[1]=(1-C)*A;
	}
	if (N==2){
		A=1/(1+1.4142135623730950488016887242097*C+pow(C,2));
		a[0]=A;
		a[1]=2*A;
		a[2]=A;

		b[0]=1;
		b[1]=(2-2*pow(C,2))*A;
		b[2]=(1-1.4142135623730950488016887242097*C+pow(C,2))*A;
	}
	if (N==3){
		A=1/(1+2*C+2*pow(C,2)+pow(C,3));
		a[0]=A;
		a[1]=3*A;
		a[2]=3*A;
		a[3]=A;

		b[0]=1;
		b[1]=(3+2*C-2*pow(C,2)-3*pow(C,3))*A;
		b[2]=(3-2*C-2*pow(C,2)+3*pow(C,3))*A;
		b[3]=(1-2*C+2*pow(C,2)-pow(C,3))*A;
	}
	if (N==4){
		A=1/(1+2.6131259*C+3.4142136*pow(C,2)+2.6131259*pow(C,3)+pow(C,4));
		a[0]=A;
		a[1]=4*A;
		a[2]=6*A;
		a[3]=4*A;
		a[4]=A;
		
		b[0]=1;
		b[1]=(4+2*2.6131259*C-2*2.6131259*pow(C,3)-4*pow(C,4))*A;
		b[2]=(6*pow(C,4)-2*3.4142136*pow(C,2)+6)*A;
		b[3]=(4-2*2.6131259*C+2*2.6131259*pow(C,3)-4*pow(C,4))*A;
		b[4]=(1-2.6131259*C+3.4142136*pow(C,2)-2.6131259*pow(C,3)+pow(C,4))*A;
	}
};

double FilterH::applyFilter(double X){
	//assumes the filter was already run so that the coefficients are able
	//recive the new input- retunrn the filter signal
	x[4]=x[3];
	x[3]=x[2];
	x[2]=x[1];
	x[1]=x[0];
	x[0]=X;
	y[4]=y[3];
	y[3]=y[2];
	y[2]=y[1];
	y[1]=y[0];
	y[0]=a[0]*x[0]+x[1]*a[1]+x[2]*a[2]+x[3]*a[3]+x[4]*a[4]-y[1]*b[1]-y[2]*b[2]-y[3]*b[3]-y[4]*b[4];
	return y[0];
};
/** Build a butterworth differentiator. T is the sample period, 
 *  cutoff is the cutoff frequency in hertz, N is the order (1,2 or 3)
 */


void FilterH::butterDifferentiator(ftype T, ftype cutoff, int N){
	//recive sample time T, cute off frequency cutoff, and the order of the filter
	//Generate the filter coefficients where a[0] is the cero order factor in the numerator
	clear_filter();
	ftype C =1/tan(M_PI*cutoff*T);
	ftype w=2*C/T;
	int i;
	for (i=0;i<MAX_FILTER_LENGTH;i++) //init filter in 0
		x[i]=1;
	
	if (N==1){
		
		a[0]= 1; 
		a[1]= -1;
		b[0]=T*(1+C)/2;
		b[1]=T*(1-C)/2;
	}
	if(N==2){
		
		a[0]= 1; 
		a[1]= 0; 
		a[2]= -1;
		b[0]= T/2*(1+1.414213562373095*C+pow(C,2));
		b[1]=T/2*(2-2*pow(C,2));
		b[2]=T/2*(1-1.414213562373095*C+pow(C,2));
	}
	if (N==3){
			
		a[0]=1;   //2*(T/2) multiplicative factor
		a[1]=1;   
		a[2]=-1;
		a[3]=-1;
		b[0]=T/2*(1+2*C+2*pow(C,2)+pow(C,3));
		b[1]=T/2*(3+2*C-2*pow(C,2)-3*pow(C,3));
		b[2]=T/2*(3-2*C-2*pow(C,2)+3*pow(C,3));
		b[3]=T/2*(1-2*C+2*pow(C,2)-pow(C,3));
	}
	if (N==4){
		a[0]=1.5;  //3*(T/2) multiplicative factor
		a[1]=1;   //2*(T/2) multiplicative factor
		a[2]=0;
		a[3]=-1;
		a[4]=-1.5;
		b[0]=T/2*(1+2.6131259*C+3.4142136*pow(C,2)+2.6131259*pow(C,3)+pow(C,4));
		b[1]=T/2*(4+2*2.6131259*C-2*2.6131259*pow(C,3)-4*pow(C,4));
		b[2]=T/2*(6*pow(C,4)-2*3.4142136*pow(C,2)+6);
		b[3]=T/2*(4-2*2.6131259*C+2*2.6131259*pow(C,3)-4*pow(C,4));
		b[4]=T/2*(1-2.6131259*C+3.4142136*pow(C,2)-2.6131259*pow(C,3)+pow(C,4));
	}
};

double FilterH::differentiator(double X){
	//assumes the filter was already run so that the coefficients are able
	//recive the new input- retunrn the filter signal
	x[4]=x[3];
	x[3]=x[2];
	x[2]=x[1];
	x[1]=x[0];
	x[0]=X;
	y[4]=y[3];
	y[3]=y[2];
	y[2]=y[1];
	y[1]=y[0];
	y[0]=(a[0]*x[0]+x[1]*a[1]+x[2]*a[2]+x[3]*a[3]+x[4]*a[4]-y[1]*b[1]-y[2]*b[2]-y[3]*b[3]-y[4]*b[4])/b[0];
	return y[0];
};
