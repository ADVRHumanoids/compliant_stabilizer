#include <compliant_stabilizer/compliantstabilizer.h>


CompliantStabilizer::CompliantStabilizer(const double sample_time, const double mass,
                                         const double ankle_height,
                                         const Eigen::Vector2d& foot_size,
                                         const double Fzmin,
                                         const Eigen::Vector3d& K, const Eigen::Vector3d& C,
                                         const Eigen::Vector3d& MaxLims,
                                         const Eigen::Vector3d& MinLims,
                                         const double samples2ODE,
                                         const double freq):
    ANKLE_X_OFFSET(DEFAULT_ANKLE_X_OFFSET),
    ANKLE_HEIGHT(ankle_height),
    FOOT_LENGTH(foot_size[0]),
    FOOT_WIDTH(foot_size[1])
{

    FzODE=Vector2d::Zero();

    m_sampletime=sample_time;
    m_g= 9.81;
    m_mass=mass; //30
    m_Fzmin= Fzmin;


    this->m_Kgains=K;
    this->m_Cgains=C;
    m_MaxLimit<<MaxLims;
    m_MinLimit<<MinLims;

    m_samples2ODE=samples2ODE;
    m_freq=freq;

    m_dZMP_bufferODE.resize(2,m_samples2ODE);
    m_cop=Vector2d::Zero();

    m_FilterCOPX.butterworth ( this->m_sampletime, m_freq, 1 );
    m_FilterCOPY.butterworth ( this->m_sampletime, m_freq, 1 );


}

CompliantStabilizer::CompliantStabilizer(const double sample_time, const double mass,
                                         const double ankle_height,
                                         const Eigen::Vector2d& foot_size):
    ANKLE_X_OFFSET(DEFAULT_ANKLE_X_OFFSET),
    ANKLE_HEIGHT(ankle_height),
    FOOT_LENGTH(foot_size[0]),
    FOOT_WIDTH(foot_size[1])
{

    FzODE=Vector2d::Zero();

    m_sampletime=sample_time;
    m_g= 9.81;
    m_mass=mass; //30
    m_Fzmin= DEFAULT_Fzmin;


    this->m_Kgains={DEFAULT_Kx,DEFAULT_Ky,0};
    this->m_Cgains={DEFAULT_Cx,DEFAULT_Cy,0};
    m_MaxLimit<<DEFAULT_MaxLimsx,DEFAULT_MaxLimsy,DEFAULT_MaxLimsz;
    m_MinLimit<<DEFAULT_MinLimsx,DEFAULT_MinLimsy,DEFAULT_MinLimsz;

    m_samples2ODE=DEFAULT_smaples2ODE;
    m_freq=DEFAULT_freq;

    m_dZMP_bufferODE.resize(2,m_samples2ODE);
    m_cop=Vector2d::Zero();

    m_FilterCOPX.butterworth ( this->m_sampletime, m_freq, 1 );
    m_FilterCOPY.butterworth ( this->m_sampletime, m_freq, 1 );


}

CompliantStabilizer::~CompliantStabilizer(){

}

void CompliantStabilizer::setAnkleOffset(const double x_offset)
{
    ANKLE_X_OFFSET = x_offset;
}

Vector3d CompliantStabilizer::update(Eigen::VectorXd FTLeft, Eigen::VectorXd FTRight,
                                     Eigen::Vector2d CoPLeft, Eigen::Vector2d CoPRight,
                                     Vector3d Lft, Vector3d Rft){


    CalcCop(FTRight, FTLeft, Rft, Lft);
    double Gmg= m_mass*m_g;
    static Vector3d deltaHip_ODE;
    double Fzl, Fzr;
    Vector2d cop_delta( 0, 0 );

    double sum0ode=0;
    double sumode[4]={0};

    Vector2d deltaZMP_L= cop_in_lft_raw.head(2) - CoPLeft;
    Vector2d deltaZMP_R= cop_in_rft_raw.head(2)- CoPRight;

    Fzl = FTLeft(2);
    Fzr = FTRight(2);
    if( Fzl< m_Fzmin ){
        Fzl= m_Fzmin;
    }
    if( Fzr< m_Fzmin ){
        Fzr= m_Fzmin;
    }

    cop_delta = deltaZMP_L * Fzl/( Fzl+Fzr ) + deltaZMP_R * Fzr/( Fzl+Fzr );

    double deltaZMPx_old= m_cop(0);
    double deltaZMPy_old= m_cop(1);
    double deltaZMPx_ODE= m_FilterCOPX.applyFilter( cop_delta( 0 ) );
    double deltaZMPy_ODE= m_FilterCOPY.applyFilter( cop_delta( 1 ) );



    /*derevative of filtered COP */
    for(int i=0;i<m_samples2ODE-1;i++){
        m_dZMP_bufferODE(0,i)=m_dZMP_bufferODE(0,i+1);
        m_dZMP_bufferODE(1,i)=m_dZMP_bufferODE(1,i+1);
    }//this shifts the data in array
    m_dZMP_bufferODE(0,m_samples2ODE-1)=(deltaZMPx_ODE-deltaZMPx_old)/this->m_sampletime;
    m_dZMP_bufferODE(1,m_samples2ODE-1)=(deltaZMPy_ODE-deltaZMPy_old)/this->m_sampletime;
    // begin of ZMP derevative filtering

    for (int j=0;j<2;j++){
        for(int i=m_samples2ODE;i>0;i--){
            sum0ode+=i;
            sumode[j] += i*m_dZMP_bufferODE(j,i-1);
        }
        m_dZMPODE[j]=sumode[j]/sum0ode;
    }

    m_cop<<deltaZMPx_ODE,deltaZMPy_ODE;

    FzODE[0]+= 1*( Fzl-m_Fzmin -FzODE[0] )*this->m_sampletime;//left
    FzODE[1]+= 1*( Fzr-m_Fzmin -FzODE[1] )*this->m_sampletime;//right

    //2.5 not 1

    double n1= ( FzODE[0]+FzODE[1] )/( Gmg*0.9 );//normalized unit based on GRF/mg

    if( n1>1.2 )     n1= 1.2;
    else if( n1<0 )    n1= 0;
    // Must normalize!!

    //---- stabilizer law: gain must < 1 to have negative feedback control
    double Kx= this->m_Kgains(0)*n1;//m_Kx= 0.2;

    double Ky= this->m_Kgains(1)*n1;//m_Ky= 0.4;
    double Cx= this->m_Cgains(0)*n1;//m_Cx=-0.01
    double Cy= this->m_Cgains(1)*n1;//m_Cy=-0.02

    /*-------------- here is the main law  -------------------*/
    deltaHip_ODE( 0 ) = Kx*deltaZMPx_ODE +Cx*m_dZMPODE[0];//delta hip x
    deltaHip_ODE( 1 ) = Ky*deltaZMPy_ODE +Cy*m_dZMPODE[1];//delta hip y

    double Fext= 1-( Fzl+Fzr )/Gmg;

    double scaleCOP=0;

    if ( ( Fzl+Fzr ) > 2*m_Fzmin )//normalized unit based on GRF/mg
    {
        scaleCOP += 5*( 1-scaleCOP )*this->m_sampletime;//if load on the ground, coefficient is 1
    }
    else
    {
        scaleCOP += 5*( 0-scaleCOP )*this->m_sampletime;// otherwise it is in the air, coefficient is 0
    }
    deltaHip_ODE( 2 ) = 0;//scaleCOP*( Fext*this->m_sampletime + this->m_Cgains( 2 )*deltaHip_ODE( 2 ) )/( this->m_Kgains(2)*this->m_sampletime+this->m_Cgains(2) );//delta hip z

    //evaluate min/max limits
    for(int i=0; i<=2; i++){
        if ( deltaHip_ODE[i]>this->m_MaxLimit(i) )
            deltaHip_ODE[i]=this->m_MaxLimit(i);

        if ( deltaHip_ODE[i]<this->m_MinLimit(i) )
            deltaHip_ODE[i]=this->m_MinLimit(i);
    }

    return deltaHip_ODE;

}
void CompliantStabilizer::setGains(double Kx, double Ky, double Cx, double Cy){
    this->m_Kgains<<Kx,Ky,0;
    this->m_Cgains<<Cx,Cy,0;

}

void CompliantStabilizer::CalcCop(VectorXd FT_foot_right, VectorXd FT_foot_left, Vector3d Rft, Vector3d Lft )
{
    /*	for COP calculation */
    double pxrDS=0;
    double pyrDS=0;
    double pxlDS=0;
    double pylDS=0;
    double pxDS=0;
    double pyDS=0;
    double pzDS=ANKLE_HEIGHT;//0.02;

    //begin computation of ZMP based on force without filtering
    // right foot co
    if(FT_foot_right[2]<0.1) FT_foot_right[2]=0.1;

    if ( FT_foot_right(2)>m_Fzmin)
    {
        pxrDS=(- FT_foot_right(4)-pzDS* FT_foot_right(0))/ FT_foot_right[2];
        pyrDS=( FT_foot_right(3)-pzDS* FT_foot_right(1))/ FT_foot_right[2];
        if (pxrDS > (0.5*FOOT_LENGTH+ANKLE_X_OFFSET))
        {
            pxrDS=(0.5*FOOT_LENGTH+ANKLE_X_OFFSET);
        }
        else if (pxrDS < -(0.5*FOOT_LENGTH-ANKLE_X_OFFSET))
        {
            pxrDS= -(0.5*FOOT_LENGTH-ANKLE_X_OFFSET);
        }
        if (pyrDS > (0.5*FOOT_WIDTH))
        {
            pyrDS=(0.5*FOOT_WIDTH);
        }
        else if(pyrDS < -(0.5*FOOT_WIDTH))
        {
            pyrDS=-(0.5*FOOT_WIDTH);
        }
    }
    else
    {
        pxrDS=0;
        pyrDS=0;
    }
    // left foot cop
    if(FT_foot_left[2]<0.1) FT_foot_left[2]=0.1;
    if ( FT_foot_left(2)>m_Fzmin)
    {
        pxlDS=(- FT_foot_left(4)-pzDS* FT_foot_left(0))/ FT_foot_left[2];
        pylDS=( FT_foot_left(3)-pzDS* FT_foot_left(1))/FT_foot_left[2];
        if (pxlDS > (0.5*FOOT_LENGTH+ANKLE_X_OFFSET))
        {
            pxlDS=(0.5*FOOT_LENGTH+ANKLE_X_OFFSET);
        }
        else if (pxlDS < -(0.5*FOOT_LENGTH-ANKLE_X_OFFSET))
        {
            pxlDS= -(0.5*FOOT_LENGTH-ANKLE_X_OFFSET);
        }
        if (pylDS > (0.5*FOOT_WIDTH))
        {
            pylDS=(0.5*FOOT_WIDTH);
        }
        else if(pylDS < -(0.5*FOOT_WIDTH))
        {
            pylDS=-(0.5*FOOT_WIDTH);
        }
    }
    else
    {
        pxlDS=0;
        pylDS=0;
    }

    //end of computation of ZMP based on force without filtering

     lcop_raw << pxlDS,pylDS,0;
     rcop_raw << pxrDS,pyrDS,0;

     double Fz_ratio_l =FT_foot_left[2]/(FT_foot_left[2]+FT_foot_right[2]);
     double Fz_ratio_r = FT_foot_right[2]/(FT_foot_left[2]+FT_foot_right[2]);

     cop_in_lft_raw = Fz_ratio_l * lcop_raw + Fz_ratio_r * (rcop_raw + Rft-Lft);
    cop_in_rft_raw = Fz_ratio_l * (lcop_raw + Lft-Rft) + Fz_ratio_r * rcop_raw;

}

