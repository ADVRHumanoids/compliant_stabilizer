
/**
 * License HERE
*/

/// \file       CompliantStabilizer.hpp
/// \brief      CompliantStabilizer that provides a COM trajectory correction
///             to the provided reference
///
/// \details    This controller takes the CoP and the desired ZMP position and provides a
///             COM motion correction to stabilize the robot.
///             This controller is a child of FeedbackController
///
/// \see        FeedbackController.hpp
///
///
/// \authors    chengxu zhou and Juan Alejandro Castano (juan.castano@iit.it)
/// \date        2017
/// \version     1.0
/// \copyright   GNU Public License
///

#ifndef COMPLIANTSTABILIZER_HPP
#define COMPLIANTSTABILIZER_HPP
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include "FilterHClass.h"
#include <compliant_stabilizer/defines.h>
#include <XBotInterface/MatLogger.hpp>

using namespace Eigen;

namespace Eigen {
    typedef Eigen::Matrix<double, 6, 1> Vector6d;
}

class CompliantStabilizer {
    ///
    /// \brief The default constructor initialize controller name and use robot control states
    ///
public:
    CompliantStabilizer(const double sample_time, const double mass,
                        const double ankle_height,
                        const Eigen::Vector2d& foot_size,
                         const double Fzmin,
                         const Eigen::Vector3d& K, const Eigen::Vector3d& C,
                         const Eigen::Vector3d& MaxLims,
                         const Eigen::Vector3d& MinLims,
                         const double samples2ODE=DEFAULT_smaples2ODE,
                         const double freq=DEFAULT_freq);

    CompliantStabilizer(const double sample_time, const double mass,
                        const double ankle_height,
                        const Eigen::Vector2d& foot_size);
    /// \brief default destructor
    ~CompliantStabilizer();

    /// \brief Apply a control action to the robot's states
    ///        according with the reference states
    /// \param[in] robotState   Present state of the robot
    /// \param[in] originalControlState Present control states to be modified by the controllers
    /// \param[out] modifiedControlState states after modification
 Vector3d update(const Eigen::Vector6d& FTLeft, const Eigen::Vector6d& FTRight,
                 const Eigen::Vector2d& CoPLeft, const Eigen::Vector2d& CoPRight,
                 const Vector3d& Lft, const Vector3d& Rft);



    /// \brief Get the gains information of the specific controller.
    /// \return A map composed of the gains name and it's value/s.
   void setGains(double Kx, double Ky, double Cx, double Cy);

   void setAnkleOffset(const double x_offset);
   

   const Eigen::Vector3d& getCoP() const;
   const Eigen::Vector3d& getCoP_scaled() const;
   
   void flush();

private:
   void CalcCop(const Eigen::Vector6d& FT_foot_right, 
                const Eigen::Vector6d& FT_foot_left, 
                const Vector3d& Rft, 
                const Vector3d& Lft);

    MatrixXd m_dZMP_bufferODE; //!< Derivative of filtered ZMP circular buffer
    int m_samples2ODE; //!< circular buffer size for the ZMP derivation
    double m_dZMPODE[2]={0}; //!< Filtered ZMP for x and y
    Vector2d m_cop;          //!< Filter delta CoP
    FilterH  m_FilterCOPX, m_FilterCOPY; //!< Butterworth COP filters

    double m_sampletime;  //!< Sampling time
    double m_freq; //!< Cutoff filters frequency

    double m_g; //!< Gravity
    double m_mass; //!< Robot mass
    double m_Fzmin; //!< Minimum force detected

    Vector2d FzODE;

    Vector3d m_Kgains; //<! Proportional control gains
    Vector3d m_Cgains; //<! Derivative control gains
    Vector3d m_MaxLimit; //<! Upper bound of the ZMP
    Vector3d m_MinLimit; //<! Lower bound of the ZMP
    double ANKLE_HEIGHT; //0.05;
    double FOOT_LENGTH; //0.2;
    double ANKLE_X_OFFSET;
    double FOOT_WIDTH; //0.1;
    Vector3d lcop_raw;
    Vector3d rcop_raw;

    Vector3d cop_in_lft_raw;
    Vector3d cop_in_rft_raw;

    XBot::MatLogger::Ptr _logger;

    Vector3d _measured_cop_raw;
    Vector3d _measured_cop_scaled;
};





#endif // COMPLIANTSTABILIZER_HPP

