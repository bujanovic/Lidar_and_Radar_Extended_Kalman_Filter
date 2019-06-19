#include <math.h>
#include <iostream>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * Predict the state
   */
   x_ = F_ * x_;
   
   P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * Update the state with Lidar by using Kalman Filter equations
   */
   
   VectorXd y_k = z - H_ * x_;
   MatrixXd S_k = H_ * P_ * H_.transpose() + R_;
   MatrixXd K_k = P_ * H_.transpose() * S_k.inverse();
   x_ = x_ + K_k * y_k;
   MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
   P_ = (I - K_k*H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * Update the state with Radar by using Extended Kalman Filter equations
   */
   VectorXd H_x = VectorXd(3);
   float px, py, vx, vy, rho, phi, rho_dot;
   px = x_[0];
   py = x_[1];
   vx = x_[2];
   vy = x_[3];
   rho = sqrt(px*px + py*py);
   phi = atan2(py, px);  // returns values between -pi and pi

   //Make sure we're not dividing by zero
   if (rho < 0.000001){
	rho = 0.000001;
   }
   rho_dot = (px*vx+py*vy)/rho;
   H_x  << rho, phi, rho_dot;
   
   VectorXd y_k = z - H_x;

   // normalize the angle between -pi to pi
 
   while(y_k(1) > M_PI){
    y_k(1) -= 2*M_PI;
   }

   while(y_k(1) < -M_PI){
    y_k(1) += 2*M_PI;
   }
     
   MatrixXd S_k = H_ * P_ * H_.transpose() + R_;
   MatrixXd K_k = P_ * H_.transpose() * S_k.inverse();
   x_ = x_ + K_k * y_k;
   MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
   P_ = (I - K_k*H_) * P_;
   
}
