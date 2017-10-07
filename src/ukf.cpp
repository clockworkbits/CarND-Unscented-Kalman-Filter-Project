#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 4.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.1 * M_PI;;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    weights_(i) = 1 / (2 * (lambda_ + n_aug_));
  }

  use_laser_ = false;
  use_radar_ = true;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_)
      || (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_)) {
    return;
    cout << "Return" << endl;
  } else {
    cout << "Process" << endl;
  }

  // Initialization

  if (!is_initialized_) {
    // first measurement
    x_ = VectorXd(5);
    x_ << 1, 1, 1, 1, 1;

    cout << "Init" << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      const double ro = meas_package.raw_measurements_(0);
      const double theta = meas_package.raw_measurements_(1);

      x_(0) = ro * cos(theta);
      x_(1) = ro * sin(theta);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }

    time_us_ = meas_package.timestamp_;

    P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  // Process measurement

  const double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  // Prediction
  Prediction(dt);

  // Update
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;

  // TODO: move it to a member variable
  MatrixXd Q = MatrixXd(2,2);
  Q << std_a_ * std_a_, 0,
      0, std_yawdd_ * std_yawdd_;

  P_aug.bottomRightCorner(Q.rows(), Q.cols()) = Q;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  // Predict sigma points

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  for (int i = 0; i < Xsig_aug.cols(); i++) {
    VectorXd aug_state = Xsig_aug.col(i);
    double p_x = aug_state(0);
    double p_y = aug_state(1);
    double v = aug_state(2);
    double psi = aug_state(3);
    double psi_dot = aug_state(4);
    double ni_a = aug_state(5);
    double ni_psi = aug_state(6);

    VectorXd pred_state = VectorXd(n_x_);

    if (psi_dot < 0.00001 && psi_dot > -0.00001) { // Psi dot is zero
      pred_state(0) = p_x + v * cos(psi) * delta_t;
      pred_state(1) = p_y + v * sin(psi) * delta_t;
    } else {
      pred_state(0) = p_x + v/psi_dot * (sin(psi + psi_dot * delta_t) - sin(psi));
      pred_state(1) = p_y + v/psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi));
    }
    pred_state(2) = v;
    pred_state(3) = psi + psi_dot * delta_t;
    pred_state(4) = psi_dot;

    // noise
    pred_state(0) += 0.5 * delta_t * delta_t * cos(psi) * ni_a;
    pred_state(1) += 0.5 * delta_t * delta_t * sin(psi) * ni_a;
    pred_state(2) += ni_a * delta_t;
    pred_state(3) += 0.5 * delta_t * delta_t * ni_psi;
    pred_state(4) += delta_t * ni_psi;


    Xsig_pred_.col(i) = pred_state;
  }

  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  //predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = NormalizeAngle(x_diff(3));
//    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
//    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  cout << "Prediction done" << x_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //transform sigma points into measurement space
  const int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  for(int i = 0; i < Xsig_pred_.cols(); i++) {
    VectorXd sigma = Xsig_pred_.col(i);
    double p_x = sigma(0);
    double p_y = sigma(1);
    double v = sigma(2);
    double psi = sigma(3);
    double psi_dot = sigma(4);

    VectorXd meas_point = VectorXd(n_z);

    meas_point(0) = sqrt(p_x * p_x + p_y * p_y);
    meas_point(1) = atan2(p_y, p_x);
    meas_point(2) = (p_x*cos(psi)*v + p_y*sin(psi)*v ) / sqrt(p_x*p_x + p_y*p_y);

    Zsig.col(i) = meas_point;
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2* n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = NormalizeAngle(z_diff(1));
//    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
//    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R <<    std_radr_*std_radr_, 0, 0,
      0, std_radphi_*std_radphi_, 0,
      0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = NormalizeAngle(z_diff(1));
//    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
//    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = NormalizeAngle(x_diff(3));
//    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
//    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  //residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = NormalizeAngle(z_diff(1));
//  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
//  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  cout << "Update done" << x_ << endl;
}

double UKF::NormalizeAngle(const double angle) {
  if (angle > two_pi_) {
    return NormalizeAngle(angle - two_pi_);
  } else if (angle < -two_pi_) {
    return NormalizeAngle(angle + two_pi_);
  } else {
    return angle;
  }
}
