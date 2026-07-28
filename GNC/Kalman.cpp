#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>
#include <cmath>
using namespace Eigen;

namespace mathtools {
inline constexpr double PI = 3.14159265358979323846;


inline Matrix3d skew(const Vector3d &v) {
    Matrix3d S;
    S << 0.0, -v.z(), v.y(),
         v.z(), 0.0, -v.x(),
        -v.y(), v.x(), 0.0;
    return S;
}

inline double wrapAngle(double a) {
    while (a > PI)  a -= 2.0 * PI;
    while (a < -PI) a += 2.0 * PI;
    return a;
}

}

struct IMU {
    Vector3d acc;
    Vector3d gyro;
};

class EKF {
    public:
        enum Mode {INERTIAL,TERMINAL};

        Mode mode;

        Matrix<double,15,15> P_inertial;

        Vector3d pos, vel;
        Quaterniond q;
        Vector3d a_bias, g_bias;
        double sigma_acc, sigma_gyro, sigma_a_bias, sigma_g_bias;
        Matrix<double,15,15> Qc;
        Vector3d g;

        Matrix<double,9,1> x_term;
        Matrix4d R_base_4;
        double sigma_jerk;
        double sigmaRange, sigmaAz, sigmaEl, sigmaRdot;
        MatrixXd P_term;     
        MatrixXd Q_term;     
        MatrixXd R_term;
        double dt = 0.001;

        EKF(Mode mode_ = INERTIAL) : mode(mode_) {
            if (mode == INERTIAL) initInertial();
            else initTerminal();
        }

        void initInertial() {
            pos.setZero();
            vel.setZero();
            q.setIdentity();
            a_bias.setZero();
            g_bias.setZero();

            P_inertial.setZero();
            P_inertial.diagonal() << 10.0,10.0,10.0,
                                    5.0,5.0,5.0,
                                    0.01,0.01,0.01,
                                    0.01,0.01,0.01,
                                    1e-3,1e-3,1e-3;

            sigma_acc    = 0.01551;
            sigma_gyro   = 0.0027596;
            sigma_a_bias = 1e-5;
            sigma_g_bias = 1e-6;

            Qc.setZero();
            Qc.block<3,3>(3,3)   = Matrix3d::Identity() * (sigma_acc*sigma_acc);
            Qc.block<3,3>(6,6)   = Matrix3d::Identity() * (sigma_gyro*sigma_gyro);
            Qc.block<3,3>(9,9)   = Matrix3d::Identity() * (sigma_a_bias*sigma_a_bias);
            Qc.block<3,3>(12,12) = Matrix3d::Identity() * (sigma_g_bias*sigma_g_bias);

            g = Vector3d(0.0, 0.0, -9.81);
        }

        void initTerminal() {
            x_term = VectorXd::Zero(9);
            P_term = MatrixXd::Identity(9,9);    // Large initial uncertainty
            P_term.block<3,3>(0,0) = Matrix3d::Identity() * 100.0;  // pos variance (m²)
            P_term.block<3,3>(3,3) = Matrix3d::Identity() * 100.0;   // vel variance (m²/s²)
            P_term.block<3,3>(6,6) = Matrix3d::Identity() * 50.0;

            sigma_jerk = 10.0;   // Jerk std dev (m/s³)

            // Continuous process noise (for constant jerk model)
            Q_term = MatrixXd::Zero(9,9);
            updateProcessNoise();

            // Default R (6x6 for pos+vel updates)
            R_term = MatrixXd::Zero(6,6);
            R_term.block<3,3>(0,0) = 5.0 * Matrix3d::Identity();  // 5 m std
            R_term.block<3,3>(3,3) = 1.0  * Matrix3d::Identity();  // 1 m/s std
        }

        void predictInertial(const IMU &imu, double dt) {
            Vector3d acc_corr = imu.acc - a_bias;
            Vector3d omega_corr = imu.gyro - g_bias;

            Vector3d dtheta = omega_corr * dt;
            double dtheta_norm = dtheta.norm();
            Quaterniond dq;
            if (dtheta_norm > 1e-8)
                dq = AngleAxisd(dtheta_norm, dtheta / dtheta_norm);
            else {
                dq.w() = 1.0;
                dq.vec() = 0.5 * dtheta;
                dq.normalize();
            }
            q = (q * dq).normalized();

            Matrix3d Rwb = q.toRotationMatrix();
            Vector3d a_world = Rwb * acc_corr + g;
            pos += vel * dt + 0.5 * a_world * dt * dt;
            vel += a_world * dt;

            Matrix<double,15,15> Fc = Matrix<double,15,15>::Zero();
            Fc.block<3,3>(0,3) = Matrix3d::Identity();
            Fc.block<3,3>(3,6) = -Rwb * mathtools::skew(acc_corr);
            Fc.block<3,3>(3,9) = -Rwb;
            Fc.block<3,3>(6,6) = -mathtools::skew(omega_corr);
            Fc.block<3,3>(6,12)= -Matrix3d::Identity();

            Matrix<double,15,15> Fd = Matrix<double,15,15>::Identity() + Fc*dt;
            Matrix<double,15,15> Qd = Qc*dt;

            P_inertial = Fd * P_inertial * Fd.transpose() + Qd;
        }

        void updateInertialPos(const Vector3d &z_p, const Matrix3d &Rpos) {
            Matrix<double,3,15> H = Matrix<double,3,15>::Zero();
            H.block<3,3>(0,0) = Matrix3d::Identity();
            Vector3d y = z_p - pos;

            Matrix3d S = H * P_inertial * H.transpose() + Rpos;
            Matrix<double,15,3> K = P_inertial * H.transpose() * S.inverse();
            Matrix<double,15,1> dx = K * y;
            applyInertialError(dx);
            P_inertial = (Matrix<double,15,15>::Identity() - K*H) * P_inertial;
        }

        void updateInertialAtt(const Quaterniond &z_q, const Matrix3d &Ratt) {
            Quaterniond q_err = z_q * q.conjugate();
            if (q_err.w() < 0.0) q_err.coeffs() *= -1.0;
            Vector3d dtheta_meas = 2.0 * q_err.vec();

            Matrix<double,3,15> H = Matrix<double,3,15>::Zero();
            H.block<3,3>(0,6) = Matrix3d::Identity();

            Vector3d y = dtheta_meas;
            Matrix3d S = H * P_inertial * H.transpose() + Ratt;
            Matrix<double,15,3> K = P_inertial * H.transpose() * S.inverse();
            Matrix<double,15,1> dx = K * y;
            applyInertialError(dx);
            P_inertial = (Matrix<double,15,15>::Identity() - K*H) * P_inertial;
        }

        void applyInertialError(const Matrix<double,15,1> &dx) {
            pos += dx.segment<3>(0);
            vel += dx.segment<3>(3);
            Vector3d dtheta = dx.segment<3>(6);
            Quaterniond dq; dq.w() = 1.0; dq.vec() = 0.5*dtheta; dq.normalize();
            q = (dq*q).normalized();
            a_bias += dx.segment<3>(9);
            g_bias += dx.segment<3>(12);
        }

        void predictTerminal(double dt_in) {
            dt = dt_in;
            updateProcessNoise();

            // State Transition Matrix (constant jerk)
            MatrixXd F = MatrixXd::Identity(9,9);
            F.block<3,3>(0,3) = dt * Matrix3d::Identity();
            F.block<3,3>(0,6) = 0.5 * dt * dt * Matrix3d::Identity();
            F.block<3,3>(3,6) = dt * Matrix3d::Identity();

            // Predict state
            x_term = F * x_term;

            // Predict covariance
            P_term = F * P_term * F.transpose() + Q_term;
        }

        // Add at top: using Eigen;
        void updateTerminalPosVel(const VectorXd& z, const MatrixXd& R,
                                const Vector3d& p_m, const Vector3d& v_m)
        {
            // Extract predicted target states
            Vector3d p_t = x_term.segment<3>(0);
            Vector3d v_t = x_term.segment<3>(3);

            // Predicted measurement (relative pos + vel)
            VectorXd z_pred(6);
            z_pred << (p_t - p_m), (v_t - v_m);

            // Measurement matrix H
            MatrixXd H = MatrixXd::Zero(6, 9);
            H.block<3,3>(0,0) = Matrix3d::Identity();
            H.block<3,3>(3,3) = Matrix3d::Identity();

            // Innovation
            VectorXd y = z - z_pred;

            // Innovation covariance
            MatrixXd S = H * P_term * H.transpose() + R;

            // Kalman gain
            MatrixXd K = P_term * H.transpose() * S.inverse();

            // Update state
            x_term = x_term + K * y;

            // Update covariance
            MatrixXd I = MatrixXd::Identity(9,9);
            P_term = (I - K * H) * P_term;
        }


        Vector3d getTerminalPos() const { return x_term.segment<3>(0); }
        Vector3d getTerminalVel() const { return x_term.segment<3>(3); }
        Vector3d getTerminalAcc() const { return x_term.segment<3>(6); }

        void updateProcessNoise() {
            double dt2 = dt*dt;
            double dt3 = dt2*dt;
            double dt4 = dt3*dt;
            double dt5 = dt4*dt;

            Matrix3d q;
            q << dt5/20.0, dt4/8.0, dt3/6.0,
                dt4/8.0,  dt3/3.0, dt2/2.0,
                dt3/6.0,  dt2/2.0, dt;

            Q_term.setZero();
            Q_term.block<3,3>(0,0) = q;
            Q_term.block<3,3>(3,3) = q;
            Q_term.block<3,3>(6,6) = q;

            Q_term *= sigma_jerk * sigma_jerk;
        }
};