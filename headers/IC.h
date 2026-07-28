#ifndef IC_H
#define IC_H

#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <array>
#include <string>
#include <cmath>
#include <iostream>
using namespace Eigen;

class StateModel{
    public:
        enum Mode {
            LR_HS_HEAD_ON_ACC_DIVE,
            LR_HS_HEAD_ON,
            SR_HS_TAIL_CHASE,
            SR_HS_LATERAL_CROSSING,
            LR_HS_LATERAL_CROSSING
        };

        enum Fidelity {
            TRIM,
            FULL_BODY
        };

        Mode mode;

        Matrix3d Rfb;
        Vector3d rel_pos_vec;
        Vector3d rel_vel_vec;
        Vector3d rel_acc_vec;
        Vector3d est_m_pos;
        Vector3d est_m_vel;
        Vector3d est_m_acc;
        Vector3d est_t_pos;
        Vector3d est_t_vel;
        Vector3d est_t_acc;
        Vector3d sight_line_rot_vec;
        Vector3d dem_acc_vec_body;
        Vector3d executed_body_a_vec;
        Quaterniond quat_vec;
        Vector3d acc_vec;
        std::array<double,3> euler_angles;
        double thrust;
        double alpha;
        double beta;
        Vector3d aero_forces;
        Vector3d aero_moments;
        Vector3d new_p;
        Vector3d new_omega;
        Vector3d new_v;
        // LOS range, range rate and transformation matrix
        double R;
        double R_rate;
        double** t_matrix;
        // initial time and time step
        double delta_t = .001; // s
        int steps = 100000;
        // Navigational constants
        double N2_i = 3;  // no unit
        double N3_i = 3;  // no unit
        // Longitudinal and lateral autopilot time-constants
        double taox_i = .15; // s 
        double taoy_i = .15; // s 
        double taoz_i = .15; // s 
        // Enviornmental parameters
        double latitude = 38.6270; // deg
        double longitude = -90.1994; // deg
        std::string country;
        std::string city;
        std::string weather;
        double wind_speed = 0; // mph
        double wind_direction = 0; // deg
        double air_pressure_sea_lvl = 101325; // Pa
        double air_viscocity = 1.81e-5; // Pa s
        double temp_lapse_rate = .0065; // K/m
        double temp_sea_lvl = 288.15; // K
        double thermal_expansion_coefficient = .0008; // 1/K
        double g = 0; // m/s^2
        double pi = 3.14159;
        double air_molar_mass = .0289644; // kg/mol
        double gas_constant = 8.31432;  // J/mol*K
        double dry_air_constant = 287.05; // J/kg*K
        double speed_of_sound = 343; // m/s
        // Missile body parameters
        double body_length = 3; // m
        double body_radius = .2; // m
        double body_mass = 30; // kg
        double nose_mass = 5; // kg
        double payload_mass = 10; // kg
        double total_mass = 80; // kg
        double body_cm = 1.5; // m
        double payload_cm = .8*2; // m
        double cross_area = pi*body_radius*body_radius; //m2
        // Fin actuator parameters
        double fin_base = .4;  // m 
        double fin_height = .1; // m
        double fin_area = .5*.4*.1; // m^2
        double front_fin_d = .05; // m
        double fin_mass = .3; // kg
        double max_ang = .436; // rad
        double max_torque = 50; // Nm
        double act_damping = .01;
        double act_tao = .02; // s
        double Kp = 200; // Nm/rad
        double Kd = 20; // Nms/rad
        // coefficients
        double C_f = .0003; // no unit
        double C_Lf_delta = 2;  // no unit 
        double C_H = 0.01; // no unit 
        double CD;  // no unit
        double CC;  // no unit
        double CL;  // no unit
        double Cl;  // no unit
        double Cm;  // no unit
        double Cn;  // no unit
        std::array<double,4> defl_angls; // rad,rad,rad,rad
        // Propulsion parameters
        double fuel_tank_length = 1; // m
        double fuel_tank_radius = .2; // m
        double initial_fuel_mass = 50; // kg
        double fuel_mass = 50; // kg
        double fuel_density =.8; // kg/m^3
        double fuel_v = 1; // m/s
        double exhaust_v = 5000; // m/s
        double cs_area = .5; // m^2
        std::array<double,3> f_cm;
        // Missile parameters
        double i_mx; // m
        double i_my; // m
        double i_mz; // m
        double i_mvx; // m/s
        double i_mvy; // m/s
        double i_mvz; // m/s
        double i_max; // m/s^2
        double i_may; // m/s^2
        double i_maz; // m/s^2
        double i_roll; // rad
        double i_pitch; // rad
        double i_yaw; // rad
        // Target parameters
        double i_tx; // m
        double i_ty; // m
        double i_tz; // m
        double i_tvx; // m/s
        double i_tvy; // m/s
        double i_tvz; // m/s
        double i_tax; // m/s^2
        double i_tay; // m/s^2
        double i_taz; // m/s^2

        Matrix3d transformation_Matrix_Fixed_to_Body(Quaterniond quat_vec){
            Matrix3d t_fb;
            t_fb << pow(quat_vec.w(),2)+pow(quat_vec.x(),2)-pow(quat_vec.y(),2)-pow(quat_vec.z(),2), 2*(quat_vec.x()*quat_vec.y()+quat_vec.w()*quat_vec.z()), 2*(quat_vec.x()*quat_vec.z()-quat_vec.w()*quat_vec.y()),
                    2*(quat_vec.x()*quat_vec.y()-quat_vec.w()*quat_vec.z()), pow(quat_vec.w(),2)-pow(quat_vec.x(),2)+pow(quat_vec.y(),2)-pow(quat_vec.z(),2), 2*(quat_vec.y()*quat_vec.z()+quat_vec.w()*quat_vec.x()),
                    2*(quat_vec.x()*quat_vec.z()+quat_vec.w()*quat_vec.y()), 2*(quat_vec.y()*quat_vec.z()-quat_vec.w()*quat_vec.x()), pow(quat_vec.w(),2)-pow(quat_vec.x(),2)-pow(quat_vec.y(),2)+pow(quat_vec.z(),2);
            return t_fb;
        }

        Matrix3d transformation_Matrix_Body_to_Fixed(Quaterniond quat_vec){
            Matrix3d t_bf;
            t_bf << pow(quat_vec.w(),2)+pow(quat_vec.x(),2)-pow(quat_vec.y(),2)-pow(quat_vec.z(),2), 2*(quat_vec.x()*quat_vec.y()-quat_vec.w()*quat_vec.z()), 2*(quat_vec.x()*quat_vec.z()+quat_vec.w()*quat_vec.y()),
                    2*(quat_vec.x()*quat_vec.y()+quat_vec.w()*quat_vec.z()), pow(quat_vec.w(),2)-pow(quat_vec.x(),2)+pow(quat_vec.y(),2)-pow(quat_vec.z(),2), 2*(quat_vec.y()*quat_vec.z()-quat_vec.w()*quat_vec.x()),
                    2*(quat_vec.x()*quat_vec.z()-quat_vec.w()*quat_vec.y()), 2*(quat_vec.y()*quat_vec.z()+quat_vec.w()*quat_vec.x()), pow(quat_vec.w(),2)-pow(quat_vec.x(),2)-pow(quat_vec.y(),2)+pow(quat_vec.z(),2);
            return t_bf;
        }

        StateModel() {}

        void setState(Mode mode_){
            mode = mode_;
            switch (mode){
                case Mode::LR_HS_HEAD_ON_ACC_DIVE:
                    innit_LR_HS_HEAD_ON_ACC_DIVE();
                    break;
                case Mode::LR_HS_HEAD_ON:
                    innit_LR_HS_HEAD_ON();
                    break;
                case Mode::SR_HS_TAIL_CHASE:
                    innit_SR_HS_TAIL_CHASE();
                    break;
                case Mode::SR_HS_LATERAL_CROSSING:
                    innit_SR_HS_LATERAL_CROSSING();
                    break;
                case Mode::LR_HS_LATERAL_CROSSING:
                    innit_LR_HS_LATERAL_CROSSING();
                    break;
                default:
                    std::cout << "Unknown State\n";
                    break;
            }
        }

        void innit_state_vectors(){
            rel_pos_vec = {};
            rel_vel_vec = {};
            rel_acc_vec = {};
            est_m_pos = {};
            est_m_vel = {};
            est_m_acc = {};
            est_t_pos = {};
            est_t_vel = {};
            est_t_acc = {};
            sight_line_rot_vec = {};
            dem_acc_vec_body = {};
            executed_body_a_vec = {};
            quat_vec = {1.0,0.0,0.0,0.0};
            defl_angls = {0.0,0.0,0.0,0.0};
            acc_vec = {};
            euler_angles = {};
            thrust = 0;
            alpha = 0;
            beta = 0;
            f_cm = {};
            aero_forces = {};
            aero_moments = {};
            new_p = {};
            new_omega = {};
            new_v = {};
        }

        void innit_LR_HS_HEAD_ON_ACC_DIVE(){
            innit_state_vectors();

            R = -1;

            i_mx = 0;
            i_my = 0; 
            i_mz = 5000; 
            i_mvx = 800;
            i_mvy = 0; 
            i_mvz = 0;
            i_max = 0; 
            i_may = 0; 
            i_maz = 0; 
            i_roll = 0; 
            i_pitch = 0; 
            i_yaw = 0; 

            i_tx = 50000; 
            i_ty = 0; 
            i_tz = 5000;
            i_tvx = -230; 
            i_tvy = 0; 
            i_tvz = -10; 
            i_tax = 0; 
            i_tay = 0; 
            i_taz = 0; 
        }

        void innit_LR_HS_HEAD_ON(){
            innit_state_vectors();

            R = -1;

            i_mx = 0;
            i_my = 0; 
            i_mz = 5000; 
            i_mvx = 600;
            i_mvy = 0; 
            i_mvz = 0;
            i_max = 0; 
            i_may = 0; 
            i_maz = 0; 
            i_roll = 0; 
            i_pitch = 0; 
            i_yaw = 0; 

            i_tx = 50000; 
            i_ty = 0; 
            i_tz = 5000;
            i_tvx = -230; 
            i_tvy = 0; 
            i_tvz = 0; 
            i_tax = 0; 
            i_tay = 0; 
            i_taz = 0;
        }

        void innit_SR_HS_TAIL_CHASE(){
            innit_state_vectors();

            R = -1;

            i_mx = 0;
            i_my = 0; 
            i_mz = 1000; 
            i_mvx = 600;
            i_mvy = 0; 
            i_mvz = 0;
            i_max = 0; 
            i_may = 0; 
            i_maz = 0; 
            i_roll = 0; 
            i_pitch = 0; 
            i_yaw = 0; 

            i_tx = 10000; 
            i_ty = 50; 
            i_tz = 1000;
            i_tvx = 220; 
            i_tvy = 20; 
            i_tvz = 20; 
            i_tax = 0; 
            i_tay = 10; 
            i_taz = 0;
        }

        void innit_LR_HS_LATERAL_CROSSING(){
            innit_state_vectors();

            R = -1;

            i_mx = 0;
            i_my = 0; 
            i_mz = 3000; 
            i_mvx = 500;
            i_mvy = 0; 
            i_mvz = 0;
            i_max = 0; 
            i_may = 0; 
            i_maz = 0; 
            i_roll = 0; 
            i_pitch = 0; 
            i_yaw = 0; 

            i_tx = 20000; 
            i_ty = -4000; 
            i_tz = 4000;
            i_tvx = 0; 
            i_tvy = 200; 
            i_tvz = -10; 
            i_tax = -5; 
            i_tay = 0; 
            i_taz = 5;
        }

        void innit_SR_HS_LATERAL_CROSSING(){
            innit_state_vectors();

            R = -1;

            i_mx = 0;
            i_my = 0; 
            i_mz = 3000; 
            i_mvx = 500;
            i_mvy = 0; 
            i_mvz = 0;
            i_max = 0; 
            i_may = 0; 
            i_maz = 0; 
            i_roll = 0; 
            i_pitch = 0; 
            i_yaw = 0; 

            i_tx = 3000; 
            i_ty = -1000; 
            i_tz = 3000;
            i_tvx = 0; 
            i_tvy = 300; 
            i_tvz = 0; 
            i_tax = 0; 
            i_tay = 0; 
            i_taz = -5;
        }
};

#endif