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
        Fidelity fid;

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
        double N1_i = 0;  // no unit
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
        double g = -9.81; // m/s^2
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
        double fin_area = .5*fin_base*fin_height; // m^2
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

        void setState(){
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

        void readConfig(const std::string& cfgfile){
            std::ifstream file(cfgfile);

            if (!file)
            {
                throw std::runtime_error("Could not open config file: " + cfgfile);
            }

            std::string heading;
            std::string value;

            while (std::getline(file, heading))
            {
                if (!std::getline(file, value))
                {
                    throw std::runtime_error(
                        "Missing value for setting: " + heading
                    );
                }

                else if (heading == "FIDELITY")
                {
                    if (value == "trim") fid = TRIM;
                    if (value == "full_body") fid = FULL_BODY;
                }

                else if (heading == "ENGAGEMENT TYPE")
                {
                    if (value == "lr_hs_head_on_acc_dive") mode = LR_HS_HEAD_ON_ACC_DIVE;
                    if (value == "lr_hs_head_on") mode = LR_HS_HEAD_ON;
                    if (value == "sr_hs_tail_chase") mode = SR_HS_TAIL_CHASE;
                    if (value == "sr_hs_lateral_crossing") mode = SR_HS_LATERAL_CROSSING;
                    if (value == "lr_hs_lateral_crossing") mode = LR_HS_LATERAL_CROSSING;
                }

                else if (heading == "TIME STEP")
                {
                    delta_t = std::stod(value);
                }

                else if (heading == "MAX STEPS")
                {
                    steps = std::stod(value);
                }

                else if (heading == "X NAVIGATIONAL CONSTANT")
                {
                    N1_i = std::stod(value);
                }

                else if (heading == "Y NAVIGATIONAL CONSTANT")
                {
                    N2_i = std::stod(value);
                }

                else if (heading == "Z NAVIGATIONAL CONSTANT")
                {
                    N3_i = std::stod(value);
                }

                else if (heading == "X LONGITUDINAL AUTOPILOT TIME CONSTANT")
                {
                    taox_i = std::stod(value);
                }

                else if (heading == "Y LATERAL AUTOPILOT TIME CONSTANT")
                {
                    taoy_i = std::stod(value);
                }

                else if (heading == "Z LATERAL AUTOPILOT TIME CONSTANT")
                {
                    taoz_i = std::stod(value);
                }

                else if (heading == "LATITUDE")
                {
                    latitude = std::stod(value);
                }

                else if (heading == "LONGITUDE")
                {
                    longitude = std::stod(value);
                }

                else if (heading == "SEA LEVEL AIR PRESSURE")
                {
                    air_pressure_sea_lvl = std::stod(value);
                }

                else if (heading == "AIR VISCOCITY")
                {
                    air_viscocity = std::stod(value);
                }

                else if (heading == "TEMPERATURE LAPSE RATE")
                {
                    temp_lapse_rate = std::stod(value);
                }

                else if (heading == "SEA LEVEL TEMPERATURE")
                {
                    temp_sea_lvl = std::stod(value);
                }

                else if (heading == "THERMAL EXPANSION COEFFICIENT")
                {
                    thermal_expansion_coefficient = std::stod(value);
                }

                else if (heading == "GRAVITY")
                {
                    g = std::stod(value);
                }

                else if (heading == "MOLAR MASS OF AIR")
                {
                    air_molar_mass = std::stod(value);
                }

                else if (heading == "GAS CONSTANT")
                {
                    gas_constant = std::stod(value);
                }

                else if (heading == "DRY AIR CONSTANT")
                {
                    dry_air_constant = std::stod(value);
                }

                else if (heading == "SPEED OF SOUND")
                {
                    speed_of_sound = std::stod(value);
                }

                else if (heading == "BODY LENGTH")
                {
                    body_length = std::stod(value);
                }

                else if (heading == "BODY RADIUS")
                {
                    body_radius = std::stod(value);
                }

                else if (heading == "BODY MASS")
                {
                    body_mass = std::stod(value);
                }

                else if (heading == "NOSE MASS")
                {
                    nose_mass = std::stod(value);
                }

                else if (heading == "PAYLOAD MASS")
                {
                    payload_mass = std::stod(value);
                }

                else if (heading == "TOTAL MASS")
                {
                    total_mass = std::stod(value);
                }

                else if (heading == "BODY CENTER OF MASS LOCATION")
                {
                    body_cm = std::stod(value);
                }

                else if (heading == "PAYLOAD CENTER OF MASS LOCATION")
                {
                    payload_cm = std::stod(value);
                }

                else if (heading == "FIN BASE LENGTH")
                {
                    fin_base = std::stod(value);
                }

                else if (heading == "FIN HEIGHT")
                {
                    fin_height = std::stod(value);
                }

                else if (heading == "FRONT FIN DISTANCE")
                {
                    front_fin_d = std::stod(value);
                }

                else if (heading == "FIN MASS")
                {
                    fin_mass = std::stod(value);
                }

                else if (heading == "MAX FIN ANGLE")
                {
                    max_ang = std::stod(value);
                }

                else if (heading == "MAX ACTUATOR TORQUE")
                {
                    max_torque = std::stod(value);
                }

                else if (heading == "ACTUATOR DAMPING")
                {
                    act_damping = std::stod(value);
                }

                else if (heading == "ACTUATOR TIME CONSTANT")
                {
                    act_tao = std::stod(value);
                }

                else if (heading == "KP")
                {
                    Kp = std::stod(value);
                }

                else if (heading == "KD")
                {
                    Kd = std::stod(value);
                }

                else if (heading == "FIN COEFFICIENT")
                {
                    C_f = std::stod(value);
                }

                else if (heading == "COEFFICIENT OF LIFT FROM FIN DEFLECTION")
                {
                    C_Lf_delta = std::stod(value);
                }

                else if (heading == "HINGE MOMENT COEFFICIENT")
                {
                    C_H = std::stod(value);
                }

                else if (heading == "FUEL TANK LENGTH")
                {
                    fuel_tank_length = std::stod(value);
                }

                else if (heading == "FUEL TANK RADIUS")
                {
                    fuel_tank_radius = std::stod(value);
                }

                else if (heading == "INITIAL FUEL MASS")
                {
                    initial_fuel_mass = std::stod(value);
                }

                else if (heading == "FUEL MASS")
                {
                    fuel_mass = std::stod(value);
                }

                else if (heading == "FUEL DENSITY")
                {
                    fuel_density = std::stod(value);
                }

                else if (heading == "FUEL VELOCITY")
                {
                    fuel_v = std::stod(value);
                }

                else if (heading == "EXHAUST VELOCITY")
                {
                    exhaust_v = std::stod(value);
                }

                else if (heading == "JET CROSS SECTIONAL AREA")
                {
                    cs_area = std::stod(value);
                }
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