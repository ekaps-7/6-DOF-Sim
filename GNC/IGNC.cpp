#pragma once 

#include <csv.h>
#include <ctime>
#include <tuple>
#include <cstdio>
#include <cmath>
#include <regex>
#include <random>
#include <IC.h>
#include <gtest/gtest_prod.h>

#include <Kalman.cpp>
#include <aero_coeff.cpp>

class Enviornment: public StateModel{
    public:
        void current_weather(double latitude,double longitude){
            char cmd[512];
            snprintf(cmd,sizeof(cmd),
                    "C:\\Users\\ericr\\AppData\\Local\\Programs\\Python\\Python313\\python.exe api.py %.6f %.6f",latitude,longitude);
            FILE* weather_data = popen(cmd,"r");
            if (!weather_data){
                std::cerr << "Failed to get weather data." << std::endl;
            }
            char buffer[100];
            std::string result = "";
            while (fgets(buffer,sizeof(buffer),weather_data) != nullptr){
                result += buffer;
            }
            pclose(weather_data);

            std::regex country_pattern("^[A-Z]{2}");
            std::regex city_pattern(",([^,]+),");
            std::regex weather_pattern(",\\[([^,]+)\\],");
            std::regex pressure_pattern("\\d{4}");
            std::regex wind_speed_pattern("\\d{1,}\\.\\d{2}");
            std::regex wind_direction_pattern("[^,]*$");
            std::smatch match;

            if (std::regex_search(result,match,country_pattern)){
                country = match[0];
            }
            if (std::regex_search(result,match,city_pattern)){
                city = match[1];
            }
            if (std::regex_search(result,match,weather_pattern)){
                weather = match[1];
            }
            if (std::regex_search(result,match,pressure_pattern)){
                std::string pressure = match[0];
                air_pressure_sea_lvl = std::stod(pressure)*100;
            }
            if (std::regex_search(result,match,wind_speed_pattern)){
                std::string speed = match[0];
                wind_speed = std::stod(speed);
            }
            if (std::regex_search(result,match,wind_direction_pattern)){
                std::string direction = match[0];
                wind_direction = std::stod(direction);
            }
        }

        double temperature_at_altitude(double mz){
            double temp = temp_sea_lvl - (temp_lapse_rate*mz);
            return temp;
        }

        double fuel_density_at_altitude(double temp){
            double fuel_dens_at_t = fuel_density*(1 - thermal_expansion_coefficient*(temp-temp_sea_lvl));
            return fuel_dens_at_t;
        }

        double air_pressure_at_altitude(double air_pressure_sea_lvl,double mz){
            double air_pres_at_alt = air_pressure_sea_lvl*pow((1-(temp_lapse_rate*mz)/temp_sea_lvl),(g*air_molar_mass)/(gas_constant*temp_lapse_rate));
            return air_pres_at_alt;
        }

        double air_density_at_altitude(double air_pressure,double air_temp){
            double air_density_at_alt = air_pressure/(dry_air_constant*air_temp);
            return air_density_at_alt;
        }
};

class Fin_Actuator: public Enviornment{
    private:
        double cap_delfection(double act_ang,double max_ang,double min_ang){
            return std::max(min_ang, std::min(act_ang, max_ang));
        }

        double control_effectiveness(double mz,double mvx,double mvy,double mvz){
            double v = sqrt(pow(mvx,2)+pow(mvy,2)+pow(mvz,2));
            double air_t = temperature_at_altitude(mz);
            double air_p = air_pressure_at_altitude(air_pressure_sea_lvl,mz);
            double air_den = air_density_at_altitude(air_p,air_t);
            double q = .5 * air_den * pow(v,2);
            double k = ((q * fin_area) / total_mass) * C_Lf_delta;
            return k;
        }

        std::array<double,3> missile_velocity_body_frame(double mvx,double mvy,double mvz,Quaterniond quat_vec){
            Matrix3d Rfb = transformation_Matrix_Fixed_to_Body(quat_vec);
            Vector3d mv = Vector3d(mvx,mvy,mvz);
            Vector3d v_b = Rfb*mv;
            return {v_b.x(),v_b.y(),v_b.z()};
        }

        double angle_of_attack(double mv_x_b,double mv_z_b){
            double alpha = atan2(mv_z_b,mv_x_b);
            return alpha;
        }

        double side_slip(double mv_x_b,double mv_y_b){
            double beta = atan2(mv_y_b,mv_x_b);
            return beta;
        }

        double hinge_moment(double mz,double mvx,double mvy,double mvz,double delta){
            double v = sqrt(pow(mvx,2)+pow(mvy,2)+pow(mvz,2));
            double air_t = temperature_at_altitude(mz);
            double air_p = air_pressure_at_altitude(air_pressure_sea_lvl,mz);
            double air_den = air_density_at_altitude(air_p,air_t);
            double q = .5 * air_den * pow(v,2);
            double H = q*fin_area*fin_height*C_H*delta;
            return H;
        }

        std::array<double,2> yaw_pitch_deflection(double k,double alpha,double beta,double cmd_ab_y,double cmd_ab_z){
            double yaw_defl = (cmd_ab_y*cos(beta) + cmd_ab_z*(sin(alpha)*sin(beta))) / k;
            double pitch_defl = (cmd_ab_z*cos(alpha)) / k;
            return {yaw_defl,pitch_defl};
        } 

        std::array<double,4> commanded_deflection_angles(std::array<double,2> y_p_defl){
            double invrs_sqrt2 = 1/sqrt(2);
            double defl_1 = y_p_defl[0]*invrs_sqrt2 + y_p_defl[1]*invrs_sqrt2;
            double defl_2 = y_p_defl[0]*-invrs_sqrt2 + y_p_defl[1]*invrs_sqrt2;
            double defl_3 = y_p_defl[0]*-invrs_sqrt2 + y_p_defl[1]*-invrs_sqrt2;
            double defl_4 = y_p_defl[0]*invrs_sqrt2 + y_p_defl[1]*-invrs_sqrt2;
            return {defl_1,defl_2,defl_3,defl_4};
        }

        std::array<double,4> update_deflection_angles(std::array<double,4> defl_angls, std::array<double,4> cmd_defl,double mz,double mvx,double mvy,double mvz){
            for (int i = 0; i < 4; i++){
                double delta_dot = (cmd_defl[i]-defl_angls[i])/act_tao;
                double tau_act = Kp*(cmd_defl[i]-defl_angls[i])-Kd*delta_dot;

                tau_act = std::max(-max_torque, std::min(tau_act,max_torque));

                double H = hinge_moment(mz,mvx,mvy,mvz,defl_angls[i]);
                double net_tau = tau_act - H;

                double delta_ddot = (net_tau-act_damping*delta_dot)/((1.0/3.0)*fin_mass*pow(fin_height,2));
                delta_dot += delta_ddot*delta_t;
                defl_angls[i] += delta_dot*delta_t;

                defl_angls[i] = cap_delfection(defl_angls[i],max_ang,-max_ang);
            }
            return defl_angls;
        }

    public:
        std::tuple<std::array<double,4>,double,double> CONTROL_SURFACE_DEFLECTIONS(double mz,double mvx,double mvy,double mvz,Quaterniond quat_vec,double cmd_ab_y,double cmd_ab_z){
            std::array<double,3> mv_b = missile_velocity_body_frame(mvx,mvy,mvz,quat_vec);
            double alpha = angle_of_attack(mv_b[0],mv_b[2]);
            double beta = side_slip(mv_b[0],mv_b[1]);
            double k = control_effectiveness(mz,mv_b[0],mv_b[1],mv_b[2]);
            std::array<double,2> y_p_defl = yaw_pitch_deflection(k,alpha,beta,cmd_ab_y,cmd_ab_z);
            std::array<double,4> cmd_defl = commanded_deflection_angles(y_p_defl);
            std::array<double,4> new_defl_angls = update_deflection_angles(defl_angls,cmd_defl,mz,mvx,mvy,mvz);
            return {new_defl_angls,alpha,beta};
        }
};

/**
 * The thrust vector control (TVC) class is currently not a part of this model but may be in the future.
 */
class TVC{
    
};

class Propulsion: public Enviornment{
    private:
        double mass_flow_rate(double fuel_density){
            double mass_flow_rate = fuel_density*cs_area*fuel_v;
            return mass_flow_rate;
        }

        double thrust_force(double fuel_density){
            double m_flow_rate = mass_flow_rate(fuel_density);
            double thrust = m_flow_rate*(exhaust_v);
            return thrust;
        }

        double get_fuel_mass(double m_flow_rate){
            double delta_m = m_flow_rate*delta_t;
            fuel_mass = fuel_mass - delta_m;
            return fuel_mass;
        }

        std::array<double,3> fuel_cm(){
            std::array<double,3> fuel_cm = {(body_length/2)*(fuel_mass/initial_fuel_mass),0,0};
            return fuel_cm;
        }

    public:
        std::tuple<double,double,std::array<double,3>> PROPULSION(double mz){
            double temp = temperature_at_altitude(mz);
            double fuel_density = fuel_density_at_altitude(temp);
            double m_flow_rate = mass_flow_rate(fuel_density);
            double thrust = thrust_force(fuel_density);
            fuel_mass = get_fuel_mass(m_flow_rate);
            f_cm = fuel_cm();
            return {thrust,fuel_mass,f_cm};
        }
};

class Aerodynamics: public Enviornment{
    public:
        double center_of_mass(std::array<double,3> fuel_cm){
            double cm = ((body_mass*body_cm)+(fuel_mass*fuel_cm[0])+(payload_mass*payload_cm)+4*(fin_mass*front_fin_d)+(nose_mass*(3/8)*body_radius))/total_mass;
            return cm;
        }

        std::tuple<double,double> q_M(Vector3d mp, Vector3d mv){
            double v = sqrt(pow(mv.x(),2)+pow(mv.y(),2)+pow(mv.z(),2));
            double air_temp = temperature_at_altitude(mp.z());
            double air_pres = air_pressure_at_altitude(air_pressure_sea_lvl,mp.z());
            double air_dens = air_density_at_altitude(air_pres,air_temp);
            double q = .5 * air_dens * pow(v,2);
            double M = v/speed_of_sound;
            return {q,M};
        }

        std::tuple<double,double,double> aero_wind_to_body(Vector3d aero_vec,double alpha,double beta){
            Matrix3d C_wb;
            C_wb << cos(alpha)*cos(beta), sin(beta), sin(alpha)*cos(beta),
                    -cos(alpha)*sin(beta), cos(beta), -sin(alpha)*sin(beta),
                    -sin(alpha), 0, cos(alpha);
            
            Matrix3d C_bw = C_wb.transpose();
            Vector3d trans_vec = C_bw*aero_vec;

            return {trans_vec.x(),trans_vec.y(),trans_vec.z()};
        }

        std::tuple<double,double,double> control_deflections(std::array<double,4> defl_angls){
            double inv_2sqrt2 = 1.0/(2.0*sqrt(2.0));

            double delta_y = (defl_angls[0]-defl_angls[1]-defl_angls[2]+defl_angls[3])*inv_2sqrt2;
            double delta_p = (defl_angls[0]+defl_angls[1]-defl_angls[2]-defl_angls[3])*inv_2sqrt2;
            double delta_r = (defl_angls[0]-defl_angls[1]+defl_angls[2]-defl_angls[3])*inv_2sqrt2;

            return {delta_y,delta_p,delta_r};
        }

        std::tuple<double,double,double> control_derivatives(double Mach, Vector3d body_rot_vec){
            double p_hat = body_rot_vec.x() * body_radius / (Mach*speed_of_sound);
            double q_hat = body_rot_vec.y() * body_length / (2.0*Mach*speed_of_sound);
            double r_hat = body_rot_vec.z() * body_length / (2.0*Mach*speed_of_sound);

            return {p_hat,q_hat,r_hat};
        }

        double aero_force_X(double q,double cross_area){
            double X = q*cross_area*CD;
            return X;
        }

        double aero_force_Y(double q,double cross_area){
            double Y = q*cross_area*CC;
            return Y;
        }

        double aero_force_Z(double q,double cross_area){
            double Z = -q*cross_area*CL;
            return Z;
        }

        double aero_moment_L(double q,double cross_area){
            double L = q*cross_area*body_radius*Cl;
            return L;
        }

        double aero_moment_M(double q,double cross_area){
            double M = q*cross_area*body_length*Cm;
            return M;
        }

        double aero_moment_N(double q,double cross_area){
            double N = q*cross_area*body_radius*Cn;
            return N;
        }

        std::tuple<Vector3d,Vector3d,std::string> AERODYNAMICS(std::array<double,3> fuel_cm,Vector3d m_vel_arr,Vector3d m_pos_arr,double alpha,double beta,
                                                                            double thrust,AeroCoefficientTables& table,std::array<double,4> defl_angls,Vector3d body_rot_vec){
            double cm = center_of_mass(fuel_cm);
            double q,Mach;
            double X,Y,Z,L,M,N;
            double delta_y,delta_p,delta_r;
            double p_hat,q_hat,r_hat;

            std::tie(q,Mach) = q_M(m_pos_arr,m_vel_arr);

            std::tie(delta_y,delta_p,delta_r) = control_deflections(defl_angls);
            std::tie(CD,CC,CL,Cl,Cm,Cn) = table.interpolate_coefficients(Mach,alpha,beta,delta_p,delta_y,delta_r);

            std::tie(p_hat,q_hat,r_hat) = control_derivatives(Mach,body_rot_vec);
            Cl += -5.0 * p_hat;
            Cm += -25.0 * q_hat;
            Cn += -14.0 * r_hat;
            
            double X_w = aero_force_X(q,cross_area);
            double Y_w = aero_force_Y(q,cross_area);
            double Z_w = aero_force_Z(q,cross_area);
            double L_w = aero_moment_L(q,cross_area);
            double M_w = aero_moment_M(q,cross_area);
            double N_w = aero_moment_N(q,cross_area);

            std::tie(X,Y,Z) = aero_wind_to_body({-X_w,-Y_w,-Z_w},alpha,beta);
            std::tie(L,M,N) = aero_wind_to_body({L_w,M_w,N_w},alpha,beta);

            Vector3d aero_forces = Vector3d(X+thrust,Y,Z);
            Vector3d aero_moments = Vector3d(L,M,N);

            std::string data = std::to_string(CD)+","+std::to_string(CC)+","+std::to_string(CL)+","
                               +std::to_string(Cl)+","+std::to_string(Cm)+","+std::to_string(Cn)+","
                               +std::to_string(X+thrust)+","+std::to_string(Y)+","+std::to_string(Z)+","
                               +std::to_string(L)+","+std::to_string(M)+","+std::to_string(N);

            return {aero_forces,aero_moments,data};
        }
};

class Dynamics: public Enviornment{
    public:
        double center_of_mass(std::array<double,3> fuel_cm){
            double cm = ((body_mass*body_cm)+(fuel_mass*fuel_cm[0])+(payload_mass*payload_cm)+4*(fin_mass*front_fin_d)+(nose_mass*(3/8)*body_radius))/total_mass;
            return cm;
        }
        
        double I_xx(){
            double I_body = .5*body_mass*pow(body_radius,2);
            double I_nose = .3*nose_mass*pow(body_radius,2);
            double I_f = (fin_mass/36.0)*(pow(fin_base,2)+pow(fin_height,2));
            double I_fin = 8*(I_f + fin_mass*pow(body_radius+(fin_base/3),2));
            double I_fuel = .5*fuel_mass*pow(fuel_tank_radius,2);
            double I_xx = I_body+I_nose+I_fin+I_fuel;
            return I_xx;
        }

        double I_yy_zz(double cm){
            double I_b = (1.0/12.0)*body_mass*(3*pow(body_radius,2)+pow(body_length,2));
            double I_body = I_b + body_mass*pow(body_cm-cm,2);
            double I_n = (3.0/20.0)*nose_mass*pow(body_radius,2);
            double I_nose = I_n + nose_mass*pow((body_length-cm+(3.0/8.0)*body_radius),2);
            double I_payload = payload_mass*pow(payload_cm-cm,2);
            double I_f = (fin_mass/36.0)*pow(fin_height,2);
            double I_fin_rear = 4*(I_f + fin_mass*pow((fin_base/3)-cm,2));
            double I_fin_front = 4*(I_f + fin_mass*pow((front_fin_d+(fin_base/3))-cm,2));
            double I_fl = (1.0/12.0)*fuel_mass*(3*pow(fuel_tank_radius,2)+pow(fuel_tank_length,2));
            double I_fuel = I_fl + fuel_mass*pow((fuel_tank_length/2)-cm,2);
            double Iyy = I_body+I_nose+I_payload+I_fin_front+I_fin_rear+I_fuel;
            return Iyy;
        }

        Vector3d velocity_transform_f_b(Vector3d v_vec,Quaterniond quat_vec){
            Matrix3d Rfb = transformation_Matrix_Fixed_to_Body(quat_vec);
            Vector3d v_vec_tranf = Rfb*v_vec;
            return v_vec_tranf;
        }

        Vector3d velocity_transform_b_f(Vector3d v_vec,Quaterniond quat_vec){
            Matrix3d Rbf = transformation_Matrix_Body_to_Fixed(quat_vec);
            Vector3d v_vec_tranf = Rbf*v_vec;
            return v_vec_tranf;
        }

        Vector3d change_in_body_rates(Vector3d body_rate,Vector3d aero_moments,double I_xx,double I_yy,double I_zz){
            Matrix3d I;
            I << I_xx, 0, 0,
                 0, I_yy, 0,
                 0, 0, I_zz;
            Vector3d dt_body_rates = I.inverse()*(aero_moments - body_rate.cross(I*body_rate));
            return dt_body_rates;
        }

        Vector3d change_in_vel_body(Vector3d body_rate,Vector3d vel_vec_b,Vector3d aero_forces,Quaterniond quat_vec){
            Vector3d g_inertial(0.0, 0.0, -9.81);
            Matrix3d Rfb = transformation_Matrix_Fixed_to_Body(quat_vec);
            Vector3d dt_vel_b = -body_rate.cross(vel_vec_b) + aero_forces/total_mass + Rfb*g_inertial;
            return dt_vel_b;
        }

        Vector3d update_vel_body(Vector3d vel_vec_b,Vector3d dt_v_b){
            Vector3d v_b = vel_vec_b + dt_v_b*delta_t;
            return v_b;
        }

        Vector3d update_body_rate(Vector3d body_rate,Vector3d dt_body_rate){
            Vector3d omega = body_rate + dt_body_rate*delta_t;
            return omega;
        }

        Vector3d update_pos_intertial(Vector3d pos_vec,Vector3d vel_vec){
            Vector3d m_pos = pos_vec + vel_vec*delta_t;
            return m_pos;
        }
    
        std::tuple<Vector3d,Vector3d,Vector3d> DYNAMICS(std::array<double,3> fuel_cm,Vector3d v_vec,
                Vector3d aero_forces,Vector3d aero_moments,Vector3d pos_vec,Quaterniond quat_vec,Vector3d body_rot_vec){
            double cm = center_of_mass(fuel_cm);
            double Ixx = I_xx();
            double Iyy = I_yy_zz(cm);
            double Izz = Iyy;
            Vector3d v_vec_b;
            v_vec_b = velocity_transform_f_b(v_vec,quat_vec);
            Vector3d dt_v_b = change_in_vel_body(body_rot_vec,v_vec_b,aero_forces,quat_vec);
            Vector3d dt_body_rate = change_in_body_rates(body_rot_vec,aero_moments,Ixx,Iyy,Izz);
            Vector3d new_v_b = update_vel_body(v_vec_b,dt_v_b);
            Vector3d new_body_rate = update_body_rate(body_rot_vec,dt_body_rate);
            Vector3d new_v_vec = velocity_transform_b_f(new_v_b,quat_vec);
            Vector3d new_pos_vec = update_pos_intertial(pos_vec,new_v_vec);
            return {new_pos_vec,new_v_vec,new_body_rate};
        }
};

class InertialSensor: public Enviornment{
    public:
        std::tuple<Vector3d,Vector3d,Vector3d,std::string> KalmanFilter(Vector3d pos_vec,Vector3d vel_vec,
                Vector3d acc_vec,EKF& ekf,int step,std::mt19937& rng,std::normal_distribution<double>& nz,Vector3d body_rot_vec){
            
            Matrix3d Rpos = Matrix3d::Identity()*(0.05*0.05); 
            Matrix3d Ratt = Matrix3d::Identity()*(0.02*0.02);

            Vector3d true_p = pos_vec;
            Vector3d true_v = vel_vec;

            Quaterniond true_q = quat_vec;

            Vector3d a_fixed = acc_vec;
            Vector3d a_body = true_q.conjugate().toRotationMatrix()*(a_fixed-ekf.g);
            Vector3d gryo_body = body_rot_vec;

            auto noisyVec = [&](const Vector3d &v, double sigma) {
                Vector3d r;
                r.x() = v.x() + sigma * nz(rng);
                r.y() = v.y() + sigma * nz(rng);
                r.z() = v.z() + sigma * nz(rng);
                return r;
            };

            IMU imu_meas;
            imu_meas.acc = noisyVec(a_body, ekf.sigma_acc);
            imu_meas.gyro = noisyVec(gryo_body,ekf.sigma_gyro);

            ekf.predictInertial(imu_meas,delta_t);

            Vector3d meas_p = Vector3d::Zero();
            bool has_meas = false;
            if (step %500 == 0) {
                Vector3d zpos = true_p + noisyVec(Vector3d::Zero(), std::sqrt(Rpos(0,0)));
                ekf.updateInertialPos(zpos, Rpos);
                meas_p = zpos;
                has_meas = true;
                } else {
                    has_meas = false;
            }

            if (step % 100 == 0) {
                Quaterniond zq = true_q.conjugate();
                Vector3d eps;
                eps.x() = 0.02 * nz(rng);
                eps.y() = 0.02 * nz(rng);
                eps.z() = 0.02 * nz(rng);
                double en = eps.norm();
                Quaterniond qnoise;
                if (en > 1e-8) qnoise = AngleAxisd(en, eps/en);
                else { qnoise.w()=1.0; qnoise.vec() = 0.5*eps; qnoise.normalize(); }
                zq = qnoise * zq;
                ekf.updateInertialAtt(zq, Ratt);
            }

            double meas_px = has_meas ? meas_p.x() : std::numeric_limits<double>::quiet_NaN();
            std::string pos_data = std::to_string(true_p.x())+","+std::to_string(ekf.pos.x())+","+std::to_string(meas_px);
            std::string vel_data = ","+std::to_string(true_v.x())+","+std::to_string(ekf.vel.x());
            std::string acc_data = ","+std::to_string(imu_meas.acc.x());

            return {ekf.pos,ekf.vel,imu_meas.acc,pos_data+vel_data+acc_data};
        }
};

class TerminalSensor : public Enviornment{
    public:
        std::tuple<Vector3d,Vector3d,Vector3d,std::string>
            KalmanFilter(Vector3d est_m_pos, Vector3d est_m_vel,
                        Vector3d t_pos_vec, Vector3d t_vel_vec,
                        Vector3d t_acc_vec, EKF& tekf, int step,
                        std::normal_distribution<double>& nz, std::mt19937& rng)
            {
                Vector3d p_m = est_m_pos;
                Vector3d v_m = est_m_vel;

                Vector3d true_p = t_pos_vec;
                Vector3d true_v = t_vel_vec;
                Vector3d true_a = t_acc_vec;

                tekf.predictTerminal(delta_t);

                int meas_interval = 50;
                if (step % meas_interval == 0) {

                    Vector3d rel_pos = (true_p - p_m) + 5.0 * Vector3d(nz(rng), nz(rng), nz(rng)); 
                    Vector3d rel_vel = (true_v - v_m) + 1.0 * Vector3d(nz(rng), nz(rng), nz(rng)); 

                    VectorXd z(6);
                    z << rel_pos, rel_vel;

                    MatrixXd R = MatrixXd::Zero(6,6);
                    R.block<3,3>(0,0) = 5.0 * Matrix3d::Identity();   
                    R.block<3,3>(3,3) = 1.0  * Matrix3d::Identity();  

                    tekf.updateTerminalPosVel(z, R, p_m, v_m);
                }

                std::string pos_data = std::to_string(true_p.x()) + "," + std::to_string(tekf.getTerminalPos().x());
                std::string vel_data = "," + std::to_string(true_v.x()) + "," + std::to_string(tekf.getTerminalVel().x());
                std::string acc_data = "," + std::to_string(true_a.x()) + "," + std::to_string(tekf.getTerminalAcc().x());

                return {tekf.getTerminalPos(),tekf.getTerminalVel(),tekf.getTerminalAcc(), pos_data + vel_data + acc_data };
            }

};

class FlightController: public Enviornment{
    private:
        Matrix4d Quaternion_Time_Evolution_Matrix(Vector3d body_rot_vec){
            Matrix4d Omega_matrix;
            Omega_matrix << 0, -(body_rot_vec.x()/2), -(body_rot_vec.y()/2), -(body_rot_vec.z()/2),
                            (body_rot_vec.x()/2), 0, (body_rot_vec.z()/2), -(body_rot_vec.y()/2),
                            (body_rot_vec.y()/2), -(body_rot_vec.z()/2), 0, (body_rot_vec.x()/2),
                            (body_rot_vec.z()/2), (body_rot_vec.y()/2), -(body_rot_vec.x()/2), 0;
            return Omega_matrix;
        }

        Vector3d Executed_Acceleration_Fixed_Frame(Vector3d executed_body_a_vec,Matrix3d Rbf){
            Vector3d acc_vec = Rbf*executed_body_a_vec;
            return acc_vec;
        }

        Vector3d Body_Rotation_Vector(Vector3d mv,Vector3d body_a_vec){
            Vector3d body_rot_vec = (mv.cross(body_a_vec))/(mv.dot(mv));
            return body_rot_vec;
        }

        std::tuple<Vector4d,Vector4d> Quaternion_Vector_Change_Rate(Matrix4d omega_matrix,Quaterniond quat_vec){
            Vector4d q = Vector4d(quat_vec.w(),quat_vec.x(),quat_vec.y(),quat_vec.z());
            Vector4d quat_vec_rate = omega_matrix*q;
            return {quat_vec_rate,q};
        }

        std::array<double,3> Euler_Angles(Quaterniond& q){
            double qw = q.w(), qx = q.x(), qy = q.y(), qz = q.z();

            double sinr_cosp = 2.0 * (qw*qx + qy*qz);
            double cosr_cosp = 1.0 - 2.0 * (qx*qx + qy*qy);
            double roll = std::atan2(sinr_cosp, cosr_cosp);

            double sinp = 2.0 * (qw*qy - qz*qx);
            if (sinp > 1.0)  sinp = 1.0;
            if (sinp < -1.0) sinp = -1.0;
            double pitch = std::asin(sinp);

            double siny_cosp = 2.0 * (qw*qz + qx*qy);
            double cosy_cosp = 1.0 - 2.0 * (qy*qy + qz*qz);
            double yaw = std::atan2(siny_cosp, cosy_cosp);

            return {roll,pitch,yaw};
        }

    public:
        std::tuple<Vector3d,std::array<double,3>,Vector3d> FLIGHT_CONTROL(Vector3d executed_body_a_vec,Quaterniond quat_vec,Vector3d missile_vel_vec,double delta_t,Vector3d body_rot_vec,Fidelity fid_){
            Matrix3d Rbf = transformation_Matrix_Body_to_Fixed(quat_vec);
            acc_vec = Executed_Acceleration_Fixed_Frame(executed_body_a_vec,Rbf);
            if (fid_ == TRIM){
                body_rot_vec = Body_Rotation_Vector(missile_vel_vec,executed_body_a_vec);
            }
            Matrix4d omega_matrix = Quaternion_Time_Evolution_Matrix(body_rot_vec);
            Vector4d quat_vec_change_rate,q;
            std::tie(quat_vec_change_rate,q) = Quaternion_Vector_Change_Rate(omega_matrix,quat_vec);
            q = q + quat_vec_change_rate*delta_t;
            quat_vec = Quaterniond(q(0), q(1), q(2), q(3));
            std::array<double,3> euler_angles = Euler_Angles(quat_vec);
            return {acc_vec,euler_angles,body_rot_vec};
        }
};

class NavigationProcessor: public Enviornment{
    private:
        Vector3d Relative_Position_Vector(Vector3d missile_pos_vec,Vector3d target_pos_vec){
            rel_pos_vec = missile_pos_vec - target_pos_vec;
            return rel_pos_vec;
        }

        Vector3d Relative_Velocity_Vector(Vector3d missile_vel_vec,Vector3d target_vel_vec){
            rel_vel_vec = missile_vel_vec - target_vel_vec;
            return rel_vel_vec;
        }

        Vector3d Relative_Acceleration_Vector(Vector3d missile_acc_vec,Vector3d target_acc_vec){
            rel_acc_vec = missile_acc_vec - target_acc_vec;
            return rel_acc_vec;
        }

        double Range(Vector3d mpos_vec,Vector3d tpos_vec){
            double R = sqrt(pow(mpos_vec.x()-tpos_vec.x(),2)+pow(mpos_vec.y()-tpos_vec.y(),2)+pow(mpos_vec.z()-tpos_vec.z(),2));
            return R;
        }

        Vector3d Sight_Line_Rotation_Vector(Vector3d pos_vec,Vector3d vel_vec){
            Vector3d sight_line_rot_vec = (pos_vec.cross(vel_vec))/(pos_vec.dot(pos_vec));
            return sight_line_rot_vec;
        }

    public:
        std::tuple<Vector3d,Vector3d,Vector3d> TRANSLATIONAL_KINEMATICS(Vector3d missile_pos_vec,Vector3d missile_vel_vec,Vector3d missile_acc_vec,
                                                                        Vector3d target_pos_vec,Vector3d target_vel_vec,Vector3d target_acc_vec){
            Vector3d rel_pos_vec = Relative_Position_Vector(missile_pos_vec,target_pos_vec);
            Vector3d rel_vel_vec = Relative_Velocity_Vector(missile_vel_vec,target_vel_vec);
            Vector3d rel_acc_vec = Relative_Acceleration_Vector(missile_acc_vec,target_acc_vec);
            return {rel_pos_vec,rel_vel_vec,rel_acc_vec};
        }

        std::tuple<double,Vector3d> SEEKER_ROTATIONAL_KINEMATICS(Vector3d rel_pos_vec,Vector3d rel_vel_vec,Vector3d mpos_vec,Vector3d tpos_vec){
            R = Range(mpos_vec,tpos_vec);
            Vector3d sight_line_rot_vec = Sight_Line_Rotation_Vector(rel_pos_vec,rel_vel_vec);
            return {R,sight_line_rot_vec};
        }
};

class GuidanceProcessor: public Enviornment{
    private:
        Matrix3d Navigational_Constant_Matrix(){
            Matrix3d N_matrix;
            N_matrix << N1_i, 0, 0,
                        0, N2_i, 0,
                        0, 0, N3_i;
            return N_matrix;
        }

        Quaterniond Quaternion_Vector(double roll, double pitch, double yaw){
            double q1_i = cos(roll/2)*cos(pitch/2)*cos(yaw/2)+sin(roll/2)*sin(pitch/2)*sin(yaw/2);
            double q2_i = sin(roll/2)*cos(pitch/2)*cos(yaw/2)-cos(roll/2)*sin(pitch/2)*sin(yaw/2);
            double q3_i = cos(roll/2)*sin(pitch/2)*cos(yaw/2)+sin(roll/2)*cos(pitch/2)*sin(yaw/2);
            double q4_i = cos(roll/2)*cos(pitch/2)*sin(yaw/2)-sin(roll/2)*sin(pitch/2)*cos(yaw/2);
            quat_vec = Quaterniond(q1_i,q2_i,q3_i,q4_i);
            return quat_vec;
        } 

        Vector3d Demanded_Body_Rotaion_Vector(Matrix3d N_matrix,Vector3d sight_line_rot_vec){
            Vector3d dem_rot_vec = N_matrix*sight_line_rot_vec;
            return dem_rot_vec;
        }

        Vector3d Demanded_Acceleration_Vector(Vector3d dem_body_rot_vec,Vector3d mv){
            Vector3d dem_acc_vec = dem_body_rot_vec.cross(mv);
            return dem_acc_vec;
        }

        Vector3d Demanded_Acceleration_Vector_Body_Frame(Vector3d dem_acc_vec, Matrix3d Rfb){
            Vector3d dem_acc_vec_body = Rfb*dem_acc_vec;
            return dem_acc_vec_body;
        }

        Vector3d Target_Acceleration_Influence(Matrix3d Rfb,Matrix3d N_matrix,Vector3d ta){
            Vector3d a_j = N_matrix*ta;
            Vector3d t_a_infl = Rfb*a_j;
            return t_a_infl;
        }

    public:
        std::tuple<Vector3d,Matrix3d,Quaterniond> APN_GUIDANCE_MODULE(Vector3d sight_line_rot_vec,double N2_i,double N3_i,Vector3d missile_vel_vec,std::array<double,3> euler_angles,Vector3d target_acc_vec){
            Matrix3d N_matrix = Navigational_Constant_Matrix();
            Vector3d dem_body_rot_vec = Demanded_Body_Rotaion_Vector(N_matrix,sight_line_rot_vec);
            Vector3d dem_acc_vec = Demanded_Acceleration_Vector(dem_body_rot_vec,missile_vel_vec);
            Quaterniond quat_vec = Quaternion_Vector(euler_angles[0],euler_angles[1],euler_angles[2]);
            Matrix3d Rfb = transformation_Matrix_Fixed_to_Body(quat_vec);
            Vector3d dem_acc_vec_body = Demanded_Acceleration_Vector_Body_Frame(dem_acc_vec,Rfb);
            Vector3d t_a_infl = Target_Acceleration_Influence(Rfb,N_matrix,target_acc_vec);
            dem_acc_vec_body = dem_acc_vec_body + t_a_infl;
            return {dem_acc_vec_body,Rfb,quat_vec};
        }
};

class Autopilot: public Enviornment{
    private:
        Matrix3d Autopilot_Time_Constant_Matrix(){
            Matrix3d Delta;
            Delta << taox_i, 0, 0,
                     0, taoy_i, 0,
                     0, 0, taoz_i;
            return Delta;
        }

        Vector3d Body_Acceleration_Vector(Vector3d ma,Matrix3d Rfb){
            Vector3d body_a_vec = Rfb*ma;
            return body_a_vec;
        }

        Vector3d Body_Acceleration_Change_Rate(Matrix3d auto_time_const_matrix,Vector3d body_a_vec,Vector3d dem_a_acc_vec){
            Vector3d body_acc_chng_rate = -auto_time_const_matrix*body_a_vec + auto_time_const_matrix*dem_a_acc_vec;
            return body_acc_chng_rate;
        }

    public:
        Vector3d AUTOPILOT_RESPONSE(double taox_i,double taoy_i,double taoz_i,Matrix3d Rfb,Vector3d dem_acc_vec_body,Vector3d missile_acc_vec,double delta_t){
            Matrix3d autopilot_matrix = Autopilot_Time_Constant_Matrix();
            Vector3d body_a_vec = Body_Acceleration_Vector(missile_acc_vec,Rfb);
            Vector3d body_a_vec_change_rate = Body_Acceleration_Change_Rate(autopilot_matrix,body_a_vec,dem_acc_vec_body);
            Vector3d executed_body_a_vec = body_a_vec + body_a_vec_change_rate;
            return executed_body_a_vec;
        }
};

class Airframe: public Enviornment{
    public:
        Fin_Actuator fin_actuator;
        TVC tvc;
        Propulsion propulsion;
        Dynamics dynamics;
        Aerodynamics aerodynamics;

        std::tuple<double,double,std::array<double,3>> JET_ENGINE(double mz){
            std::tie(thrust,fuel_mass,f_cm) = propulsion.PROPULSION(mz);
            return {thrust,fuel_mass,f_cm};
        }

        std::tuple<std::array<double,4>,double,double> CONTROL_SURFACES(double mz,double mvx,double mvy,double mvz,Quaterniond quat_vec,double cmd_ab_y,double cmd_ab_z){
            std::tie(defl_angls,alpha,beta) = fin_actuator.CONTROL_SURFACE_DEFLECTIONS(mz,mvx,mvy,mvz,quat_vec,cmd_ab_y,cmd_ab_z);
            return {defl_angls,alpha,beta};
        }

        std::tuple<Vector3d,Vector3d,std::string> AERO(Vector3d pos_vec,Vector3d vel_vec,double thrust,Vector3d body_rot_vec,AeroCoefficientTables& table,std::array<double,4> defl_angls){
            std::string aerodata;
            std::tie(aero_forces,aero_moments,aerodata) = aerodynamics.AERODYNAMICS(f_cm,vel_vec,pos_vec,alpha,beta,thrust,table,defl_angls,body_rot_vec);
            return {aero_forces,aero_moments,aerodata};
        }

        std::tuple<Vector3d,Vector3d,Vector3d> DYN(Vector3d v_vec,Vector3d p_vec,Quaterniond quat_vec,Vector3d body_rot_vec){
            std::tie(new_p,new_omega,new_v) = dynamics.DYNAMICS(f_cm,v_vec,aero_forces,aero_moments,p_vec,quat_vec,body_rot_vec);
            return {new_p,new_v,new_omega};
        }
};

class Target{
    public:
        double t_x;
        double t_y;
        double t_z;
        double t_vx;
        double t_vy;
        double t_vz;
        double t_ax;
        double t_ay;
        double t_az;

        Target(double t_x, double t_y,double t_z,double t_vx,double t_vy,double t_vz,double t_ax,double t_ay,double t_az){
            this->t_x = t_x;
            this->t_y = t_y;
            this->t_z = t_z;
            this->t_vx = t_vx;
            this->t_vy = t_vy;
            this->t_vz = t_vz;
            this->t_ax = t_ax;
            this->t_ay = t_ay;
            this->t_az = t_az;
        }

    std::tuple<double,double,double,double,double,double> update_params(double delta_t){
        double new_vx = t_vx + t_ax*delta_t;
        double new_vy = t_vy + t_ay*delta_t;
        double new_vz = t_vz + t_az*delta_t;
        double new_x = t_x + new_vx*delta_t;
        double new_y = t_y + new_vy*delta_t;
        double new_z = t_z + new_vz*delta_t;
        return {new_vx,new_vy,new_vz,new_x,new_y,new_z};
    }
};

class Missile: public Enviornment{
    private:
        Airframe airframe;
        InertialSensor inertialsensor;
        TerminalSensor terminalsensor;
        FlightController flightcontroller;
        NavigationProcessor navigationprocessor;
        GuidanceProcessor guidanceprocessor;
        Autopilot autopilot;

    public:
        double m_x,m_y,m_z;
        double m_vx,m_vy,m_vz;
        double m_ax,m_ay,m_az;
        double roll,pitch,yaw;
        Vector3d body_rot_vec;

        Missile(double m_x, double m_y,double m_z,double m_vx,double m_vy,double m_vz,double m_ax,double m_ay,double m_az,double roll,double pitch,double yaw){
            this->m_x = m_x;
            this->m_y = m_y;
            this->m_z = m_z;
            this->m_vx = m_vx;
            this->m_vy = m_vy;
            this->m_vz = m_vz;
            this->m_ax = m_ax;
            this->m_ay = m_ay;
            this->m_az = m_az;
            this->roll = roll;
            this->pitch = pitch;
            this->yaw = yaw;
            Vector3d body_rot_vec = Vector3d(0.0,0.0,0.0);
        }

        std::tuple<double,double,double,double,double> COLLISION_COURSE_GEOMETRY(VectorXd mv,VectorXd tv,VectorXd rel_pos_vec,VectorXd rel_vel_vec,double R){
            double psi_t = atan2(tv(1),tv(0));
            double theta_t = atan2(tv(2),sqrt(tv(0)*tv(0)+tv(1)*tv(1)));

            double psi_tm = atan2(rel_pos_vec(1),rel_pos_vec(0));
            double theta_tm = atan2(rel_pos_vec(2),sqrt(rel_pos_vec(0)*rel_pos_vec(0)+rel_pos_vec(1)*rel_pos_vec(1)));

            double beta_tm = acos(cos(theta_t)*cos(psi_t)*cos(theta_tm)*cos(psi_tm)+cos(theta_t)*sin(psi_t)*cos(theta_tm)*sin(psi_tm)+sin(theta_t)*sin(theta_tm));
            double beta_cc_mt = asin((tv.norm()/mv.norm())*sin(beta_tm));

            double vc_cc_mt = mv.norm()*cos(beta_cc_mt)-tv.norm()*cos(beta_tm);
            double T_go = R/vc_cc_mt;

            double theta_cc_mt = asin((vc_cc_mt/mv.norm())*sin(theta_tm)+(tv.norm()/mv.norm())*sin(theta_t));
            double psi_cc_mt = atan2(((vc_cc_mt/mv.norm())*cos(theta_tm)*sin(psi_tm)+(tv.norm()/mv.norm())*cos(theta_t)*sin(psi_t)),((vc_cc_mt/mv.norm())*cos(theta_tm)*cos(psi_tm)+(tv.norm()/mv.norm())*cos(theta_t)*cos(psi_t)));

            double psi_los_dot = (rel_pos_vec(0)*rel_vel_vec(1)-rel_pos_vec(1)*rel_vel_vec(0))/(rel_pos_vec(0)*rel_pos_vec(0)+rel_pos_vec(1)*rel_pos_vec(1));
            double theta_los_dot = ((rel_pos_vec(0)*rel_pos_vec(0)+rel_pos_vec(1)*rel_pos_vec(1))*rel_vel_vec(2)-rel_pos_vec(2)*(rel_pos_vec(0)*rel_vel_vec(0)+rel_pos_vec(1)*rel_vel_vec(1)))/(sqrt(rel_pos_vec(0)*rel_pos_vec(0)+rel_pos_vec(1)*rel_pos_vec(1))*(rel_pos_vec(0)*rel_pos_vec(0)+rel_pos_vec(1)*rel_pos_vec(1)+rel_pos_vec(2)*rel_pos_vec(2)));
            
            return {T_go,theta_cc_mt*(180/3.14159),psi_cc_mt*(180/3.14159),theta_los_dot,psi_los_dot};
        }

        std::tuple<double,double,double,double,double,double,double,double,double> update_params(Vector3d acc_vec,double delta_t){
            double new_ax = acc_vec.x();
            double new_ay = acc_vec.y();
            double new_az = acc_vec.z();
            double new_vx = m_vx + new_ax*delta_t;
            double new_vy = m_vy + new_ay*delta_t;
            double new_vz = m_vz + new_az*delta_t;
            double new_x = m_x + (new_vx*delta_t) + (new_vx*pow(delta_t,2));
            double new_y = m_y + (new_vy*delta_t) + (new_vy*pow(delta_t,2));
            double new_z = m_z + (new_vz*delta_t) + (new_vz*pow(delta_t,2));
            return {new_ax,new_ay,new_az,new_vx,new_vy,new_vz,new_x,new_y,new_z};
        }

        std::tuple<Vector3d,Vector3d,Vector3d,Vector3d,Vector3d,Vector3d,std::string,std::string> SENSORS(Vector3d mpos_vec,Vector3d mvel_vec,
                Vector3d macc_vec,Vector3d tpos_vec,Vector3d tvel_vec,Vector3d tacc_vec,
                EKF& ekf,int step,std::mt19937& rng,std::normal_distribution<double>& nz,EKF& tekf,Vector3d body_rot_vec){

            std::string ekfdata;
            std::string tekfdata;
            std::tie(est_m_pos,est_m_vel,est_m_acc,ekfdata) = inertialsensor.KalmanFilter(mpos_vec,mvel_vec,macc_vec,ekf,step,rng,nz,body_rot_vec);
            std::tie(est_t_pos,est_t_vel,est_t_acc,tekfdata) = terminalsensor.KalmanFilter(est_m_pos,est_m_vel,tpos_vec,tvel_vec,tacc_vec,tekf,step,nz,rng);
            return {est_m_pos,est_m_vel,est_m_acc,est_t_pos,est_t_vel,est_t_acc,ekfdata,tekfdata};
        }

        std::tuple<double,Vector3d,std::array<double,3>,Quaterniond,Vector3d,Vector3d> INTELLIGENT_INTEGRATED_GUIDANCE_NAVIGATION_AND_CONTROL_SYSTEM(Vector3d missile_pos_vec,Vector3d missile_vel_vec,Vector3d missile_acc_vec,Vector3d target_pos_vec,Vector3d target_vel_vec,Vector3d target_acc_vec,
                        Vector3d mtrue_p_vec,Vector3d mtrue_v_vec,Vector3d mtrue_a_vec,Vector3d ttrue_p_vec,Vector3d body_rot_vec,Fidelity fid_){
            std::tie(rel_pos_vec,rel_vel_vec,rel_acc_vec) = navigationprocessor.TRANSLATIONAL_KINEMATICS(missile_pos_vec,missile_vel_vec,missile_acc_vec,target_pos_vec,target_vel_vec,target_acc_vec);
            std::tie(R,sight_line_rot_vec) = navigationprocessor.SEEKER_ROTATIONAL_KINEMATICS(rel_pos_vec,rel_vel_vec,mtrue_p_vec,ttrue_p_vec);
            std::tie(dem_acc_vec_body,Rfb,quat_vec) = guidanceprocessor.APN_GUIDANCE_MODULE(sight_line_rot_vec,N2_i,N3_i,missile_vel_vec,euler_angles,target_acc_vec);
            executed_body_a_vec = autopilot.AUTOPILOT_RESPONSE(taox_i,taoy_i,taoz_i,Rfb,dem_acc_vec_body,mtrue_a_vec,delta_t);
            std::tie(acc_vec,euler_angles,body_rot_vec) = flightcontroller.FLIGHT_CONTROL(executed_body_a_vec,quat_vec,mtrue_v_vec,delta_t,body_rot_vec,fid_);
            return {R,acc_vec,euler_angles,quat_vec,executed_body_a_vec,body_rot_vec};
        }

        std::tuple<double,double,std::array<double,4>,Vector3d,Vector3d,Vector3d,Vector3d,Vector3d,std::string> AIRFRAME(Vector3d missile_pos_vec,Vector3d missile_vel_vec,Vector3d executed_body_a_vec,Quaterniond quat_vec,Vector3d body_rot_vec,AeroCoefficientTables& table){
            std::string aerodata;
            std::tie(thrust,fuel_mass,f_cm) = airframe.JET_ENGINE(missile_pos_vec[2]);
            std::tie(defl_angls,alpha,beta) = airframe.CONTROL_SURFACES(missile_pos_vec[2],missile_vel_vec[0],missile_vel_vec[1],missile_vel_vec[2],quat_vec,executed_body_a_vec[1],executed_body_a_vec[2]);
            std::tie(aero_forces,aero_moments,aerodata) = airframe.AERO(missile_pos_vec,missile_vel_vec,thrust,body_rot_vec,table,defl_angls);
            std::tie(new_p,new_omega,new_v) = airframe.DYN(missile_vel_vec,missile_pos_vec,quat_vec,body_rot_vec);
            return {thrust,fuel_mass,defl_angls,aero_forces,aero_moments,new_p,new_v,new_omega,aerodata};
        }
};

class Simulation: public Enviornment{
    public:
        std::mt19937 rng;
        std::normal_distribution<double> nz;

        Simulation() : rng(std::random_device{}()), nz(0.0, 1.0) {}

        void run(const std::string& cfgfile){
            std::cout<<"Initializing Simulation\n"<<std::endl;
            readConfig(cfgfile);
            std::cout<<"Reading Config File\n"<<std::endl;
            setState();
            std::cout<<"State set\n"<<std::endl;
    
            CSV_Writer csv("C:\\Software Development\\6 DOF Sim\\Data\\GNC_Data2.csv");
            csv.write_line("Time,Range,Mx,My,Mz,Tx,Ty,Tz,q1,q2,q3,q4,q,Thrust,fuel mass,delta1,delta2,delta3,delta4");
            CSV_Writer csv2("C:\\Software Development\\6 DOF SIM\\Data\\Kalman_data.csv");
            csv2.write_line("time,tpx,px,mpx,tvx,vx,mbax");
            CSV_Writer csv3("C:\\Software Development\\6 DOF SIM\\Data\\Target_Kalman_data.csv");
            csv3.write_line("time,tpx,px,tvx,vx,tax,ax");
            CSV_Writer csv4("C:\\Software Development\\6 DOF SIM\\Data\\Aero_data.csv");
            csv4.write_line("time,CD,CC,CL,Cl,Cm,Cn,X,Y,Z,L,M,N");
            CSV_Writer csv5("C:\\Software Development\\6 DOF SIM\\Data\\Collision_Course_data.csv");
            csv5.write_line("time,tgo,thetacc,psicc,theta_los_dot,psi_los_dot");

            std::cout<<"Initializing CSV Files\n"<<std::endl;

            Missile missile(i_mx,i_my,i_mz,i_mvx,i_mvy,i_mvz,i_max,i_may,i_maz,i_roll,i_pitch,i_yaw);
            Target target(i_tx,i_ty,i_tz,i_tvx,i_tvy,i_tvz,i_tax,i_tay,i_taz);

            std::cout<<"Missile and Target Objects Initialized\n"<<std::endl;

            EKF ekf(EKF::INERTIAL);
            ekf.pos = Vector3d(i_mx,i_my,i_mz);
            ekf.vel = Vector3d(i_mvx,i_mvy,i_mvz);

            EKF tekf(EKF::TERMINAL);
            tekf.x_term.segment<3>(0) = Vector3d(i_tx, i_ty, i_tz); 
            tekf.x_term.segment<3>(3) = Vector3d(i_tvx, i_tvy, i_tvz);  
            tekf.x_term.segment<3>(6) = Vector3d::Zero();

            std::cout<<"Kalman Filters Initialized\n"<<std::endl;

            AeroCoefficientTables tables;

            int step = 0;

            double intercept_radius = 0.1;   
            int patience_steps = 100;
            int steps_since_min = 0;

            double min_R = 1000000;
            bool cpa_reached = false;

            std::cout<<"World Data\n"<<std::endl;

            current_weather(latitude,longitude);
            std::cout<<"Country: " << country << std::endl;
            std::cout<<"City: " << city << std::endl;
            std::cout<<"Weather: "<<weather<<std::endl;
            std::cout<<"Air Pressure at Sea Level: " << air_pressure_sea_lvl << std::endl;
            std::cout<<"Wind Speed: " << wind_speed << std::endl;
            std::cout<<"Wind Direction: " << wind_direction << std::endl;

            std::cout<<"\n"<<std::endl;

            while (step < steps && !cpa_reached){
                double time = step * delta_t;
                std::string ekfdata;
                std::string tekfdata;
                std::string aerodata;

                Vector3d missile_pos_vec = {missile.m_x,missile.m_y,missile.m_z};
                Vector3d missile_vel_vec = {missile.m_vx,missile.m_vy,missile.m_vz};
                Vector3d missile_acc_vec = {missile.m_ax,missile.m_ay,missile.m_az};
                Vector3d target_pos_vec = {target.t_x,target.t_y,target.t_z};
                Vector3d target_vel_vec = {target.t_vx,target.t_vy,target.t_vz};
                Vector3d target_acc_vec = {target.t_ax,target.t_ay,target.t_az};
                Vector3d body_rot_vec = missile.body_rot_vec;

                std::tie(est_m_pos,est_m_vel,est_m_acc,est_t_pos,est_t_vel,est_t_acc,ekfdata,tekfdata) = missile.SENSORS(missile_pos_vec,missile_vel_vec,missile_acc_vec,target_pos_vec,target_vel_vec,target_acc_vec,ekf,step,rng,nz,tekf,body_rot_vec);
                std::tie(R,acc_vec,euler_angles,quat_vec,executed_body_a_vec,body_rot_vec) = missile.INTELLIGENT_INTEGRATED_GUIDANCE_NAVIGATION_AND_CONTROL_SYSTEM(est_m_pos,est_m_vel,est_m_acc,est_t_pos,est_t_vel,est_t_acc,missile_pos_vec,missile_vel_vec,missile_acc_vec,target_pos_vec,body_rot_vec,fid);
                std::tie(thrust,fuel_mass,defl_angls,aero_forces,aero_moments,new_p,new_v,new_omega,aerodata) = missile.AIRFRAME(missile_pos_vec,missile_vel_vec,executed_body_a_vec,quat_vec,body_rot_vec,tables);
                std::tie(T_go,theta_cc_mt,psi_cc_mt,theta_los_dot,psi_los_dot) = missile.COLLISION_COURSE_GEOMETRY(missile_vel_vec,target_vel_vec,rel_pos_vec,rel_vel_vec,R);
                
                if (fid == TRIM){
                    std::tie(missile.m_ax,missile.m_ay,missile.m_az,missile.m_vx,missile.m_vy,missile.m_vz,missile.m_x,missile.m_y,missile.m_z) = missile.update_params(acc_vec,delta_t);
                }
                else if (fid == FULL_BODY){
                    missile.m_x = new_p.x();
                    missile.m_y = new_p.y();
                    missile.m_z = new_p.z();
                    missile.m_vx = new_v.x();
                    missile.m_vy = new_v.y();
                    missile.m_vz = new_v.z();

                    missile.body_rot_vec = new_omega;
                }

                missile.roll = euler_angles[0];
                missile.pitch = euler_angles[1];
                missile.yaw = euler_angles[2];

                std::tie(target.t_vx,target.t_vy,target.t_vz,target.t_x,target.t_y,target.t_z) = target.update_params(delta_t);

                step += 1;

                std::string data1 = std::to_string(time)+","+std::to_string(R)+","+std::to_string(missile.m_x)+","+std::to_string(missile.m_y)+","+std::to_string(missile.m_z);
                std::string data2 = ","+std::to_string(target.t_x)+","+std::to_string(target.t_y)+","+std::to_string(target.t_z);
                std::string data3 = ","+std::to_string(quat_vec.w())+","+std::to_string(quat_vec.x())+","+std::to_string(quat_vec.y())+","+std::to_string(quat_vec.z())+","+std::to_string(sqrt(pow(quat_vec.w(),2)+pow(quat_vec.x(),2)+pow(quat_vec.y(),2)+pow(quat_vec.z(),2)));
                std::string data4 = ","+std::to_string(thrust)+","+std::to_string(fuel_mass);
                std::string data5 = ","+std::to_string(defl_angls[0])+","+std::to_string(defl_angls[1])+","+std::to_string(defl_angls[2])+","+std::to_string(defl_angls[3]);
                std::string data = data1+data2+data3+data4+data5;
                csv.write_line(data);
                csv2.write_line(std::to_string(time)+","+ekfdata);
                csv3.write_line(std::to_string(time)+","+tekfdata);
                csv4.write_line(std::to_string(time)+","+aerodata);
                csv5.write_line(std::to_string(time)+","+std::to_string(T_go)+","+std::to_string(theta_cc_mt)+","+std::to_string(psi_cc_mt)+","+std::to_string(theta_los_dot)+","+std::to_string(psi_los_dot));

                if (R <= intercept_radius) {
                    std::cout<<"Engagement Results\n"<<std::endl;
                    std::cout << "Intercept (within radius) at t=" << time << "s, range=" << R << " m\n";
                    cpa_reached = true;
                    break;
                }

                if (R < min_R) {
                    min_R = R;
                    steps_since_min = 0;          
                } else {
                    ++steps_since_min;              
                }

                if (steps_since_min >= patience_steps || !std::isfinite(missile.m_x)) {
                    std::cout<<"Engagement Results\n"<<std::endl;
                    std::cout << "Closest point of approach detected at t=" << time
                            << " s (min_range=" << min_R <<"m)\n";
                    cpa_reached = true;
                    break;
                }
            }
            if (!cpa_reached) {
                std::cout<<"Engagement Results\n"<<std::endl;
                std::cout << "Simulation ended by max_steps or other stop condition. min_range=" << min_R << " m\n";
            }
            csv.close_file();
            csv2.close_file();
            csv3.close_file();
            csv4.close_file();
            csv5.close_file();
        }

        void run_Monte_Carlo(Mode mode_,Fidelity fid_,int n){
            CSV_Writer csv("C:\\Software Development\\6 DOF Sim\\Data\\Monte_Carlo_Data.csv");
            csv.write_line("miss");

            for (int i = 0; i < n; ++i){
                setState();

                Missile missile(i_mx,i_my,i_mz,i_mvx,i_mvy,i_mvz,i_max,i_may,i_maz,i_roll,i_pitch,i_yaw);
                Target target(i_tx,i_ty,i_tz,i_tvx,i_tvy,i_tvz,i_tax,i_tay,i_taz);

                EKF ekf(EKF::INERTIAL);
                ekf.pos = Vector3d(i_mx,i_my,i_mz);
                ekf.vel = Vector3d(i_mvx,i_mvy,i_mvz);

                EKF tekf(EKF::TERMINAL);
                tekf.x_term.segment<3>(0) = Vector3d(i_tx, i_ty, i_tz); 
                tekf.x_term.segment<3>(3) = Vector3d(i_tvx, i_tvy, i_tvz);  
                tekf.x_term.segment<3>(6) = Vector3d::Zero();

                AeroCoefficientTables tables;

                int step = 0;

                double intercept_radius = 0.1;   
                int patience_steps = 100;
                int steps_since_min = 0;
                double min_R = 1000000;
                bool cpa_reached = false;

                while (step < steps && !cpa_reached){
                    std::string ekfdata;
                    std::string tekfdata;
                    std::string aerodata;

                    Vector3d missile_pos_vec = {missile.m_x,missile.m_y,missile.m_z};
                    Vector3d missile_vel_vec = {missile.m_vx,missile.m_vy,missile.m_vz};
                    Vector3d missile_acc_vec = {missile.m_ax,missile.m_ay,missile.m_az};
                    Vector3d target_pos_vec = {target.t_x,target.t_y,target.t_z};
                    Vector3d target_vel_vec = {target.t_vx,target.t_vy,target.t_vz};
                    Vector3d target_acc_vec = {target.t_ax,target.t_ay,target.t_az};
                    Vector3d body_rot_vec = missile.body_rot_vec;

                    std::tie(est_m_pos,est_m_vel,est_m_acc,est_t_pos,est_t_vel,est_t_acc,ekfdata,tekfdata) = missile.SENSORS(missile_pos_vec,missile_vel_vec,missile_acc_vec,target_pos_vec,target_vel_vec,target_acc_vec,ekf,step,rng,nz,tekf,body_rot_vec);
                    std::tie(R,acc_vec,euler_angles,quat_vec,executed_body_a_vec,body_rot_vec) = missile.INTELLIGENT_INTEGRATED_GUIDANCE_NAVIGATION_AND_CONTROL_SYSTEM(est_m_pos,est_m_vel,est_m_acc,est_t_pos,est_t_vel,est_t_acc,missile_pos_vec,missile_vel_vec,missile_acc_vec,target_pos_vec,body_rot_vec,fid);
                    std::tie(thrust,fuel_mass,defl_angls,aero_forces,aero_moments,new_p,new_v,new_omega,aerodata) = missile.AIRFRAME(missile_pos_vec,missile_vel_vec,executed_body_a_vec,quat_vec,body_rot_vec,tables);

                    if (fid == TRIM){
                        std::tie(missile.m_ax,missile.m_ay,missile.m_az,missile.m_vx,missile.m_vy,missile.m_vz,missile.m_x,missile.m_y,missile.m_z) = missile.update_params(acc_vec,delta_t);
                    }
                    else if (fid == FULL_BODY){
                        missile.m_x = new_p.x();
                        missile.m_y = new_p.y();
                        missile.m_z = new_p.z();
                        missile.m_vx = new_v.x();
                        missile.m_vy = new_v.y();
                        missile.m_vz = new_v.z();

                        missile.body_rot_vec = new_omega;
                    }

                    missile.roll = euler_angles[0];
                    missile.pitch = euler_angles[1];
                    missile.yaw = euler_angles[2];

                    std::tie(target.t_vx,target.t_vy,target.t_vz,target.t_x,target.t_y,target.t_z) = target.update_params(delta_t);

                    step += 1;

                    if (R <= intercept_radius) {
                        csv.write_line(std::to_string(R));
                        cpa_reached = true;
                        break;
                    }

                    if (R < min_R) {
                        min_R = R;
                        steps_since_min = 0;          
                    } else {
                        ++steps_since_min;              
                    }

                    if (steps_since_min >= patience_steps || !std::isfinite(missile.m_x)) {
                        csv.write_line(std::to_string(min_R));
                        cpa_reached = true;
                        break;
                    }
                }
                if (!cpa_reached) {
                    csv.write_line(std::to_string(min_R));
                }
            }
            csv.close_file();
        }

        void run_Monte_Optimized(){
            
        }
};