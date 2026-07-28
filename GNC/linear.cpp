#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <IGNC.cpp>

#include <iomanip>
#include <limits>

using namespace Eigen;

class Prop: public Enviornment{
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

class Fin_Act: public Enviornment{       
    private:
        double angle_of_attack(double mv_x_b,double mv_z_b){
            double alpha = atan2(mv_z_b,mv_x_b);
            return alpha;
        }

        double side_slip(double mv_x_b,double mv_y_b){
            double beta = atan2(mv_y_b,mv_x_b);
            return beta;
        }

    public:
        std::tuple<double,double> CONTROL_SURFACE_DEFLECTIONS(double mz,double mvx,double mvy,double mvz){
            double alpha = angle_of_attack(mvx,mvz);
            double beta = side_slip(mvx,mvy);
            return {alpha,beta};
        }
};

class Aero: public Enviornment{
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
            double L = q*cross_area*2*body_radius*Cl;
            return L;
        }

        double aero_moment_M(double q,double cross_area){
            double M = q*cross_area*body_length*Cm;
            return M;
        }

        double aero_moment_N(double q,double cross_area){
            double N = q*cross_area*2*body_radius*Cn;
            return N;
        }

        std::tuple<Vector3d,Vector3d> AERODYNAMICS(std::array<double,3> fuel_cm,Vector3d m_vel_arr,Vector3d m_pos_arr,double alpha,double beta,
                                                                            double thrust,AeroCoefficientTables& table,Vector3d defl_angls,Vector3d body_rot_vec){
            double cm = center_of_mass(fuel_cm);
            double q,Mach;
            double X,Y,Z,L,M,N;
            double delta_y,delta_p,delta_r;
            double p_hat,q_hat,r_hat;

            std::tie(q,Mach) = q_M(m_pos_arr,m_vel_arr);

            delta_y = defl_angls[0];
            delta_p = defl_angls[1];
            delta_r = defl_angls[2];
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

            return {aero_forces,aero_moments};
        }
};

class Dyn: public Enviornment{
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

        std::tuple<Vector3d,Vector3d,Vector3d> DYNAMICS(std::array<double,3> fuel_cm,Vector3d v_vec,
                Vector3d aero_forces,Vector3d aero_moments,Vector3d pos_vec,Quaterniond quat_vec,Vector3d body_rot_vec){
            double cm = center_of_mass(fuel_cm);
            double Ixx = I_xx();
            double Iyy = I_yy_zz(cm);
            double Izz = Iyy;
            Vector3d dt_v_b = change_in_vel_body(body_rot_vec,v_vec,aero_forces,quat_vec);
            Vector3d dt_body_rate = change_in_body_rates(body_rot_vec,aero_moments,Ixx,Iyy,Izz);
            Vector3d p_dot = velocity_transform_b_f(v_vec,quat_vec);
            return {dt_v_b,dt_body_rate,p_dot};
        }
};

class Airfr: public Enviornment{
    public:
        Fin_Act fin_actuator;
        Prop propulsion;
        Dyn dynamics;
        Aero aerodynamics;

        std::tuple<double,double,std::array<double,3>> JET_ENGINE(double mz){
            std::tie(thrust,fuel_mass,f_cm) = propulsion.PROPULSION(mz);
            return {thrust,fuel_mass,f_cm};
        }

        std::tuple<double,double> CONTROL_SURFACES(double mz,double mvx,double mvy,double mvz){
            std::tie(alpha,beta) = fin_actuator.CONTROL_SURFACE_DEFLECTIONS(mz,mvx,mvy,mvz);
            return {alpha,beta};
        }

        std::tuple<Vector3d,Vector3d> AERO(Vector3d pos_vec,Vector3d vel_vec,double thrust,Vector3d body_rot_vec,AeroCoefficientTables& table,Vector3d defl_angls){
            std::tie(aero_forces,aero_moments) = aerodynamics.AERODYNAMICS(f_cm,vel_vec,pos_vec,alpha,beta,thrust,table,defl_angls,body_rot_vec);
            return {aero_forces,aero_moments};
        }

        std::tuple<Vector3d,Vector3d,Vector3d> DYN(Vector3d v_vec,Vector3d p_vec,Quaterniond quat_vec,Vector3d body_rot_vec){
            Vector3d dt_v_b;
            Vector3d dt_body_rate;
            Vector3d p_dot;
            std::tie(dt_v_b,dt_body_rate,p_dot) = dynamics.DYNAMICS(f_cm,v_vec,aero_forces,aero_moments,p_vec,quat_vec,body_rot_vec);
            return {dt_v_b,dt_body_rate,p_dot};
        }
};

Quaterniond Quaternion_Vector(double roll, double pitch, double yaw){
    Quaterniond quat_vec;
    double q1_i = cos(roll/2)*cos(pitch/2)*cos(yaw/2)+sin(roll/2)*sin(pitch/2)*sin(yaw/2);
    double q2_i = sin(roll/2)*cos(pitch/2)*cos(yaw/2)-cos(roll/2)*sin(pitch/2)*sin(yaw/2);
    double q3_i = cos(roll/2)*sin(pitch/2)*cos(yaw/2)+sin(roll/2)*cos(pitch/2)*sin(yaw/2);
    double q4_i = cos(roll/2)*cos(pitch/2)*sin(yaw/2)-sin(roll/2)*sin(pitch/2)*cos(yaw/2);
    quat_vec = Quaterniond(q1_i,q2_i,q3_i,q4_i);
    return quat_vec;
}

std::tuple<Vector3d,Vector3d,Vector3d,Vector3d> AIRFRAME(Airfr airframe,Vector3d missile_pos_vec,Vector3d missile_vel_vec,Vector3d body_rot_vec,Vector3d euler,Vector3d defl_angls,AeroCoefficientTables& table){
    double thrust,fuel_mass,alpha,beta;
    std::array<double,3> f_cm;
    Vector3d aero_forces,aero_moments,dt_v_b,dt_body_rate,p_dot;
    Quaterniond quat_vec = Quaternion_Vector(euler[0],euler[1],euler[2]);
    std::tie(thrust,fuel_mass,f_cm) = airframe.JET_ENGINE(missile_pos_vec[2]);
    std::tie(alpha,beta) = airframe.CONTROL_SURFACES(missile_pos_vec[2],missile_vel_vec[0],missile_vel_vec[1],missile_vel_vec[2]);
    std::tie(aero_forces,aero_moments) = airframe.AERO(missile_pos_vec,missile_vel_vec,thrust,body_rot_vec,table,defl_angls);
    std::tie(dt_v_b,dt_body_rate,p_dot) = airframe.DYN(missile_vel_vec,missile_pos_vec,quat_vec,body_rot_vec);

    double phi   = euler(0);
    double theta = euler(1);

    Matrix3d T;

    T <<
    1, sin(phi)*tan(theta), cos(phi)*tan(theta),
    0, cos(phi),           -sin(phi),
    0, sin(phi)/cos(theta), cos(phi)/cos(theta);

    Vector3d euler_dot = T * body_rot_vec;

    return {dt_v_b,dt_body_rate,p_dot,euler_dot};
}

VectorXd State_Derivative(VectorXd& x,VectorXd& u,AeroCoefficientTables& table){
    VectorXd x_dot(12);

    Airfr airframe;
    Vector3d dt_v_b,dt_body_rate,p_dot,euler_dot;
    Vector3d mpos = {x(6),x(7),x(8)};
    Vector3d mvel = {x(0),x(1),x(2)};
    Vector3d brot = {x(3),x(4),x(5)};
    Vector3d euler = {x(9),x(10),x(11)};
    Vector3d defl = {u(0),u(1),u(2)};
    
    std::tie(dt_v_b,dt_body_rate,p_dot,euler_dot) = AIRFRAME(airframe,mpos,mvel,brot,euler,defl,table);

    x_dot <<
        dt_v_b(0),
        dt_v_b(1),
        dt_v_b(2),

        dt_body_rate(0),
        dt_body_rate(1),
        dt_body_rate(2),

        p_dot(0),
        p_dot(1),
        p_dot(2),

        euler_dot(0),
        euler_dot(1),
        euler_dot(2);

    return x_dot;
}

std::tuple<MatrixXd,MatrixXd> Calculate_Jacobians(VectorXd& x,VectorXd& u,AeroCoefficientTables table){
    int n = x.size();
    int m = u.size();

    MatrixXd A(n,n);
    MatrixXd B(n,m);

    for (int i=0; i<n; i++){
        double h = std::max(1e-6, 1e-6 * std::abs(x(i)));

        VectorXd xp = x;
        VectorXd xm = x;

        xp(i) += h;
        xm(i) -= h;

        VectorXd fp = State_Derivative(xp,u,table);
        VectorXd fm = State_Derivative(xm,u,table);

        A.col(i) = (fp-fm)/(2.0*h);
    }

    for(int i=0;i<m;i++){
        double h = std::max(1e-6, 1e-6 * std::abs(u(i)));

        VectorXd up = u;
        VectorXd um = u;

        up(i)+=h;
        um(i)-=h;

        VectorXd fp = State_Derivative(x,up,table);
        VectorXd fm = State_Derivative(x,um,table);

        B.col(i)= (fp-fm)/(2.0*h);

        double tol = 1e-10;

        A = A.unaryExpr([tol](double x)
        {
            return (std::abs(x) < tol) ? 0.0 : x;
        });

        B = B.unaryExpr([tol](double x)
        {
            return (std::abs(x) < tol) ? 0.0 : x;
        });
            }

    return {A,B};
}

std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd> Long_Lat_Matrices(MatrixXd A,MatrixXd B){
    std::vector<int> lon = {0,2,4,10};
    std::vector<int> lat = {1,3,5,9};

    MatrixXd A_lon(4,4);
    MatrixXd B_lon(4,1);
    MatrixXd A_lat(4,4);
    MatrixXd B_lat(4,2);

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++)
        {
            A_lon(i,j)=A(lon[i],lon[j]);
            A_lat(i,j)=A(lat[i],lat[j]);
        }
    }

    for(int i=0;i<4;i++){
        B_lon(i,0)=B(lon[i],1);
        B_lat(i,0)=B(lat[i],0);
        B_lat(i,1)=B(lat[i],2);
    }

    return {A_lon,A_lat,B_lon,B_lat};
}

MatrixXd Modal_Control_Matrix(MatrixXcd& B_modal,bool norm = true){
    MatrixXd C(B_modal.rows(), B_modal.cols());

    for(int i=0;i<B_modal.rows();i++)
    {
        for(int j=0;j<B_modal.cols();j++)
        {
            C(i,j)=std::abs(B_modal(i,j));
        }
    }

    if (norm){
        for(int i=0;i<C.rows();i++)
        {
            double rowNorm=C.row(i).norm();

            if(rowNorm>1e-12)
                C.row(i)/=rowNorm;
        }
    }

    return C;
}

void PrintModalControlEffectiveness(MatrixXd& C,const std::vector<std::string>& modeNames,const std::vector<std::string>& controlNames){
    using std::cout;
    using std::setw;
    using std::left;
    using std::fixed;
    using std::setprecision;

    cout << "\n=============================================================\n";
    cout << "Modal Control Effectiveness Matrix\n";
    cout << "=============================================================\n";

    cout << left << setw(22) << "Mode";

    for(const auto& control : controlNames)
        cout << setw(15) << control;

    cout << "\n-------------------------------------------------------------\n";

    for(int i=0;i<C.rows();i++)
    {
        cout << left << setw(22) << modeNames[i];

        for(int j=0;j<C.cols();j++)
        {
            cout << setw(15)
                 << fixed << setprecision(4)
                 << C(i,j);
        }

        cout << '\n';
    }

    cout << "=============================================================\n";
}

void Level_Flight_Stability_Analysis(VectorXd& x,VectorXd& u,AeroCoefficientTables& tables){
    MatrixXd A, B, A_lon, A_lat, B_lon, B_lat;

    std::tie(A, B) = Calculate_Jacobians(x, u, tables);

    std::cout << "A matrix\n" << A << '\n';
    std::cout << '\n';
    std::cout << "B matrix\n" << B << '\n';
    std::cout << '\n';

    std::tie(A_lon,A_lat,B_lon,B_lat) = Long_Lat_Matrices(A,B);

    std::cout << "A Longitudinal matrix\n" << A_lon << '\n';
    std::cout << '\n';
    std::cout << "B Longitudinal matrix\n" << B_lon << '\n';
    std::cout << '\n';
    std::cout << "A Lateral matrix\n" << A_lat << '\n';
    std::cout << '\n';
    std::cout << "B lateral matrix\n" << B_lat << '\n';
    std::cout << '\n';

    EigenSolver<MatrixXd> lon(A_lon);
    EigenSolver<MatrixXd> lat(A_lat);

    auto lon_eigvals = lon.eigenvalues();
    auto lon_eigvecs = lon.eigenvectors();
    auto lat_eigvals = lat.eigenvalues();
    auto lat_eigvecs = lat.eigenvectors();

    MatrixXcd lon_inv = lon_eigvecs.inverse();
    MatrixXcd lat_inv = lat_eigvecs.inverse();

    MatrixXcd lon_modal = lon_inv*B_lon;
    MatrixXcd lat_modal = lat_inv*B_lat;

    MatrixXd lon_modal_control = Modal_Control_Matrix(lon_modal,false);
    MatrixXd lat_modal_control = Modal_Control_Matrix(lat_modal);

    MatrixXcd lon_participation_factors(lon_eigvecs.rows(), lon_eigvecs.cols());
    for(int mode = 0; mode < lon_eigvecs.cols(); mode++)
    {
        for(int state = 0; state < lon_eigvecs.rows(); state++)
        {
            lon_participation_factors(state, mode) =
                lon_eigvecs(state, mode) * lon_inv(mode, state);
        }
    }

    MatrixXcd lat_participation_factors(lat_eigvecs.rows(), lat_eigvecs.cols());
    for(int mode = 0; mode < lat_eigvecs.cols(); mode++)
    {
        for(int state = 0; state < lat_eigvecs.rows(); state++)
        {
            lat_participation_factors(state, mode) =
                lat_eigvecs(state, mode) * lat_inv(mode, state);
        }
    }

    std::vector<std::string> lonNames = {"u","w","q","theta"};
    std::vector<std::string> latNames = {"v","p","r","phi"};

    std::cout << "Longitudinal Eigen Analysis\n";
    std::cout << "\n";
    std::cout << "=================================================================================================================\n";
    std::cout << std::left
            << std::setw(24) << "Mode"
            << std::setw(24) << "Eigenvalue"
            << std::setw(10) << "wn"
            << std::setw(10) << "wd"
            << std::setw(10) << "zeta"
            << std::setw(12) << "Tau (s)"
            << std::setw(12) << "Period (s)\n";
    std::cout << "=================================================================================================================\n";

    for (int i = 0; i < lon_eigvals.size(); ++i) {
        std::complex<double> lambda = lon_eigvals(i);

        double sigma = lambda.real();
        double wd    = std::abs(lambda.imag());

        double wn = std::sqrt(sigma*sigma + wd*wd);

        double zeta = (wn > 1e-10) ? -sigma/wn : 1.0;

        double tau = (std::abs(sigma) > 1e-10) ? -1.0/sigma : std::numeric_limits<double>::infinity();

        double period = (wd > 1e-10) ? 2.0*3.14159/wd : std::numeric_limits<double>::infinity();

        Eigen::VectorXcd v = lon_eigvecs.col(i);

        v /= v.cwiseAbs().maxCoeff();

        std::stringstream dominant;
        bool first = true;

        const double threshold = 0.01;  

        for (int j = 0; j < v.size(); ++j)
        {
            double mag = std::abs(v(j));

            if (mag < threshold)
                continue;

            if (!first)
                dominant << ", ";

            dominant << lonNames[j];

            first = false;
        }

        std::string mode;

        if (wd > 1e-3)
        {
            if (wn > 5.0)
                mode = "Short Period";
            else
                mode = "Phugoid";
        }
        else
        {
            mode = "Real Mode";
        }

        std::stringstream eigString;

        eigString << std::fixed << std::setprecision(3)
                << sigma;

        if(lambda.imag() >= 0)
            eigString << " + ";
        else
            eigString << " - ";

        eigString << wd << "i";

        std::cout << std::left
                << std::setw(24) << mode
                << std::setw(24) << eigString.str()
                << std::setw(10) << std::fixed << std::setprecision(3) << wn
                << std::setw(10) << wd
                << std::setw(10) << zeta;

        if (abs(tau) < 1e6)
            std::cout << std::setw(12) << tau;
        else
            std::cout << std::setw(12) << "-";

        if (abs(period) < 1e6)
            std::cout << std::setw(12) << period;
        else
            std::cout << std::setw(12) << "-";

        std::cout << '\n';
    }
    std::cout << "=================================================================================================================\n";

    std::cout << '\n';
    std::cout << "=============================================================\n";
    std::cout << "Longitudinal Participation Factors\n";
    std::cout << "=============================================================\n";

    const char* states[] = {"u", "w", "q", "theta"};

    std::cout << std::setw(12) << "";

    for(int mode = 0; mode < lon_participation_factors.cols(); mode++)
        std::cout << std::setw(12) << ("Mode " + std::to_string(mode+1));

    std::cout << "\n-------------------------------------------------------------\n";

    for(int state = 0; state < lon_participation_factors.rows(); state++)
    {
        std::cout << std::setw(12) << states[state];

        for(int mode = 0; mode < lon_participation_factors.cols(); mode++)
        {
            std::cout << std::setw(12)
                    << std::fixed << std::setprecision(3)
                    << std::abs(lon_participation_factors(state,mode));
        }

        std::cout << '\n';
    }

    std::vector<std::string> lonModes =
    {
        "Short Period",
        "Short Period",
        "Real Mode",
        "Real Mode"
    };

    std::vector<std::string> lonControls =
    {
        "Pitch Defl."
    };

    PrintModalControlEffectiveness(lon_modal_control, lonModes, lonControls);

    std::cout << '\n';
    std::cout << "Lateral Eigen Analysis\n";
    std::cout << "\n";
    std::cout << "=================================================================================================================\n";
    std::cout << std::left
            << std::setw(24) << "Mode"
            << std::setw(24) << "Eigenvalue"
            << std::setw(10) << "wn"
            << std::setw(10) << "wd"
            << std::setw(10) << "zeta"
            << std::setw(12) << "Tau (s)"
            << std::setw(12) << "Period (s)\n";
    std::cout << "=================================================================================================================\n";

    for (int i = 0; i < lat_eigvals.size(); ++i) {
        std::complex<double> lambda = lat_eigvals[i];

        double sigma = lambda.real();
        double wd    = std::abs(lambda.imag());

        double wn = std::sqrt(sigma*sigma + wd*wd);

        double zeta = (wn > 1e-10) ? -sigma/wn : 1.0;

        double tau = (std::abs(sigma) > 1e-10) ? -1.0/sigma : std::numeric_limits<double>::infinity();

        double period = (wd > 1e-10) ? 2.0*3.14159/wd : std::numeric_limits<double>::infinity();

        Eigen::VectorXcd v = lat_eigvecs.col(i);

        v /= v.cwiseAbs().maxCoeff();

        std::string mode;

        if (wd > 1e-3)
        {
            mode = "Dutch Roll";
        }
        else
        {
            if (std::abs(sigma) > 1.0)
                mode = "Roll";
            else
                mode = "Spiral";
        }

        std::stringstream eigString;

        eigString << std::fixed << std::setprecision(3)
                << sigma;

        if(lambda.imag() >= 0)
            eigString << " + ";
        else
            eigString << " - ";

        eigString << wd << "i";

        std::cout << std::left
                << std::setw(24) << mode
                << std::setw(24) << eigString.str()
                << std::setw(10) << std::fixed << std::setprecision(3) << wn
                << std::setw(10) << wd
                << std::setw(10) << zeta;

        if (abs(tau) < 1e6)
            std::cout << std::setw(12) << tau;
        else
            std::cout << std::setw(12) << "-";

        if (abs(period) < 1e6)
            std::cout << std::setw(12) << period;
        else
            std::cout << std::setw(12) << "-";
        std::cout << '\n';
    }
    std::cout << "=================================================================================================================\n";

    std::cout << '\n';
    std::cout << "=============================================================\n";
    std::cout << "Lateral Participation Factors\n";
    std::cout << "=============================================================\n";

    const char* latstates[] = {"v", "p", "r", "phi"};

    std::cout << std::setw(12) << "";

    for(int mode = 0; mode < lat_participation_factors.cols(); mode++)
        std::cout << std::setw(12) << ("Mode " + std::to_string(mode+1));

    std::cout << "\n-------------------------------------------------------------\n";

    for(int state = 0; state < lat_participation_factors.rows(); state++)
    {
        std::cout << std::setw(12) << latstates[state];

        for(int mode = 0; mode < lat_participation_factors.cols(); mode++)
        {
            std::cout << std::setw(12)
                    << std::fixed << std::setprecision(3)
                    << std::abs(lat_participation_factors(state,mode));
        }

        std::cout << '\n';
    }

    std::vector<std::string> latModes =
    {
        "Dutch Roll",
        "Dutch Roll",
        "Roll",
        "Spiral"
    };

    std::vector<std::string> latControls =
    {
        "Yaw Defl.",
        "Roll Defl."
    };

    PrintModalControlEffectiveness(lat_modal_control, latModes, latControls);
}

void Turning_Stability_Analysis(VectorXd& x,VectorXd& u,AeroCoefficientTables& tables){
    MatrixXd A, B;

    std::tie(A,B) = Calculate_Jacobians(x, u, tables);

    std::cout << "A matrix\n" << A << '\n';
    std::cout << '\n';
    std::cout << "B matrix\n" << B << '\n';
    std::cout << '\n';

    EigenSolver<MatrixXd> Amat(A);

    auto eigvals = Amat.eigenvalues();
    auto eigvecs = Amat.eigenvectors();

    MatrixXcd vec_inv = eigvecs.inverse();

    MatrixXcd modal = vec_inv*B;

    MatrixXd turn_modal = Modal_Control_Matrix(modal);

    MatrixXcd participation_factors(eigvecs.rows(), eigvecs.cols());
    for(int mode = 0; mode < eigvecs.cols(); mode++)
    {
        for(int state = 0; state < eigvecs.rows(); state++)
        {
            participation_factors(state, mode) =
                eigvecs(state, mode) * vec_inv(mode, state);
        }
    }

    std::vector<std::string> stateNames = {"u","v","w","p","q","r","x","y","z","phi","theta","psi"};

    std::cout << "Eigen Analysis\n";
    std::cout << "\n";
    std::cout << "===============================================================================================================================================\n";
    std::cout << std::left
            << std::setw(24) << "Mode"
            << std::setw(24) << "Eigenvalue"
            << std::setw(10) << "wn"
            << std::setw(10) << "wd"
            << std::setw(10) << "zeta"
            << std::setw(12) << "Tau (s)"
            << std::setw(12) << "Period (s)\n";
    std::cout << "===============================================================================================================================================\n";

    for (int i = 0; i < eigvals.size(); i++)
    {
        std::complex<double> lambda = eigvals(i);

        double sigma = lambda.real();
        double wd    = std::abs(lambda.imag());
        double wn    = std::sqrt(sigma * sigma + wd * wd);

        double zeta = (wn > 1e-8) ? -sigma / wn : 1.0;
        double tau = (std::abs(sigma) > 1e-8) ? -1.0 / sigma
                                            : std::numeric_limits<double>::infinity();
        double period = (wd > 1e-8) ? 2.0 * 3.14159 / wd
                                    : std::numeric_limits<double>::infinity();

        Eigen::VectorXcd eig = eigvecs.col(i);
        eig /= eig.cwiseAbs().maxCoeff();

        std::string mode;

        double U     = std::abs(eig(0));
        double V     = std::abs(eig(1));
        double W     = std::abs(eig(2));

        double P     = std::abs(eig(3));
        double Q     = std::abs(eig(4));
        double R     = std::abs(eig(5));

        double X     = std::abs(eig(6));
        double Y     = std::abs(eig(7));
        double Z     = std::abs(eig(8));

        double Phi   = std::abs(eig(9));
        double Theta = std::abs(eig(10));
        double Psi   = std::abs(eig(11));

        bool oscillatory = wd > 1e-3;

        double transMax = std::max({U,V,W});
        double rateMax  = std::max({P,Q,R});
        double angleMax = std::max({Phi,Theta,Psi});
        double posMax   = std::max({X,Y,Z});

        if (posMax > transMax &&
            posMax > rateMax &&
            posMax > angleMax)
        {
            if (std::abs(sigma) < 1e-4)
                mode = "Position Integrator";
            else if (oscillatory)
                mode = "Trajectory Osc.";
            else
                mode = "Trajectory";
        }
        else if (oscillatory &&
                W == transMax &&
                (W + Q) > (V + P))
        {
            mode = "Short Period";
        }
        else if (oscillatory &&
                V == transMax)
        {
            mode = "Dutch Roll";
        }
        else if (oscillatory &&
                (U + Theta) > (W + Q))
        {
            mode = "Phugoid";
        }
        else if (!oscillatory &&
                P == rateMax)
        {
            mode = "Roll Subsidence";
        }
        else if (!oscillatory &&
                Phi == angleMax)
        {
            mode = "Spiral";
        }
        else if (!oscillatory &&
                U == transMax)
        {
            mode = "Longitudinal";
        }
        else
        {
            mode = "Mixed";
        }

        std::stringstream eigString;

        eigString << std::fixed << std::setprecision(3)
                << sigma;

        if(lambda.imag() >= 0)
            eigString << " + ";
        else
            eigString << " - ";

        eigString << wd << "i";

        std::cout << std::left
                << std::setw(24) << mode
                << std::setw(24) << eigString.str()
                << std::setw(10) << std::fixed << std::setprecision(3) << wn
                << std::setw(10) << wd
                << std::setw(10) << zeta;

        if (abs(tau) < 1e6)
            std::cout << std::setw(12) << tau;
        else
            std::cout << std::setw(12) << "-";

        if (abs(period) < 1e6)
            std::cout << std::setw(12) << period;
        else
            std::cout << std::setw(12) << "-";

        std::cout << '\n';
    }
    std::cout << "===============================================================================================================================================\n";

    std::cout << '\n';
    std::cout << "=======================================================================================================================================================================================\n";
    std::cout << "Participation Factors\n";
    std::cout << "=======================================================================================================================================================================================\n";

    const char* turnstates[] = {"u","v","w","p","q","r","x","y","z","phi","theta","psi"};

    std::cout << std::setw(12) << "";

    for(int mode = 0; mode < participation_factors.cols(); mode++)
        std::cout << std::setw(12) << ("Mode " + std::to_string(mode+1));

    std::cout << "\n---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";

    for(int state = 0; state < participation_factors.rows(); state++)
    {
        std::cout << std::setw(12) << turnstates[state];

        for(int mode = 0; mode < participation_factors.cols(); mode++)
        {
            std::cout << std::setw(12)
                    << std::fixed << std::setprecision(3)
                    << std::abs(participation_factors(state,mode));
        }

        std::cout << '\n';
    }

    std::vector<std::string> Modes =
    {
        "Mode 1",
        "Mode 2",
        "Mode 3",
        "Mode 4",
        "Mode 5",
        "Mode 6",
        "Mode 7",
        "Mode 8",
        "Mode 9",
        "Mode 10",
        "Mode 11",
        "Mode 12"
    };

    std::vector<std::string> Controls =
    {
        "Yaw Defl.",
        "Pitch Defl.",
        "Roll Defl."
    };

    PrintModalControlEffectiveness(turn_modal, Modes, Controls);

}

enum class TrimMode
{
    LEVEL,
    PRESCRIBED_AOA,
    COORDINATED_TURN,
    CLIMB,
    FULL_6DOF
};

struct TrimConditions
{
    double mach = 0.7;
    double altitude = 1000.0;

    double gamma = 0.0;
    double phi   = 0.0;

    double alpha_command = 0.0;
    bool fix_alpha = false;

    double beta_command = 0.0;
    bool fix_beta = false;

    double turn_rate = 0.0;

};

struct TrimState
{
    double alpha = 0;
    double beta  = 0;
    double phi   = 0;
    double theta = 0;

    double deltaPitch = 0;
    double deltaYaw   = 0;
    double deltaRoll  = 0;
};

struct TrimResult
{
    TrimState vars;

    VectorXd state;
    VectorXd control;

    bool converged = false;
    int iterations = 0;
    double residual = 0;
};

class TrimEngine
{
public:

    TrimEngine(AeroCoefficientTables& t)
        : tables(t) {}

    TrimResult Solve(const TrimConditions& c, TrimMode mode)
    {
        TrimState v = InitialGuess(c, mode);

        double lambda = 1e-3;

        TrimResult out;

        for(int iter=0; iter<100; iter++)
        {
            VectorXd R = Residual(v, c, mode);
            double err = R.norm();

            std::cout << "Iter " << iter
                      << " | err = " << err
                      << std::endl;

            if(err < 1e-6)
            {
                out.converged = true;
                break;
            }

            MatrixXd J = Jacobian(v, c, mode);

            MatrixXd A = J.transpose()*J
                       + lambda*MatrixXd::Identity(J.cols(), J.cols());

            VectorXd b = -J.transpose()*R;

            VectorXd dz = A.ldlt().solve(b);

            TrimState v_new = v;
            Apply(v_new, dz, mode);

            double err_new = Residual(v_new, c, mode).norm();

            if(err_new < err)
            {
                v = v_new;
                lambda *= 0.3;
            }
            else
            {
                lambda *= 2.0;
            }

            out.iterations = iter;
            out.residual = err;
        }

        out.vars = v;
        out.state = BuildState(v, c);
        out.control = BuildControl(v);

        return out;
    }

private:

    AeroCoefficientTables& tables;

    TrimState InitialGuess(const TrimConditions& c, TrimMode mode)
    {
        TrimState v;

        double V = c.mach * 343.0;

        v.alpha = 2.0 * 3.14159/180.0;
        v.theta = v.alpha;

        if(mode == TrimMode::COORDINATED_TURN)
            v.beta = c.beta_command;

        return v;
    }

    VectorXd BuildState(const TrimState& v, const TrimConditions& c)
    {
        VectorXd x(12);

        double V = c.mach * 343.0;

        x(0) = V * cos(v.alpha) * cos(v.beta);
        x(1) = V * sin(v.beta);
        x(2) = V * sin(v.alpha) * cos(v.beta);

        x(3) = -c.turn_rate*sin(v.theta);
        x(4) = c.turn_rate*sin(v.phi)*cos(v.theta);
        x(5) = c.turn_rate*cos(v.phi)*cos(v.theta);

        x(6) = 0;
        x(7) = 0;
        x(8) = c.altitude;

        x(9)  = v.phi;
        x(10) = v.theta;
        x(11) = 0;

        return x;
    }

    VectorXd BuildControl(const TrimState& v)
    {
        VectorXd u(3);
        u << v.deltaYaw, v.deltaPitch, v.deltaRoll;
        return u;
    }

    VectorXd Residual(const TrimState& v,
                      const TrimConditions& c,
                      TrimMode mode)
    {
        VectorXd x = BuildState(v,c);
        VectorXd u = BuildControl(v);

        VectorXd xd = State_Derivative(x,u,tables);

        std::vector<double> R;

        auto add = [&](double v){ R.push_back(v); };

        add(xd(2)); // Wdot
        add(xd(4)); // Qdot

        if(mode == TrimMode::FULL_6DOF)
        {
            add(xd(1)); // Vdot
            add(xd(3)); // Pdot
            add(xd(5)); // Rdot
            add(v.theta - v.alpha - c.gamma);
        }

        if(mode == TrimMode::COORDINATED_TURN)
        {
            add(xd(1)); // Vdot
            add(xd(3)); // Pdot
            add(xd(5)); // Rdot
            add(v.theta - v.alpha - c.gamma);
        }

        if(mode == TrimMode::LEVEL)
        {
            add(v.theta - v.alpha);
        }

        if(mode == TrimMode::CLIMB)
        {
            add(v.theta - v.alpha - c.gamma);
        }

        if(mode == TrimMode::PRESCRIBED_AOA)
        {
            add(v.alpha - c.alpha_command);
        }

        VectorXd Rv(R.size());
        for(int i=0;i<R.size();i++) Rv(i)=R[i];

        return Rv;
    }

    MatrixXd Jacobian(const TrimState& v,
                      const TrimConditions& c,
                      TrimMode mode)
    {
        const double eps = 1e-6;

        VectorXd R0 = Residual(v,c,mode);

        int n = 7; // max variables
        MatrixXd J(R0.size(), n);

        for(int i=0;i<n;i++)
        {
            TrimState vp = v;

            ((double*)&vp)[i] += eps;

            VectorXd Rp = Residual(vp,c,mode);

            J.col(i) = (Rp - R0)/eps;
        }

        return J;
    }

    void Apply(TrimState& v,
               const VectorXd& dz,
               TrimMode mode)
    {
        v.alpha      += dz(0);
        v.beta       += dz(1);
        v.phi        += dz(2);
        v.theta      += dz(3);
        v.deltaPitch += dz(4);
        v.deltaYaw   += dz(5);
        v.deltaRoll  += dz(6);
    }
};

int main(){
    AeroCoefficientTables tables;
    double d2r = 3.14159/180.0;

    TrimConditions c;
    c.mach = 2;
    c.altitude = 3000;
    c.gamma = 15.0*d2r;
    //c.phi = 0.0*d2r;
    c.alpha_command = 0.0*d2r;
    c.beta_command = 0.0*d2r;
    c.turn_rate = 2.0*(1/57.3);
    c.fix_alpha = false;
    c.fix_beta = false;

    TrimEngine solver(tables);
    TrimMode mode = TrimMode::COORDINATED_TURN;

    auto res = solver.Solve(c, mode);

    std::cout << "alpha = " << res.vars.alpha*180/3.14159 << "\n";
    std::cout << "beta = " << res.vars.beta*180/3.14159<< "\n";
    std::cout << "phi = " << res.vars.phi*180/3.14159 << "\n";
    std::cout << "theta = " << res.vars.theta*180/3.14159 << "\n";
    std::cout << '\n';
    std::cout << "Trimmed State\n";
    std::cout << res.state << '\n';
    std::cout << "Trimmed Controls\n";
    std::cout << res.control << '\n';
    std::cout << "\n";

    if (mode != TrimMode::COORDINATED_TURN){
        Level_Flight_Stability_Analysis(res.state, res.control,tables);
    } else {
        Turning_Stability_Analysis(res.state, res.control,tables);
    }

    return 0;
}