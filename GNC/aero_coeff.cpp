#include <vector>
#include <iostream>
#include <array>
#include <csv.h>
#include <cmath>

struct CacheEntry{
    int i0;
    double x0,x1;
    double inv_dx;
};

class Interpolation {
    public:
        std::vector<std::vector<double>> axes;
        std::vector<double> data;
        std::vector<int> sizes;
        std::vector<CacheEntry> cache;

        Interpolation() {}

        Interpolation(const std::vector<std::vector<double>>& axes_,
                    const std::vector<double>& data_,
                    const std::vector<int>& sizes_)
            : axes(axes_), data(data_), sizes(sizes_)
        {
            cache.resize(axes.size());

            for (int d = 0; d < axes.size(); d++) {
                cache[d].i0 = 0;
                cache[d].x0 = axes[d][0];
                cache[d].x1 = axes[d][1];
                cache[d].inv_dx = 1.0 / (axes[d][1] - axes[d][0]);
            }
        }

        int findInterval(const std::vector<double>& axis, double x) {
            if (x <= axis.front()) return 0;
            if (x >= axis.back()) return axis.size() - 2;

            auto it = std::lower_bound(axis.begin(), axis.end(), x);
            return std::distance(axis.begin(), it) - 1;
        }

        void updateCache(int d, double x) {
            auto& c = cache[d];

            if (x >= c.x0 && x <= c.x1) return;

            int i0 = findInterval(axes[d], x);

            c.i0 = i0;
            c.x0 = axes[d][i0];
            c.x1 = axes[d][i0 + 1];
            c.inv_dx = 1.0 / (c.x1 - c.x0);
        }

        int getFlatIndex(const std::vector<int>& indices) {
            int index = 0;
            int stride = 1;

            for (int d = sizes.size() - 1; d >= 0; --d) {
                index += indices[d] * stride;
                stride *= sizes[d];
            }

            return index;
        }

        double interpolate(const std::vector<double>& point) {
            int dims = axes.size();

            std::vector<double> t(dims);
            std::vector<int> baseIdx(dims);

            for (int d = 0; d < dims; d++) {
                double x = point[d];

                if (x < axes[d].front()) x = axes[d].front();
                if (x > axes[d].back())  x = axes[d].back();

                updateCache(d, x);

                auto& c = cache[d];
                baseIdx[d] = c.i0;

                double val = (x - c.x0) * c.inv_dx;

                if (val < 0.0) val = 0.0;
                if (val > 1.0) val = 1.0;

                t[d] = val;
            }

            bool exactPoint = true;
            std::vector<int> exactIdx(dims);

            for (int d = 0; d < dims; d++) {
                auto& c = cache[d];

                if (std::abs(point[d] - c.x0) < 1e-12) {
                    exactIdx[d] = c.i0;
                }
                else if (std::abs(point[d] - c.x1) < 1e-12) {
                    exactIdx[d] = c.i0 + 1;
                }
                else {
                    exactPoint = false;
                    break;
                }
            }

            if (exactPoint) {
                return data[getFlatIndex(exactIdx)];
            }

            double result = 0.0;
            int corners = 1 << dims;

            for (int cMask = 0; cMask < corners; cMask++) {
                std::vector<int> idx = baseIdx;
                double weight = 1.0;

                for (int d = 0; d < dims; d++) {
                    if (cMask & (1 << d)) {
                        idx[d]++;
                        weight *= t[d];
                    } else {
                        weight *= (1.0 - t[d]);
                    }
                }

                result += weight * data[getFlatIndex(idx)];
            }

            return result;
        }
};

class AeroCoefficientTables {

private:
    std::vector<std::vector<double>> axes;

    std::vector<int> sizes = {10,10,10};

    std::vector<double> CD;
    std::vector<double> CL;
    std::vector<double> CC;

    std::vector<double> Cl;
    std::vector<double> Cm;
    std::vector<double> Cn;

    std::vector<double> CL_delta_p;
    std::vector<double> Cm_delta_p;

    std::vector<double> CC_delta_y;
    std::vector<double> Cn_delta_y;

    std::vector<double> Cl_delta_r;

public:

    Interpolation CD_interp;
    Interpolation CL_interp;
    Interpolation CC_interp;

    Interpolation Cl_interp;
    Interpolation Cm_interp;
    Interpolation Cn_interp;

    Interpolation CL_delta_p_interp;
    Interpolation Cm_delta_p_interp;

    Interpolation CC_delta_y_interp;
    Interpolation Cn_delta_y_interp;

    Interpolation Cl_delta_r_interp;

    AeroCoefficientTables()
        : axes(3)
    {
        generateAxes();
        generateTables();

        CD_interp = Interpolation(axes, CD, sizes);
        CL_interp = Interpolation(axes, CL, sizes);
        CC_interp = Interpolation(axes, CC, sizes);

        Cl_interp = Interpolation(axes, Cl, sizes);
        Cm_interp = Interpolation(axes, Cm, sizes);
        Cn_interp = Interpolation(axes, Cn, sizes);

        CL_delta_p_interp = Interpolation(axes, CL_delta_p, sizes);
        Cm_delta_p_interp = Interpolation(axes, Cm_delta_p, sizes);
        CC_delta_y_interp = Interpolation(axes, CC_delta_y, sizes);
        Cn_delta_y_interp = Interpolation(axes, Cn_delta_y, sizes);
        Cl_delta_r_interp = Interpolation(axes, Cl_delta_r, sizes);
    }

    void generateAxes(){

        int nMach  = sizes[0];
        int nAlpha = sizes[1];
        int nBeta  = sizes[2];

        axes[0].resize(nMach);
        axes[1].resize(nAlpha);
        axes[2].resize(nBeta);

        for(int i=0; i<nMach; i++)
            axes[0][i] = 0.3 + 2.7*i/(nMach-1);

        for(int i=0; i<nAlpha; i++)
            axes[1][i] = (-25.0 + 50.0*i/(nAlpha-1))*3.14159/180.0;;

        for(int i=0; i<nBeta; i++)
            axes[2][i] = (-25.0 + 50.0*i/(nBeta-1))*3.14159/180.0;;
    }

    void generateTables(){

        int nMach  = sizes[0];
        int nAlpha = sizes[1];
        int nBeta  = sizes[2];

        int totalSize = nMach * nAlpha * nBeta;

        CD.resize(totalSize);

        CL.resize(totalSize);
        CC.resize(totalSize);

        Cl.resize(totalSize);
        Cm.resize(totalSize);
        Cn.resize(totalSize);

        CL_delta_p.resize(totalSize);
        Cm_delta_p.resize(totalSize);

        CC_delta_y.resize(totalSize);
        Cn_delta_y.resize(totalSize);

        Cl_delta_r.resize(totalSize);

        int idx = 0;

        for(int iMach=0;iMach<nMach;iMach++)
        for(int iAlpha=0;iAlpha<nAlpha;iAlpha++)
        for(int iBeta=0;iBeta<nBeta;iBeta++)
        {
            double M = axes[0][iMach];
            double a = axes[1][iAlpha];
            double b = axes[2][iBeta];

            double CD0;

            if(M < 0.9)
                CD0 = 0.035;

            else if(M < 1.2)
                CD0 = 0.035 + 0.12*(M-0.9);

            else
                CD0 = 0.071 + 0.015*(M-1.2);

            CD[idx] = CD0 + 0.60*a*a + 0.80*b*b;
            CL[idx] = 3.5*a - 0.15*a*a*a;
            CC[idx] = -2.2*b;
            Cl[idx] = -0.15*b;
            Cm[idx] = -1.2*a + 0.05*b*b;
            Cn[idx] = 0.90*b;

            CL_delta_p[idx] = 3.8 * (1.0 - 0.15*std::abs(a));
            Cm_delta_p[idx] = -2.8 * (1.0 - 0.10*std::abs(a));
            CC_delta_y[idx] = 2.8 * (1.0 - 0.10*std::abs(b));
            Cn_delta_y[idx] = 2.2 * (1.0 - 0.10*std::abs(b));
            Cl_delta_r[idx] = 2.5;

            idx++;
        }
    }

    std::tuple<double,double,double,double,double,double>interpolate_coefficients(
        double Mach,
        double alpha_deg,
        double beta_deg,
        double delta_pitch_deg,
        double delta_yaw_deg,
        double delta_roll_deg)
    {
        std::vector<double> point = {Mach,alpha_deg,beta_deg};

        double dp = delta_pitch_deg;
        double dy = delta_yaw_deg;
        double dr = delta_roll_deg;

        double CD_b = CD_interp.interpolate(point);
        double CL_b = CL_interp.interpolate(point);
        double CC_b = CC_interp.interpolate(point);
        double Cl_b = Cl_interp.interpolate(point);
        double Cm_b = Cm_interp.interpolate(point);
        double Cn_b = Cn_interp.interpolate(point);

        double CL = CL_b + CL_delta_p_interp.interpolate(point) * dp;
        double CC = CC_b + CC_delta_y_interp.interpolate(point) * dy;
        double Cl = Cl_b + Cl_delta_r_interp.interpolate(point) * dr;
        double Cm = Cm_b + Cm_delta_p_interp.interpolate(point) * dp;
        double Cn = Cn_b + Cn_delta_y_interp.interpolate(point) * dy;

        return {CD_b,CC,CL,Cl,Cm,Cn};
    }
};