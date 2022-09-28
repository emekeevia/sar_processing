#include <iostream>
#include <vector>
#include <complex.h>
#include <complex>
#include <cmath>
#include <string>
#include <fftw3.h>
#include "extra_tools.h"
#include "compression_algorithm.h"
#include "tester.h"

struct Satellite{
    // Sensor parameters (ERS satellite)
    double fs = 4.5e+7;  //Range Sampling Frequency [Hz]
    double K_r = 4.8e+12;     // FM Rate Range Chirp [1/s^2] --> up-chirp
    double tau_p = 4.5511111111111114e-5;       // Chirp duration [s]
    double V = 100.0;                // Effective satellite velocity [m/s]
    double Lambda = 0.03;            // Length of carrier wave [m]
    double R_0 = 1.0e+7;              // Range to center of antenna footprint [m]
    double ta = 1.95;                     // Aperture time [s]
    double prf = 1000.0;// Pulse Repitition Frequency [Hz]
    double c = 3.0e8;
    double f_c = 1.0e+10; //Central frequency
    double alpha = 0.52359877559829893; //squint angle of sight сделать из этого вектор
    double squint = M_PI/2.0;
    //угол скольжения - угол между направлением на Надир и на аппарат
    vector<vector<double>> pos;

    void make_coordinates(size_t size_az){
        double R0_apCenter = 10000.0;
        double Xc = R0_apCenter*cos(alpha);
        //y
        double dy = V/prf;
        double Y_c = Xc*tan(M_PI/2.0-squint);
        vector<double> y = fill_up(-size_az/2, size_az/2, size_az)*dy + vector<double>(size_az,Y_c);
        //z
        double h = 5000.0;
        //pos
        pos = vector<vector<double>>(size_az, vector<double>(3));
        for(size_t i = 0; i < size_az; i++){
            pos[i][0] = Xc;
            pos[i][1] = y[i];
            pos[i][2] = h;
        }
    }

};

double complex_part(std::complex<double> c, string s){
    if(s == "r"){
        return c.real();
    }else if(s == "i"){
        return c.imag();
    }
    return 0.0;
}

vector<std::complex<double>> interp(vector<double>& KY,vector<std::complex<double>>& data,const vector<double>& KY_model, string s){
    size_t size = data.size();
    size_t pos;
    vector<std::complex<double>> temp(size);
    vector<double>::iterator low;
    for(size_t i = 0; i < size;i++){
        low = lower_bound(KY.begin(), KY.end(), KY_model[i]);
        pos = low - KY.begin();
        if(pos < size - 1){
            temp[i] = KY[pos] + KY_model[i]*(complex_part(data[pos+1], s) - complex_part(data[pos+1], s))/(KY[pos+1] - KY[pos]);
        }else{
            temp[i] = data[size-1];
        }
    }
    return temp;
}

template<typename T>
void simple_shift(vector<T>& v){
    vector<T> tmp(v.begin()+(v.size()/2), v.end());
    for(size_t i = 0; i < v.size()/2;i++){
        tmp.push_back(v[i]);
    }
    v = tmp;
}

void shift(vector<vector<std::complex<double>>> &Raw_data){
    size_t az_size = Raw_data.size();
    for(size_t i = 0; i < az_size;i++){
        simple_shift(Raw_data[i]);
    }
    simple_shift(Raw_data);
}

void RVP_correct(vector<vector<std::complex<double>>> &Raw_data, double fs, double K_r, size_t size_range, size_t size_azimuth, double tau_p){
    double c = 3e8;
    double dr = c/(2*tau_p * K_r);
    vector<double> f_r = fill_up(-size_range/2, size_range/2, static_cast<int>(size_range))*(2.0*K_r*dr/c);//range frequency

    //cout << "f_r " << f_r << endl;

    vector<vector<std::complex<double>>> RD(size_azimuth, vector<std::complex<double>>(size_range));
    vector<std::complex<double>> phs_compensation(size_range);
    double norm_range = 1.0/size_range;
    for(size_t i = 0; i < size_range;i++){
        phs_compensation[i] = exp(-M_PI*(f_r[i]*f_r[i]/K_r)*static_cast<complex<double>>(I));
    }
    //cout << "phs_compensation: " << phs_compensation << endl;
    fftw_plan plan_f, plan_b;
    shift(Raw_data);
    for(size_t i = 0; i < size_azimuth;i++) {
        plan_f = fftw_plan_dft_1d(size_range, (fftw_complex *) &Raw_data[i][0],
                                  (fftw_complex *) &RD[i][0], FFTW_FORWARD, FFTW_ESTIMATE); //making draft plan

        fftw_execute(plan_f); // Fourier Transform

    }
    shift(RD);

    for(size_t i = 0; i < size_azimuth;i++) {
        RD[i] = RD[i] * phs_compensation;
    }
    shift(RD);
    for(size_t i = 0; i < size_azimuth;i++) {
        plan_b = fftw_plan_dft_1d(size_range, (fftw_complex*) &RD[i][0],
                                  (fftw_complex*) &Raw_data[i][0], FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_b);
        Raw_data[i] = Raw_data[i] * norm_range;
    }
    shift(Raw_data);
    fftw_destroy_plan(plan_f);
    fftw_destroy_plan(plan_b);
}

void PHS_to_const_ref(vector<vector<std::complex<double>>> &Raw_data, vector<vector<double>>& pos, double tau_p,double  fs, double K_r, double f_c, size_t size_az, double c = 3.0e8){
    vector<double> r(size_az);
    double min_r = 1000000000000.0;

    for(size_t i = 0; i < size_az; i++){
        r[i] = sqrt(pos[i][0]*pos[i][0] + pos[i][1] * pos[i][1] + pos[i][2] * pos[i][2]);
        if(r[i] < min_r){
            min_r = r[i];
        }
    }

    vector<double> DR = r + vector<double>(size_az,-min_r);
    vector<double> tau = fill_up(-tau_p/2, tau_p/2, 1/fs);// time axis in range
    for(size_t i = 0; i < size_az; i++){
        for(size_t j = 0; j < tau.size(); j++){
            Raw_data[i][j] = Raw_data[i][j]*exp(-4*M_PI*K_r *(f_c/K_r + tau[j])* DR[i]*static_cast<complex<double>>(I)/c);//Вставить время по дальности
        }

    }

}
void Azimuth_FFT(vector<vector<std::complex<double>>> &Raw_data, size_t size_range, size_t size_azimuth){
    fftw_plan plan_f;
    vector<complex<double>> temp(size_azimuth);
    for(size_t j = 0; j < size_range;j++){
        for(size_t i = 0; i < size_azimuth;i++){
            temp[i] = Raw_data[i][j];
        }
        plan_f = fftw_plan_dft_1d(size_azimuth, (fftw_complex*) &temp[0],
                                  (fftw_complex*) &temp[0], FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan_f);
        for(size_t i = 0; i < size_azimuth;i++){
            Raw_data[i][j] = temp[i];
        }
    }
    fftw_destroy_plan(plan_f);
}

void Matching_filter(vector<vector<std::complex<double>>>& data, size_t size_azimuth, size_t size_range, Satellite& sat, vector<vector<double>>& KY, double KY_min, double KY_max){
    vector<vector<std::complex<double>>> filter(size_azimuth, vector<std::complex<double>>(size_range));
    vector<double> KR(size_range);
    vector<double> KX(size_azimuth);

    double  r_c, sqrt_KX_KR, KX_max;


    vector<double> t = fill_up(-sat.ta/2, sat.ta/2, 1/sat.prf);// time axis in azimuth
    vector<double> tau = fill_up(-sat.tau_p/2, sat.tau_p/2, 1/sat.fs);// time axis in range

    for(size_t i = 0; i < size_azimuth;i++){
        KX[i] = 4*M_PI*sat.f_c * cos(sat.alpha)/sat.c;
        r_c = 0.0;//must depend on azimuth
        for(size_t j = 0;j < size_range;j++){
            KR[j] = 4*M_PI*sat.K_r/sat.c * ((sat.f_c/sat.K_r) + tau[j] - 2*r_c/sat.c);//частично реализовано в PHS_to_const_ref
            if(KR[j] * KR[j] - KX[i] * KX[i] > 0.0){
                sqrt_KX_KR = sqrt(KR[j] * KR[j] - KX[i] * KX[i]);
            }else{
                sqrt_KX_KR = 0.0;
            }
            KY[i][j] = sqrt_KX_KR;

            if(sqrt_KX_KR < KY_min){
                KY_min = sqrt_KX_KR;
            }
            if(sqrt_KX_KR > KY_max && KX[i] > KX_max){
                KY_max = sqrt_KX_KR;
                KX_max = KX[i];
            }
            data[i][j] *= exp(KR[j] * r_c + r_c * sqrt_KX_KR * static_cast<complex<double>>(I));
        }
    }
}
void RMA(){ //Range migration algorithm or omega-K algorithm
    vector<vector<std::complex<double>>> Raw_data = read_file(false, "/home/ivan/CLionProjects/Omega-K/Example_with_rectangle.txt");
    size_t size_azimuth = Raw_data.size();
    size_t size_range = Raw_data[0].size();

    Satellite new_sat;
    new_sat.make_coordinates(size_azimuth);
    //(1) Compression
    //В sim_demo.py нет стадии свёртки с опорным сигналом, но скорее всего она тут нужна
    //SAR(Raw_data, new_sat.fs, new_sat.K_r, new_sat.tau_p, new_sat.V, new_sat.Lambda,
    //    new_sat.R_0, new_sat.ta, new_sat.prf, size_azimuth, size_range);


    //(2) RVP correction. Checked!
    RVP_correct(Raw_data, new_sat.fs, new_sat.K_r, size_range, size_azimuth, new_sat.tau_p);
    equality(Raw_data, "/home/ivan/CLionProjects/Omega-K/Example_with_rectangle_after_RVP_correct.txt","RVP", "2");
    //(2.1) Const ref
    PHS_to_const_ref(Raw_data, new_sat.pos, new_sat.tau_p,new_sat.fs, new_sat.K_r, new_sat.f_c, size_azimuth);
    equality(Raw_data, "/home/ivan/CLionProjects/Omega-K/Example_with_rectangle_after_constant_ref.txt","Constant_ref", "2");
    //(3) Azimuth FFT
    Azimuth_FFT(Raw_data, size_range, size_azimuth);
    equality(Raw_data, "/home/ivan/CLionProjects/Omega-K/Example_with_rectangle_after_azimuth_fft.txt","Azimuth FFT", "2");
    //(4) Matching filter
    vector<vector<double>> KY(size_azimuth, vector<double>(size_range));
    double KY_min = 0.0, KY_max = 0.0;
    Matching_filter(Raw_data, size_azimuth, size_range, new_sat, KY, KY_min, KY_max);

    //(5) Stolt interpolation
    //double B_eff = K_r * tau_p;
    //double P_o = cos(psi);
    //double theta_p = (static_cast<double>(size_azimuth)/static_cast<double>(size_range))/((f_c/B_eff) - 0.5);
    //double r_one = (4*M_PI/c)*(f_c - B_eff/2.0)*P_o;
    //double Delta_X = 2 * r_one * tan(theta_p/2.0);
    //double KR_max = *std::max_element(KR.begin(), KR.end());

    vector<double> KY_model = fill_up(KY_min, KY_max, size_range);
    fftw_complex *New_phase;
    New_phase = (fftw_complex*) fftw_malloc(size_azimuth*size_range* sizeof(fftw_complex));
    //fftw_complex New_phase[size_range][size_azimuth];
    fftw_complex *Result;
    Result = (fftw_complex*) fftw_malloc(size_azimuth*size_range* sizeof(fftw_complex));
    //vector<vector<std::complex<double>>> New_phase(size_azimuth, vector<std::complex<double>>(size_range, 0.0));
    vector<std::complex<double>> temp(size_range);
    for(size_t i = 0;i < size_azimuth; i++){
        temp = interp(KY[i], Raw_data[i], KY_model,  "r") +
                interp(KY[i], Raw_data[i], KY_model, "i") * static_cast<std::complex<double>>(I);
        for(size_t j = 0; j < size_range;j++){
            New_phase[j+size_range*i][0] = temp[j].real();
            New_phase[j+size_range*i][1] = temp[j].imag();
        }
    }
    //(6) 2D-IFFT

    fftw_plan plan = fftw_plan_dft_2d(size_azimuth,size_range,
                                      New_phase, Result,
                                    FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    Write_in_file(Raw_data, "Test_file_for_omega_k");
}

int main() {
    RMA();
    return 0;
}
