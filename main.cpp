//#include <iostream>
//#include <fstream>
#include <vector>
#include <complex>
//#include <iomanip>
#include <cmath>
#include <string>
#include <fftw3.h>
#include <algorithm>
#include "extra_tools.h"
#include "compression_algorithm.h"
#include "tester.h"
#include "Time_duration_file.h"

struct Satellite{
    // Sensor parameters (ERS satellite)
    //DIRSIG
    double fs = 5.0e+7;  //Range Sampling Frequency [Hz]
    double K_r = 15956460791541.389;     // FM Rate Range Chirp [1/s^2] --> up-chirp
    double tau_p = 1.0e-5;       // Chirp duration [s]
    double V = 200.0;                // Effective satellite velocity [m/s]
    double Lambda = 0.02;            // Length of carrier wave [m]

    double ta = 1.95;                     // Aperture time [s]
    double prf = 400.0;// Pulse Repitition Frequency [Hz]
    double c = 3.0e8;
    double f_c = 1.5e+10; //Central frequency
    double alpha = M_PI/4.0; //squint angle of sight сделать из этого вектор
    double squint = M_PI/2.0;
    double L;
    //угол скольжения - угол между направлением на Надир и на аппарат
    vector<vector<double>> pos;
    vector<double> r;
    double R_0 = 15000/cos(alpha);              // Range to center of antenna footprint [m]
    double min_r = 15000/cos(alpha);
    /*
    double fs = 4.5e+7;  //Range Sampling Frequency [Hz]
    double K_r = 4.8e+12;     // FM Rate Range Chirp [1/s^2] --> up-chirp
    double tau_p = 4.5511111111111114e-5;       // Chirp duration [s]
    double V = 100.0;                // Effective satellite velocity [m/s]
    //double Lambda = 0.03;            // Length of carrier wave [m]
    //double R_0 = 1.0e+7;              // Range to center of antenna footprint [m]
    double ta = 1.95;                     // Aperture time [s]
    double prf = 1000.0;// Pulse Repitition Frequency [Hz]
    double c = 3.0e8;
    double f_c = 1.0e+10; //Central frequency
    double alpha = 0.52359877559829893; //squint angle of sight сделать из этого вектор
    double squint = M_PI/2.0;
    double L;
    //угол скольжения - угол между направлением на Надир и на аппарат
    vector<vector<double>> pos;
    vector<double> r;
    double min_r = 1000000000000.0;*/
    size_t size_range;
    size_t size_azimuth;
    Satellite(size_t size_az, size_t size_r):size_azimuth(size_az), size_range(size_r){
        r = vector<double>(size_azimuth);
    }

    void make_coordinates(){
        double R0_apCenter = 10000.0;
        double Xc = R0_apCenter*cos(alpha);
        //y
        double dy = V/prf;
        double Y_c = Xc*tan(M_PI/2.0-squint);
        vector<double> y = fill_up(-static_cast<int>(size_azimuth/2), static_cast<int>(size_azimuth/2), static_cast<int>(size_azimuth))*dy ;
        y = y + vector<double>(size_azimuth,Y_c);
        //z
        double h = 5000.0;
        //pos
        pos = vector<vector<double>>(size_azimuth, vector<double>(3));
        for(size_t i = 0; i < size_azimuth; i++){
            pos[i][0] = Xc;
            pos[i][1] = y[i];
            pos[i][2] = h;
        }
        L = sqrt((pos[0][0] - pos[size_azimuth-1][0])*(pos[0][0] - pos[size_azimuth-1][0]) +
                         (pos[0][1] - pos[size_azimuth-1][1])*(pos[0][1] - pos[size_azimuth-1][1]) +
                         (pos[0][2] - pos[size_azimuth-1][2])*(pos[0][2] - pos[size_azimuth-1][2]));
        for(size_t i = 0; i < size_azimuth; i++){
            r[i] = sqrt(pos[i][0]*pos[i][0] + pos[i][1] * pos[i][1] + pos[i][2] * pos[i][2]);
            if(r[i] < min_r){
                min_r = r[i];
            }
        }
    }

};


void RVP_correct(vector<vector<std::complex<double>>> &Raw_data, Satellite& sat){
    double c = 299792458.0;
    double dr = c/(2*sat.tau_p * sat.K_r);
    vector<double> f_r = fill_up(-sat.size_range/2, sat.size_range/2, static_cast<int>(sat.size_range))*(2.0*sat.K_r*dr/c);//range frequency

    vector<vector<std::complex<double>>> RD(sat.size_azimuth, vector<std::complex<double>>(sat.size_range));
    vector<std::complex<double>> phs_compensation(sat.size_range);
    double norm_range = 1.0/sat.size_range;
    for(size_t i = 0; i < sat.size_range;i++){
        phs_compensation[i] = exp(-M_PI*(f_r[i]*f_r[i]/sat.K_r)*static_cast<complex<double>>(I));
    }

    fftw_plan plan_f, plan_b;
    shift(Raw_data);
    for(size_t i = 0; i < sat.size_azimuth;i++) {
        plan_f = fftw_plan_dft_1d(sat.size_range, (fftw_complex *) &Raw_data[i][0],
                                  (fftw_complex *) &RD[i][0], FFTW_FORWARD, FFTW_ESTIMATE); //making draft plan

        fftw_execute(plan_f); // Fourier Transform

    }
    shift(RD);

    for(size_t i = 0; i < sat.size_azimuth;i++) {
        RD[i] = RD[i] * phs_compensation;
    }
    shift(RD);
    for(size_t i = 0; i < sat.size_azimuth;i++) {
        plan_b = fftw_plan_dft_1d(sat.size_range, (fftw_complex*) &RD[i][0],
                                  (fftw_complex*) &Raw_data[i][0], FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_b);
        Raw_data[i] = Raw_data[i] * norm_range;
    }
    shift(Raw_data);
    fftw_destroy_plan(plan_f);
    fftw_destroy_plan(plan_b);
}

void PHS_to_const_ref(vector<vector<std::complex<double>>> &Raw_data, Satellite& sat){

    vector<double> DR = sat.r + vector<double>(sat.size_azimuth,-sat.min_r);
    vector<double> tau = fill_up(-sat.tau_p/2.0, sat.tau_p/2.0, 1.0/sat.fs);// time axis in range
    vector<vector<std::complex<double>>> e(sat.size_azimuth, vector<complex<double>>(tau.size()));

    for(size_t i = 0; i < sat.size_azimuth; i++){
        for(size_t j = 0; j < tau.size(); j++){
            e[i][j] = exp(-4.0*M_PI*sat.K_r *(sat.f_c/sat.K_r + tau[j])* DR[i]*static_cast<complex<double>>(I)/sat.c);
        }
    }
    for(size_t i = 0; i < sat.size_azimuth; i++){
        for(size_t j = 0; j < tau.size(); j++){
            Raw_data[i][j] = Raw_data[i][j]*e[i][j];
        }
    }
}

void Azimuth_FFT(vector<vector<std::complex<double>>> &Raw_data, size_t size_range, size_t size_azimuth){
    fftw_plan plan_f;
    vector<complex<double>> temp(size_azimuth);

    shift(Raw_data);

    for(size_t j = 0; j < size_range;j++){
        for(size_t i = 0; i < size_azimuth;i++){
            temp[i] = Raw_data[i][j];
        }
        //plan_f = fftw_plan_dft_1d(size_azimuth, (fftw_complex*) &temp[0],
        //                          (fftw_complex*) &temp[0], FFTW_FORWARD, FFTW_ESTIMATE); //было
        plan_f = fftw_plan_dft_1d(size_azimuth, (fftw_complex*) &temp[0],
                                  (fftw_complex*) &temp[0], FFTW_BACKWARD, FFTW_ESTIMATE); //сделано на основе RITSAR
        fftw_execute(plan_f);
        for(size_t i = 0; i < size_azimuth;i++){
            Raw_data[i][j] = complex<double>(temp[i].real()/size_azimuth, temp[i].imag()/size_azimuth);//https://habr.com/ru/company/otus/blog/449996/
        }
    }
    shift(Raw_data);
    fftw_destroy_plan(plan_f);
}

void Matching_filter(vector<vector<std::complex<double>>>& data, Satellite& sat, vector<double>& KR, vector<double>& KX, double& KY_min, double& KY_max){
    vector<vector<std::complex<double>>> filter(sat.size_azimuth, vector<std::complex<double>>(sat.size_range));
    KR = vector<double>(sat.size_range);


    double phase_mf;
    double KY;
    vector<double> t = fill_up(-sat.ta/2, sat.ta/2, 1/sat.prf);// time axis in azimuth
    vector<double> tau = fill_up(-static_cast<int>(sat.size_range/2), static_cast<int>(sat.size_range/2), static_cast<int>(sat.size_range)) * (1.0/sat.fs);// time axis in range

    KX = fill_up(-static_cast<int>(sat.size_azimuth/2), static_cast<int>(sat.size_azimuth/2), static_cast<int>(sat.size_azimuth)) * (2.0 * M_PI/sat.L);
    for(size_t i = 0; i < sat.size_azimuth;i++){
        //r_c = 0.0;//must depend on azimuth
        for(size_t j = 0;j < sat.size_range;j++){
            KR[j] = 4*M_PI/sat.c * (sat.f_c + sat.K_r*tau[j]);//частично реализовано в PHS_to_const_ref
            KY = sqrt(KR[j]*KR[j] - KX[i]*KX[i]);
            if(i == 0 && j == 0){
                KY_min = KY;
            }
            if(KY > KY_max){
                KY_max = KY;
            }
            if(KY < KY_min){
                KY_min = KY;
            }
            phase_mf = -sat.min_r*KR[j] + sat.min_r * KY;
            data[i][j] = data[i][j]*exp(phase_mf*static_cast<complex<double>>(I));
        }
    }
}

void Stolt_interpolation(vector<vector<std::complex<double>>>& data,Satellite& sat, double KY_min, double KY_max, vector<double>& KR, vector<double>& KX, int tay = 17, int upsample = 6){
    vector<double> KY = fill_up(KY_min, KY_max, sat.size_range);
    vector<vector<complex<double>>> S(sat.size_azimuth, vector<complex<double>>(sat.size_range, 0.0));
    vector<vector<complex<double>>> S_new(sat.size_azimuth-1, vector<complex<double>>(sat.size_range-1, 0.0));
    vector<double> KY_temp(sat.size_range);
    vector<double> data_real(sat.size_range);
    vector<double> data_imag(sat.size_range);
    vector<double> S_real(sat.size_range);
    vector<double> S_imag(sat.size_range);

    //KY_temp = K_x
    //KX = K_y
    //KY = K_xi


    for(size_t i = 0; i < sat.size_azimuth; i++){
        fill_up_KY(KY_temp, KR, KX[i]);
        data_real = real(data[i]);
        data_imag = imag(data[i]);
        interp1d( KY_temp,data_real, KY, S_real);
        interp1d( KY_temp,data_imag, KY, S_imag);
        for(size_t j = 0; j < sat.size_range; j++){
            S[i][j] += complex<double>(S_real[j], S_imag[j]);

        }

    }

    for(size_t i = 0; i < sat.size_azimuth-1;i++){
        for(size_t j = 0; j < sat.size_range-1;j++){
            S_new[i][j] = S[i][j];
        }
    }

    vector<double> win_x = taylor(S_new[0].size(), tay);
    vector<int> new_v_x = {(int)S_new.size(), 1};
    vector<vector<double>> win_x_new = tile_x(win_x, new_v_x);

    vector<double> win_y = taylor(S_new.size(), tay);
    vector<int> new_v_y = {1, (int)S_new[0].size()};
    vector<vector<double>> win_y_new = tile_y(win_y, new_v_y);


    vector<vector<double>> win = win_x_new*win_y_new;


    vector<vector<complex<double>>> w_c(win.size(), vector<complex<double>>(win[0].size()));

    for(size_t i = 0;i < win.size();i++){
        for(size_t j = 0;j < win[0].size();j++){
            w_c[i][j] = complex<double>(win[i][j]);
        }

    }

    vector<vector<complex<double>>> S_win = S_new*w_c;

    int length = pow(2, (int)(log2(S_new.size()*upsample))+1);
    int pad_x = length-(int)S_win[0].size();
    int pad_y = length-(int)S_win.size();

    vector<vector<complex<double>>> S_pad = pad_constant(S_win,{{pad_y/2, pad_y/2},{pad_x/2,pad_x/2}});

    data = S_pad;
}

void IFFT_2D(vector<vector<std::complex<double>>>& Raw_data){
    size_t size_azimuth = Raw_data.size();
    size_t size_range = Raw_data[0].size();
    shift(Raw_data);
    fftw_complex *New_phase;
    New_phase = (fftw_complex*) fftw_malloc(Raw_data.size()*Raw_data[0].size()* sizeof(fftw_complex));
    //fftw_complex New_phase[size_range][size_azimuth];
    fftw_complex *Result;
    Result = (fftw_complex*) fftw_malloc(Raw_data.size()*Raw_data[0].size()* sizeof(fftw_complex));
    //vector<vector<std::complex<double>>> New_phase(size_azimuth, vector<std::complex<double>>(size_range, 0.0));
    for(size_t i = 0;i < size_azimuth; i++){
        for(size_t j = 0; j < size_range;j++){
            New_phase[j+size_range*i][0] = Raw_data[i][j].real();
            New_phase[j+size_range*i][1] = Raw_data[i][j].imag();
        }
    }


    fftw_plan plan = fftw_plan_dft_2d((int)size_azimuth,(int)size_range,
                                      New_phase, Result,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);

    for(size_t i = 0;i < size_azimuth; i++){
        for(size_t j = 0; j < size_range;j++){
            Raw_data[i][j] = std::complex<double>(Result[j+size_range*i][0],Result[j+size_range*i][1]);
        }
    }

    fftw_destroy_plan(plan);
    shift(Raw_data);
    //equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_2D_FFT.txt","2D FFT", "2");
}
void RMA(){ //Range migration algorithm or omega-K algorithm
    vector<vector<std::complex<double>>> Raw_data = read_file(false, "/home/ivanemekeev/CLionProjects/SAR-data/DIRSIG_input.txt");
    size_t size_azimuth = Raw_data.size();
    size_t size_range = Raw_data[0].size();

    Satellite new_sat(size_azimuth, size_range);
    new_sat.make_coordinates();
    //(1) Compression

    //В sim_demo.py нет стадии свёртки с опорным сигналом, но скорее всего она тут не нужна
    //SAR(Raw_data, new_sat.fs, new_sat.K_r, new_sat.tau_p, new_sat.V, new_sat.Lambda,
    //    new_sat.R_0, new_sat.ta, new_sat.prf, size_azimuth, size_range);


    //(2) RVP correction.
    RVP_correct(Raw_data, new_sat);
    //simple_equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_RVP_correct.txt","RVP", "2");

    //(2.1) Const ref
    PHS_to_const_ref(Raw_data, new_sat);
    //simple_equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_constant_ref.txt","Constant_ref", "2");

    //(3) Azimuth FFT
    Azimuth_FFT(Raw_data, size_range, size_azimuth);
    //simple_equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_azimuth_fft.txt","Azimuth FFT", "2");

    //(4) Matching filter
    vector<vector<double>> KY(size_azimuth, vector<double>(size_range));
    double KY_min = 0.0, KY_max = 0.0;
    vector<double> KX, KR;
    Matching_filter(Raw_data, new_sat, KR, KX, KY_min, KY_max);
    //simple_equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_matching_filter.txt","Matching filter", "2");

    //(5) Stolt interpolation
    Stolt_interpolation(Raw_data, new_sat, KY_min, KY_max, KR, KX, 13,  2);
    //simple_equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_before_2D_FFT.txt","Stolt interpolation", "2");

    //(6) 2D-IFFT
    IFFT_2D(Raw_data);
    Write_in_file(Raw_data, "DIRSIG_output");
}

int main() {
    {
        TimeDuration t("RMA");
        RMA();
    }
    return 0;
}
