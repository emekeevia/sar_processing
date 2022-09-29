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
    double L;
    //угол скольжения - угол между направлением на Надир и на аппарат
    vector<vector<double>> pos;
    vector<double> r;
    double min_r = 1000000000000.0;
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

void RVP_correct(vector<vector<std::complex<double>>> &Raw_data, Satellite& sat){
    double c = 3e8;
    double dr = c/(2*sat.tau_p * sat.K_r);
    vector<double> f_r = fill_up(-sat.size_range/2, sat.size_range/2, static_cast<int>(sat.size_range))*(2.0*sat.K_r*dr/c);//range frequency

    //cout << "f_r " << f_r << endl;

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

void Matching_filter(vector<vector<std::complex<double>>>& data, Satellite& sat, vector<double>& KR, vector<double>& KX){
    vector<vector<std::complex<double>>> filter(sat.size_azimuth, vector<std::complex<double>>(sat.size_range));
    KR = vector<double>(sat.size_range);


    double phase_mf;

    vector<double> t = fill_up(-sat.ta/2, sat.ta/2, 1/sat.prf);// time axis in azimuth
    vector<double> tau = fill_up(-sat.tau_p/2, sat.tau_p/2, 1/sat.fs);// time axis in range
    KX = fill_up(-static_cast<int>(sat.size_azimuth/2), static_cast<int>(sat.size_azimuth/2), static_cast<int>(sat.size_azimuth)) * (2.0 * M_PI/sat.L);
    for(size_t i = 0; i < sat.size_azimuth;i++){
        //r_c = 0.0;//must depend on azimuth
        for(size_t j = 0;j < sat.size_range;j++){
            KR[j] = 4*M_PI/sat.c * (sat.f_c + sat.K_r*tau[j]);//частично реализовано в PHS_to_const_ref
            phase_mf = -sat.min_r*KR[j] + sat.min_r * sqrt(KR[j]*KR[j] - KX[i]*KX[i]);
            data[i][j] = data[i][j]*exp(phase_mf*static_cast<complex<double>>(I));
        }
    }
}

void fill_up_KY(vector<double>& KY, vector<double>& KR,double KX_i, Satellite& sat){
    double temp;
    for(size_t i = 0; i < sat.size_range;i++){
        temp = KR[i] * KR[i] - KX_i*KX_i;
        if(temp >= 0){
            KY[i] = sqrt(temp);
        }else{
            KY[i] = 0.0;
        }
    }
}

void interp1d(vector<std::complex<double>>& data,vector<double>& KY_temp, vector<double>& KY, vector<complex<double>>& w, Satellite& sat){

    //(1) смотрим есть ли элементы из KY в KY_temp
    double KY_temp_max = *std::max_element(KY_temp.begin(), KY_temp.end());
    double KY_temp_min = *std::min_element(KY_temp.begin(), KY_temp.end());
    size_t start = lower_bound(KY.begin(), KY.end(), KY_temp_min) - KY.begin();
    size_t finish = upper_bound(KY.begin(), KY.end(), KY_temp_max) - KY.begin();
    if(finish-start > 0){
        //(2) если элементы есть,то находим м/у какими точками KY_temp они лежат
        for(size_t i = start; i < finish;i++){
            size_t left = lower_bound(KY_temp.begin(), KY_temp.end(), KY[i]) - KY_temp.begin();//индекс, который указывает на первый элемент меньше заданного
            complex<double> k = (data[left+1] - data[left])/(KY_temp[left+1] - KY_temp[left]);
            w = w + vector<complex<double>>(sat.size_range,(data[left]+(KY[i] - KY_temp[left])*k));
        }
    }else if(finish == start){

        size_t left = lower_bound(KY_temp.begin(), KY_temp.end(), KY[start]) - KY_temp.begin();//индекс, который указывает на первый элемент меньше заданного
        complex<double> k = (data[left+1] - data[left])/(KY_temp[left+1] - KY_temp[left]);
        w = w + vector<complex<double>>(sat.size_range,(data[left]+(KY[start] - KY_temp[left])*k));
    }

}


void Stolt_interpolation(vector<vector<std::complex<double>>>& data,Satellite& sat, double KY_min, double KY_max, vector<double>& KR, vector<double>& KX){
    vector<double> KY = fill_up(KY_min, KY_max, sat.size_range);
    vector<vector<complex<double>>> S(sat.size_azimuth, vector<complex<double>>(sat.size_range, 0.0));
    vector<double> KY_temp;
    for(size_t i = 0; i < sat.size_azimuth; i++){
        fill_up_KY(KY_temp, KR, KX[i], sat);
        interp1d(data[i], KY_temp, KY, S[i], sat);
    }
    equality(S, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_Stolt_interpolation_first_part.txt","Stolt interpolation", "2");
    /*
      S = np.nan_to_num(S)
    [p1,p2] = phs_inscribe(np.abs(S))
    S_new = S[p1[1]:p2[1],
              p1[0]:p2[0]]

    #Create window
    win_x = sig.taylor(S_new.shape[1],taylor)
    win_x = np.tile(win_x, [S_new.shape[0],1])

    win_y = sig.taylor(S_new.shape[0],taylor)
    win_y = np.array([win_y]).T
    win_y = np.tile(win_y, [1,S_new.shape[1]])

    win = win_x*win_y

    #Apply window
    S_win = S_new*win

    #Pad Spectrum
    length = 2**(int(np.log2(S_new.shape[0]*upsample))+1)
    pad_x = length-S_win.shape[1]
    pad_y = length-S_win.shape[0]
    S_pad = np.pad(S_win,((pad_y//2, pad_y//2),(pad_x//2,pad_x//2)), mode = 'constant')
     */


}
void RMA(){ //Range migration algorithm or omega-K algorithm
    vector<vector<std::complex<double>>> Raw_data = read_file(false, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle.txt");
    size_t size_azimuth = Raw_data.size();
    size_t size_range = Raw_data[0].size();

    Satellite new_sat(size_azimuth, size_range);
    new_sat.make_coordinates();
    //(1) Compression
    //В sim_demo.py нет стадии свёртки с опорным сигналом, но скорее всего она тут нужна
    //SAR(Raw_data, new_sat.fs, new_sat.K_r, new_sat.tau_p, new_sat.V, new_sat.Lambda,
    //    new_sat.R_0, new_sat.ta, new_sat.prf, size_azimuth, size_range);


    //(2) RVP correction. Checked!
    RVP_correct(Raw_data, new_sat);
    equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_RVP_correct.txt","RVP", "2");
    //(2.1) Const ref
    PHS_to_const_ref(Raw_data, new_sat);
    equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_constant_ref.txt","Constant_ref", "2");
    //(3) Azimuth FFT
    Azimuth_FFT(Raw_data, size_range, size_azimuth);
    equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_azimuth_fft.txt","Azimuth FFT", "2");
    //(4) Matching filter
    vector<vector<double>> KY(size_azimuth, vector<double>(size_range));
    double KY_min = 0.0, KY_max = 0.0;
    vector<double> KX, KR;
    Matching_filter(Raw_data, new_sat, KR, KX);
    equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_matching_filter.txt","Matching filter", "2");

    //(5) Stolt interpolation
    Stolt_interpolation(Raw_data, new_sat, KY_min, KY_max, KR, KX);
    //equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_Stolt_interpolation_last_part.txt","Stolt interpolation", "2");

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
