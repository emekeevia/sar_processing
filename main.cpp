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
#include "CubicSlpine.h"
//#include "compression_algorithm.h"
#include "tester.h"
#include "Time_duration_file.h"

struct Satellite{
    // Sensor parameters (ERS satellite)
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


void RVP_correct(vector<vector<std::complex<double>>> &Raw_data, Satellite& sat){
    double c = 3e8;
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

void fill_up_KY(vector<double>& KY, vector<double>& KR,double KX_i, Satellite& sat){
    double temp;
    for(size_t i = 0; i < KR.size();i++){
        temp = KR[i] * KR[i] - (KX_i * KX_i);
        KY[i] = sqrt(temp);
    }
}


std::vector<int> argsort(std::vector<double>& v){
    std::vector<int> temp(v.size(), 0);
    for(int i = 0; i < v.size();i++){
        temp[i] = i;
    }
    quick_sort(v, temp,0, v.size() - 1);
    return temp;
}


void interp1d(std::vector<double>& x,std::vector<double>& y, std::vector<double>& x_new, std::vector<double>& y_new,  string kind = "linear"){

    std::vector<int> ind = argsort(x);
    std::vector<double> data_sorted(x.size());

    for(int i = 0;i < y.size();i++){
        data_sorted[i] = y[ind[i]];
    }

    std::vector<bool> entry(x_new.size(), true);
    std::vector<int> x_new_indices(x_new.size());
    for(int i = 0;i < x_new.size();i++){
        size_t start = lower_bound(x.begin(), x.end(), x_new[i]) - x.begin();
        if(x_new[i] < x[0]) {
            entry[i] = false;
        }
        if(x_new[i] == x[start]){
            start++;
        }
        if(start > x.size() - 1){
            start = x.size() - 1;
            entry[i] = false;
        }

        x_new_indices[i]=(int)start;
    }

    std::vector<int> lo = x_new_indices - 1;
    std::vector<int> hi = x_new_indices;

    double x_lo , x_hi , y_lo , y_hi ;
    double slope ;


    if(kind == "linear"){
        for(size_t i = 0; i < lo.size();i++){
            x_lo = x[lo[i]];
            x_hi = x[hi[i]];
            y_lo = data_sorted[lo[i]];
            y_hi = data_sorted[hi[i]];

            if(x_hi == x_lo){
                y_new[i] = 0.0;
            }else{
                if(entry[i]){
                    slope = (y_hi - y_lo)/(x_hi - x_lo);
                    y_new[i] = (slope * (x_new[i] - x_lo)) + y_lo;
                }else{
                    y_new[i] = 0.0;
                }

            }

        }
    }else if(kind == "kubic"){
        CubicSpline spline(x, data_sorted);
        for(size_t i = 0; i < lo.size();i++){
            y_new[i] = spline.interpolate(x_new[i]);
        }
    }

}




vector<double> taylor(size_t nsamples, size_t S_L=43){
    vector<double> xi = fill_up(-0.5, 0.5, nsamples);

    double A = (1.0/M_PI) * acosh(pow(10.0, ((double)S_L*1.0)/20.0));
    int n_bar = (int)(2*A*A + 0.5) + 1;
    double sigma_p = (double)n_bar/sqrt(A*A+((double)n_bar-0.5)*((double)n_bar-0.5));
    vector<int> m = fill_up(1, n_bar);
    vector<int> n = fill_up(1,n_bar);
    vector<double> F_m(n_bar-1, 0.0);
    for(auto i: m){
        double num = 1.0;
        double den = 1.0;
        for(auto j: n){
            num = num*pow(-1, i+1)*(1-i*i*1.0/(sigma_p*sigma_p)/(A*A+(j-0.5)*(j-0.5)));
            if(i!=j){
                den = den*(1-i*i*1.0/(j*j));
            }
        }
        F_m[i-1] = num/den;
    }


    vector<double> w(nsamples, 1.0);
    for(auto i : m) {
        for(size_t j = 0;j < nsamples;j++){
            w[j] += F_m[i - 1] * cos(2 * M_PI * i * xi[j]);
        }
    }
    w = w/(*max_element(w.begin(), w.end()));

    return w;

}

template<typename T>
vector<vector<T>> tile_x(vector<T>& A, vector<int> reps){
    vector<vector<T>> temp(reps[0], vector<T>(reps[1]*A.size(), 0.0));
    for(size_t i = 0; i < reps[0];i++){
        for(size_t j = 0; j < reps[1]*A.size();j++) {
            temp[i][j] = A[j % A.size()];
        }
    }
    return temp;
}
template<typename T>
vector<vector<T>> tile_y(vector<T>& A, vector<int> reps){
    vector<vector<T>> temp(reps[0]*A.size(), vector<T>(reps[1], 0.0));
    for(size_t i = 0; i < reps[0]*A.size();i++){
        for(size_t j = 0; j < reps[1];j++) {
            temp[i][j] = A[i % A.size()];
        }
    }
    return temp;
}




vector<vector<complex<double>>> pad_constant(vector<vector<std::complex<double>>>& array,vector<vector<int>> pad_width){

    vector<vector<complex<double>>> new_array(array.size() +pad_width[0][0] + pad_width[0][1], vector<complex<double>>(array[0].size() +pad_width[1][0] + pad_width[1][1], 0.0));
    for(size_t i = 0; i < array.size();i++){
        for(size_t j = 0; j < array[0].size() ;j++){
            new_array[i+pad_width[0][0]][j+pad_width[1][0]] = array[i][j];
        }
    }
    //equality(new_array, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_expansion.txt","Stolt interpolation part 3 after expansion", "2");

    vector<vector<double>> values(2, vector<double>(2, 0.0));
    for(size_t k = 0; k < 2;k++){
        vector<vector<complex<double>>> roi;
        if(k == 0){
            for(size_t i = 0;i < new_array.size();i++){
                vector<complex<double>> temp_temp;
                for(size_t j = pad_width[1][0]; j< new_array.size()-pad_width[1][0]; j++){
                    temp_temp.push_back(new_array[i][j]);
                }
                roi.push_back(temp_temp);
            }
        }else{
            roi = new_array;
        }

        if(k == 0){
            //left_slice
            vector<int> rows_l = {0, pad_width[k][0]};
            vector<int> cols_l = {0, (int)roi[0].size()};
            set_array(roi, values[k], rows_l, cols_l);

            //right_slice
            vector<int> rows_r = {(int)roi.size()-pad_width[k][1], (int)roi.size()};
            vector<int> cols_r = {0, (int)roi[0].size()};
            set_array(roi, values[k], rows_r, cols_r);
        }else{
            //left_slice
            vector<int> rows_l = {0, (int)roi.size()};
            vector<int> cols_l = {0, pad_width[k][0]};

            set_array(roi, values[k], rows_l, cols_l);
            //right_slice
            vector<int> rows_r = {0, (int)roi.size()};
            vector<int> cols_r = {(int)roi[0].size() - pad_width[k][1], (int)roi[0].size()};

            set_array(roi, values[k], rows_r, cols_r);

        }


        //k = axis
        //pad_width[k] = width_pair
        //values[k] = value_pair

    }
    return new_array;
}

void Stolt_interpolation(vector<vector<std::complex<double>>>& data,Satellite& sat, double KY_min, double KY_max, vector<double>& KR, vector<double>& KX, int upsample = 6){
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
        fill_up_KY(KY_temp, KR, KX[i], sat);
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

    vector<double> win_x = taylor(S_new[0].size(), 17);
    vector<int> new_v_x = {(int)S_new.size(), 1};
    vector<vector<double>> win_x_new = tile_x(win_x, new_v_x);

    vector<double> win_y = taylor(S_new.size(), 17);
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
void RMA(){ //Range migration algorithm or omega-K algorithm
    vector<vector<std::complex<double>>> Raw_data = read_file(false, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle.txt");
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
    simple_equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_RVP_correct.txt","RVP", "2");
    //(2.1) Const ref
    PHS_to_const_ref(Raw_data, new_sat);
    simple_equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_constant_ref.txt","Constant_ref", "2");
    //(3) Azimuth FFT
    Azimuth_FFT(Raw_data, size_range, size_azimuth);
    simple_equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_azimuth_fft.txt","Azimuth FFT", "2");
    //(4) Matching filter
    vector<vector<double>> KY(size_azimuth, vector<double>(size_range));
    double KY_min = 0.0, KY_max = 0.0;
    vector<double> KX, KR;
    Matching_filter(Raw_data, new_sat, KR, KX, KY_min, KY_max);
    simple_equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_matching_filter.txt","Matching filter", "2");

    //(5) Stolt interpolation

    //Raw_data = read_file(false, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_matching_filter.txt");

    Stolt_interpolation(Raw_data, new_sat, KY_min, KY_max, KR, KX, 2);
    simple_equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_before_2D_FFT.txt","Stolt interpolation", "2");
    //(6) 2D-IFFT
    Raw_data = read_file(false, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_before_2D_FFT.txt");
    size_azimuth = Raw_data.size();
    size_range = Raw_data[0].size();
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
    simple_equality(Raw_data, "/home/ivanemekeev/CLionProjects/SAR-data/Example_with_rectangle_after_2D_FFT.txt","2D FFT", "2");
    //Write_in_file(Raw_data, "Test_file_for_omega_k");
}

int main() {
    {
        TimeDuration t("RMA");
        RMA();
    }
    return 0;
}
