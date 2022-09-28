#include "tester.h"

double relative_er(double one, double two){
    if(one > 5.0e-4 && two > 5.0e-4){
        return abs(one - two)/max(abs(one), abs(two));
    }else {
        return abs(one - two);
    }
}

double metrik_inf(complex<double> my,complex<double> example){
    double my_r = my.real();
    double my_i = my.imag();
    double example_r = example.real();
    double example_i = example.imag();
    return max(relative_er(my_r, example_r), relative_er(my_i, example_i));
}
double metrik_2(complex<double> my,complex<double> example){
    double my_abs = abs(my);
    double example_abs = abs(example);
    if(max(my_abs, example_abs) < 1.0e-5){
        return abs(my_abs - example_abs);
    }
    return abs(my_abs - example_abs)/max(my_abs, example_abs);
}

void equality(vector<vector<complex<double>>>& in_mas, string file_with_comp_data, string step_name, string metrik){
    vector<vector<complex<double>>> sample_mas = read_file(false, file_with_comp_data);
    size_t azimuth_size = in_mas.size();
    size_t range_size = in_mas[0].size();
    double epsilon = 0.015;
    double in_real,in_imag, sample_real, sample_imag;
    bool flag = false;
    double max_rel_error = 0.0;
    double rel_error = 0.0;
    int x, y;

    for(size_t i = 0; i < azimuth_size;i++){
        for(size_t j = 0; j < range_size;j++){
            if(metrik == "2"){
                rel_error = metrik_2(in_mas[i][j], sample_mas[i][j]);
            }else if(metrik == "inf"){
                rel_error = metrik_inf(in_mas[i][j], sample_mas[i][j]);
            }

            if(rel_error>= epsilon){
                flag = true;
                if(rel_error > max_rel_error){
                    max_rel_error = rel_error;
                    x = i;
                    y = j;
                }
            }
        }
    }
    if(flag){
        cerr << step_name << ": error " << max_rel_error << "\n";
        cerr << x << " " << y << "\n";
        cerr << "My: " << in_mas[x][y] << "\n";
        cerr << "Sample: " << sample_mas[x][y] << "\n";
    }else{
        cerr << step_name << ": OK" << "\n";
    }
}