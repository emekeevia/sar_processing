#include "extra_tools.h"


vector<double> fill_up(const double& start,const double& finish, const double& step){
    vector<double> mass;
    for(double temp = start; temp <= finish;temp+=step){
        mass.push_back(temp);
    }
    return mass;
}

vector<double> real(vector<std::complex<double>>& data){
    vector<double> temp(data.size());
    for(size_t i = 0; i < data.size();i++){
        temp[i] = data[i].real();
    }
    return temp;
}

vector<double> imag(vector<std::complex<double>>& data){
    vector<double> temp(data.size());
    for(size_t i = 0; i < data.size();i++){
        temp[i] = data[i].imag();
    }
    return temp;
}

vector<int> fill_up(const int& start, const int& finish){
    vector<int> mass;
    for(double temp = start; temp < finish;temp++){
        mass.push_back(temp);
    }
    return mass;
}

vector<double> fill_up(const double& start,const double& finish, size_t N){
    vector<double> mass(N);
    double step = (finish - start)/(N-1);
    for(size_t i = 0; i < N;i++){
        mass[i] = start + i*step;
    }
    return mass;
}
vector<double> fill_up(const int& start, const int& finish, int N){
    vector<double> mass;
    double step = (static_cast<double>(finish) - static_cast<double>(start) + 1.0)/static_cast<double>(N);
    for(int i = 0; i < N;i++){
        mass.push_back(static_cast<double>(start + i*step));
    }
    return mass;
}

void conjugate(vector<complex<double>>& sopr, const vector<complex<double>>& orig){
    size_t size = orig.size();
    for(size_t i = 0; i < size;i++){
        sopr[i] = conj(orig[i]);
    }
}

vector<complex<double>> operator*(const vector<complex<double>>& v1, const vector<complex<double>>& v2){
    size_t size = v1.size();
    if(size != v2.size()){
        throw("Not equal size!");
    }
    vector<complex<double>> temp(size);
    for(size_t i = 0; i < size;i++){
        temp[i] = v1[i]*v2[i];
    }
    return temp;
}

void shift(vector<vector<std::complex<double>>> &Raw_data){
    size_t az_size = Raw_data.size();
    for(size_t i = 0; i < az_size;i++){
        iFFTshift(Raw_data[i]);
    }
    iFFTshift(Raw_data);
}

size_t partition(std::vector<double>& v, std::vector<int>& temp, size_t l, size_t r){
    double pivot = v[(l + r) / 2];
    size_t i = l;
    size_t j = r;
    while (i <= j){
        while (v[i] < pivot){
            i++;
        }
        while (v[j] > pivot){
            j--;
        }
        if (i >= j){
            break;
        }

        swap(v[i], v[j]);
        swap(temp[i++], temp[j--]);
    }

    return j;
}


void quick_sort(std::vector<double>& v, std::vector<int>& temp, size_t l, size_t r){
    if(l < r){
        size_t q = partition(v,temp, l, r);
        quick_sort(v, temp, l, q);
        quick_sort(v, temp, q + 1, r);
    }

}

void set_array(vector<vector<std::complex<double>>>& array, vector<double> values, vector<int> bound_r, vector<int> bound_c){
    for(size_t i = bound_r[0];i < bound_r[1];i++){
        for(size_t j = bound_c[0]; j < bound_c[1];j++){
            array[i][j] = values[0];
        }
    }
}

std::vector<std::vector<std::complex<double>>> read_file(bool pen_writing, string path){
    //Test_Image.mat
    std::vector<std::vector<std::complex<double>>> v;
    std::string line, elem;
    std::vector<std::complex<double>> temp;
    std::complex<double> ch;
    double real, im;
    char sign;
    std::ifstream in_; // окрываем файл для чтения
    in_.open(path, std::ios::in);
    if (in_.is_open()) {
        while (getline(in_, line)) {
            std::stringstream ss, micro_ss;
            ss << line;
            while (ss >> elem) {

                if (pen_writing) {
                    micro_ss << elem;
                    micro_ss >> real;
                    micro_ss.get(sign);
                    micro_ss >> im;
                    if (sign == '-') {
                        im *= -1;
                    }
                    micro_ss.str(std::string());
                    ch = complex<double>(real, im);
                } else {
                    micro_ss << elem;
                    micro_ss >> ch;
                }
                temp.push_back(ch);
            }
            v.push_back(temp);
            temp.clear();
        }
    }
    in_.close();

    return v;
}

void Write_in_file(vector<vector<complex<double>>>& v, string file_name){
    std::ofstream out;          // поток для записи
    out.open(file_name + ".txt");
    size_t azimuth_size = v.size();
    size_t range_size = v[0].size();

    for(size_t i = 0; i < azimuth_size;i++){
        for(size_t j = 0; j < range_size;j++){
            if(j != range_size - 1){
                out << abs(v[i][j]) << ' ';
            }else {
                out << abs(v[i][j]);
            }
        }
        out << '\n';
    }
    out.close();
}






