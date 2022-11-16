#pragma once

#include <iostream>
#include <cmath>
#include <array>
#include <string>
#include <sstream>
#include <vector>
#include <fftw3.h>
#include <complex>
#include <complex.h>
#include <algorithm>
#include <fstream>

using namespace std;



vector<double> fill_up(const double& start,const double& finish, const double& step);
vector<double> fill_up(const double& start,const double& finish, size_t N);
vector<double> fill_up(const int& start,const int& finish, int N);
vector<int> fill_up(const int& start, const int& finish);

vector<double> real(vector<std::complex<double>>& data);
vector<double> imag(vector<std::complex<double>>& data);

void shift(vector<vector<std::complex<double>>> &Raw_data);

size_t partition(std::vector<double>& v, std::vector<int>& temp, size_t l, size_t r);
void quick_sort(std::vector<double>& v, std::vector<int>& temp, size_t l, size_t r);

void set_array(vector<vector<std::complex<double>>>& array, vector<double> values, vector<int> bound_r, vector<int> bound_c);

void conjugate(vector<complex<double>>& sopr, const vector<complex<double>>& orig);

vector<complex<double>> operator*(const vector<complex<double>>& v1, const vector<complex<double>>& v2);



std::vector<std::vector<std::complex<double>>> read_file(bool pen_writing, string path);

void Write_in_file(vector<vector<complex<double>>>& v, string file_name);


template<typename T>
void iFFTshift(vector<T>& v){
    size_t middle;
    if(v.size()%2 == 0){
        middle = v.size()/2;
    }else{
        middle = v.size()/2 + 1;
    }
    rotate(v.begin(),v.begin() + middle,v.end());
}

template<typename T>
T operator+(const T& v1, const T& v2){
	if(v1.size() != v2.size()){
		throw -1;
	}
	size_t count = v2.size();
	T temp(count);
	for(size_t i = 0; i < count;i++){
		temp[i] = v1[i] + v2[i];
	}
	return temp;
}

template<typename T>
std::vector<T> operator-(std::vector<T>& v1, T a){
    std::vector<T> temp(v1.size());
    for(size_t i = 0; i < v1.size();i++){
        temp[i] = v1[i]-a;
    }
    return temp;
}

template<typename T>
std::vector<T> operator/(std::vector<T>& v1, T a){
    std::vector<T> temp(v1.size());
    for(size_t i = 0; i < v1.size();i++){
        temp[i] = v1[i]/a;
    }
    return temp;
}

template<typename T>
vector<T> operator-(vector<T>& v1, vector<T>& v2){
	if(v1.size() != v2.size()){
		throw -1;
	}
	size_t count = v1.size();
	T temp(count);
	for(size_t i = 0; i < count;i++){
		temp[i] = v1[i]-v2[i];
	}
	return temp;
}

template<typename T, int l>
array<T,l> operator-(array<T,l>& v1, array<T,l>& v2){
	size_t count = v1.size();
	array<T,l> temp;
	for(size_t i = 0; i < count;i++){
		temp[i] = v1[i]-v2[i];
	}
	return temp;
}



template<typename T, typename D>
T operator*(const T& v, D a){
	size_t size = v.size();
    T temp(size);
	for(size_t i = 0; i < size;i++){
		temp[i] = v[i] * a;
	}
	return temp;
}

template<typename T>
T operator*(const T& v1, const T& v2){
    T temp(v1.size());
    if(v1.size() != v2.size()){
        throw("Not equal size!");
    }
    for(size_t i = 0; i < v1.size();i++){
        temp[i] = v1[i]*v2[i];
    }
    return temp;
}



template<typename T>
ostream& operator<<(ostream& os, const vector<T>& v){
	os << "{";
	bool flag = false;
	for(auto a: v){
		if(!flag){
			os << a;
			flag = true;
		}else{
			os << ", " << a;
		}
	}
	os << "}";
	return os;
}

template<typename T, long long unsigned int l>
ostream& operator<<(ostream& os, const array<T,l>& v){
	os << "[";
	bool flag = false;
	for(auto a: v){
		if(!flag){
			os << a;
			flag = true;
		}else{
			os << ", " << a;
		}

	}
	os << "]";
	return os;
}


