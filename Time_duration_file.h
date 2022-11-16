
#ifndef OMEGA_K_TIME_DURATION_FILE_H
#define OMEGA_K_TIME_DURATION_FILE_H


#include <ctime>
#include <iostream>
#include <string>

class TimeDuration{
    clock_t end;
    clock_t start;
    std::string prog_name;
    public:
        TimeDuration(std::string n):prog_name(n){
            start = clock();
        }
        ~TimeDuration(){
            end = clock();
            double seconds = (double)(end - start) / CLOCKS_PER_SEC;
            std::cerr << "Work time of " + prog_name + ": " << seconds << " s" << std::endl;
        }
};

#endif //OMEGA_K_TIME_DURATION_FILE_H
