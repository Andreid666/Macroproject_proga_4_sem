#include <iostream>
#include <thread>
#include <chrono>
#include <string>
#include <ctime>
#include <vector>

class Calc{
public:
    std::vector<long long int> v;
Calc(){}

void foo(int n){
    //std::cout << n << " start\n";
    long long x = 0;
    for (long long int i=0; i < 1000000000/*9000000000*/; i++){
        x += 1 + v[n];
    }
    v[n] = x;
    //std::cout << n << " finish\n";
}

void proc(){
    time_t* speeds = new time_t[12];
    for (int i=0; i<12; i++){
        speeds[i] = hhh(i+1);
    }
    for (int i=0; i<12; i++){
        std::cout << speeds[i] << ' ';
    }
    std::cout << '\n';
}


time_t hhh(int cores){
    v = std::vector<long long int>(cores);
    time_t time = std::time(0);
    std::thread* helper = new std::thread[cores];
    for (int j=0; j<cores; j++){
        helper[j] = std::thread(&Calc::foo, this, j);
    }
    std::cout << "Тык\n";
    for (int i=0; i<cores; i++){
        helper[i].join();
    }
    time = std::time(0) - time;
    std::cout << cores << " ядер: " << time << "s\n"; 
    delete [] helper;
    /*for (int i=0; i<cores; i++){
        std::cout << v[i] << ' ';
    }
    std::cout << '\n';*/
    return time;
}
};
int main(){
    //time_t time = std::time(0);
    //foo(0);
    //speeds[0] = (std::time(0) - time);
    Calc c = Calc();
    c.proc();
}