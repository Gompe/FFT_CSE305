#pragma once

#include <iostream>
#include <chrono>
#include <cstring>
#include <functional>

class Timer
{
public:
    Timer(){
        startTime = std::chrono::high_resolution_clock::now();
    }

    Timer(std::string title) : title(title){
        startTime = std::chrono::high_resolution_clock::now();
    }

    ~Timer(){
        using namespace std::chrono_literals;
        std::chrono::duration<float> duration = std::chrono::high_resolution_clock::now() - startTime;

        std::cout << "Timer " << title << ": " << (int) (duration / 1.0us) << "us" << std::endl;
    }
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    std::string title;
};

void timeFunction(std::function<void()> f, std::string title = ""){
    Timer timer = Timer(title);
    f();
}