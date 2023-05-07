#include <thread>
#include <iostream>
#include <cmath>
#include <stack>
#include <limits>
#include <semaphore.h>
#include <iomanip>


#define LOCAL_TO_GLOBAL_TRANSFER_COUNT 32
#define MAX_GLOBAL_STACK_SIZE 16

int NUM_THREADS = 4;


constexpr double start_point = 0.0;
constexpr double end_point   = 0.99999999;


struct IntegrateSegmentInfo {
    double start;
    double end;

    double f_start;
    double f_end;

    double area;
};


struct ThreadInfo {
    sem_t* p_semStackAccess;
    sem_t* p_semSumAccess;
    sem_t* p_semSegmentPresent;

    std::stack<IntegrateSegmentInfo>* p_globalStack;

    double* p_totalSum;

    double (*func)(const double); 
    double epsilon;

    int* p_activeThreads;
};


double f(const double x) { return sin(1/(1-x)); }



double IntegrateSingle(const double start, const double end, double (*func)(const double), const double epsilon) {
    double totalSum = 0;

    double A = start;
    double B = end;
    
    double fA = func(A);
    double fB = func(B);

    std::stack<IntegrateSegmentInfo> localStack;

    double sAB = (fA + fB) * (end - start) / 2;

    while(true) {
        const double C = (A + B) / 2;

        const double fC = func(C);

        const double sAC = (fA + fC) * (C - A) / 2;
        const double sCB = (fC + fB) * (B - C) / 2;

        const double sACB = sAC + sCB;


        const double diff = std::abs(sAB - sACB);
        if(diff >= epsilon * std::abs(sACB) && diff > std::numeric_limits<double>::epsilon()) {
            localStack.push(IntegrateSegmentInfo{A, C, fA, fC, sAC});

            A = C;
            fA = fC;
            sAB = sCB;
        }
        else {
            totalSum += sACB;

            if(localStack.empty())
                break;

            auto newSegment = localStack.top();

            A = newSegment.start;
            B = newSegment.end;
            fA = newSegment.f_start;
            fB = newSegment.f_end;
            sAB = newSegment.area;

            localStack.pop();
        }
    }

    return totalSum;
}



void threadFunc(ThreadInfo* info) {
    auto& globalStack = *(info->p_globalStack);
    auto& numActiveThreads = *(info->p_activeThreads);

    auto func = info->func;
    auto epsilon = info->epsilon;

    double localSum = 0.0;

    std::stack<IntegrateSegmentInfo> localStack;

    double sAB;

    double A;
    double B;
    double fA;
    double fB;

    while(true) {
        sem_wait(info->p_semSegmentPresent);

        sem_wait(info->p_semStackAccess);

        auto newSegment = globalStack.top();

        A = newSegment.start;
        B = newSegment.end;
        fA = newSegment.f_start;
        fB = newSegment.f_end;
        sAB = newSegment.area;

        globalStack.pop();

        if(!globalStack.empty())
            sem_post(info->p_semSegmentPresent);

        if(A <= B)
            ++numActiveThreads;

        sem_post(info->p_semStackAccess);

        if(A > B)
            break;


        while(true) {

            const double C = (A + B) / 2;

            const double fC = func(C);

            const double sAC = (fA + fC) * (C - A) / 2;
            const double sCB = (fC + fB) * (B - C) / 2;

            const double sACB = sAC + sCB;

            const double diff = std::abs(sAB - sACB);
            if(diff >= epsilon * std::abs(sACB) && diff > std::numeric_limits<double>::epsilon()) {
                localStack.push(IntegrateSegmentInfo{A, C, fA, fC, sAC});

                A = C;
                fA = fC;
                sAB = sCB;
            }
            else {
                localSum += sACB;

                if(localStack.empty())
                    break;

                auto newSegment = localStack.top();

                A = newSegment.start;
                B = newSegment.end;
                fA = newSegment.f_start;
                fB = newSegment.f_end;
                sAB = newSegment.area;

                localStack.pop();
            }

            if(localStack.size() > LOCAL_TO_GLOBAL_TRANSFER_COUNT) {
                sem_wait(info->p_semStackAccess);
                if(globalStack.empty()) {
                    while(localStack.size() > 1 && globalStack.size() < MAX_GLOBAL_STACK_SIZE) {
                        globalStack.push(localStack.top());
                        localStack.pop();
                    }
                    sem_post(info->p_semSegmentPresent);
                }
                sem_post(info->p_semStackAccess);
            }
        }


        sem_wait(info->p_semStackAccess);
        
        --numActiveThreads;

        if(!numActiveThreads && globalStack.empty()) {
            for(int i = 0; i < NUM_THREADS; ++i)
                globalStack.push(IntegrateSegmentInfo{2.0, 1.0, 0.0, 0.0, 0.0});

            sem_post(info->p_semSegmentPresent);
        }
        sem_post(info->p_semStackAccess);
    }

    sem_wait(info->p_semSumAccess);
    *(info->p_totalSum) += localSum;
    sem_post(info->p_semSumAccess);
    
    return;
}



double IntegrateThreads(const double start, const double end, double (*func)(const double), const double epsilon) {
    double totalSum = 0;

    std::stack<IntegrateSegmentInfo> globalStack;

    ThreadInfo infos[NUM_THREADS];
    std::thread threads[NUM_THREADS];

    double fA = func(start);
    double fB = func(end);

    globalStack.push(IntegrateSegmentInfo{start, end, fA, fB, (fA + fB) * (end - start) / 2});

    sem_t semStackAccess;
    sem_t semSumAccess;
    sem_t semSegmentPresent;

    sem_init(&semStackAccess, 0, 1);
    sem_init(&semSumAccess, 0, 1);
    sem_init(&semSegmentPresent, 0, 1);

    int activeThreads = 0;

    for(int i = 0; i < NUM_THREADS; ++i) {
        infos[i].p_semStackAccess = &semStackAccess;
        infos[i].p_semSumAccess = &semSumAccess;
        infos[i].p_semSegmentPresent = &semSegmentPresent;
        infos[i].p_globalStack = &globalStack;
        infos[i].p_totalSum = &totalSum;
        infos[i].func = func;
        infos[i].epsilon = epsilon;
        infos[i].p_activeThreads = &activeThreads;

        threads[i] = std::thread(threadFunc, &infos[i]);
    }


    for(int i = 0; i < NUM_THREADS; ++i)
        threads[i].join();


    sem_destroy(&semStackAccess);
    sem_destroy(&semSumAccess);
    sem_destroy(&semSegmentPresent);

    return totalSum;
}



int main(int argc, char* argv[])
{
    double integralValue = 0.0;

    if(argc < 3) {
        printf("Usage: %s <num_threads> <accuracy>\n", argv[0]);
        return -1;
    }

    NUM_THREADS = atoi(argv[1]);
    double accuracy = atof(argv[2]);

    if(NUM_THREADS == 1)
        integralValue = IntegrateSingle(start_point, end_point, f, accuracy);
    else
        integralValue = IntegrateThreads(start_point, end_point, f, accuracy);
    
    std::cout << "Integral value = " << std::setprecision(16) << integralValue << std::endl;

    return 0;
}
