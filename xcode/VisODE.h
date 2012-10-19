//
//  VisODE.h
//  ParameterTester
//
//  Created by Peter on 19.10.12.
//
//

#ifndef ParameterTester_VisODE_h
#define ParameterTester_VisODE_h

class StepCounter{
public:
    StepCounter();
    
    void runDone();
    void resetCounter();
    
    double  stepSize;
    int     stepsPerRun;
    int     passedSteps;
    double  passedTime;
    
    
};

StepCounter::StepCounter() {

    stepSize    = 0.01;
    stepsPerRun = 1;
    passedSteps = 0;
    passedTime  = 0;
    

}

void StepCounter::runDone(){

    passedSteps += stepsPerRun;
    passedTime  = passedSteps * stepSize;

}

void StepCounter::resetCounter(){

    passedSteps = 0;
    passedTime  = 0;

}

#endif
