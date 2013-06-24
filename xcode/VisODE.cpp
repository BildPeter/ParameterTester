//
//  VisODE.cpp
//  ParameterTester
//
//  Created by Peter A. Kolski on 31.01.12.
//  Copyright (c) 2012 Peter A. Kolski. All rights reserved.
//

// --- Boost (has to be before Cinder)
#include <boost/numeric/odeint.hpp>
#include <boost/lexical_cast.hpp>

// --- C++
#include <algorithm>

// --- Cinder
#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"
#include "cinder/ImageIO.h"
#include "cinder/gl/Texture.h"
#include "cinder/Perlin.h"
#include "cinder/Channel.h"
#include "cinder/Vector.h"
#include "cinder/Utilities.h"
#include "cinder/params/Params.h"

// --- Peter Extension
#include "stepCounter.h"
#include "peterVis.h"

using namespace boost::numeric::odeint;

using namespace ci;
using namespace ci::app;
using namespace std;

/* --------------------------------------------------------------------------------------------------------
 AIM:
 
 TODO:
 
 - Rescale the axis
 
 -------------------------------------------------------------------------------------------------------- */

typedef vector< double >    State;
float   beta    = 0.42;
float   lamda   = 0.097;


void sirDynamic( const State &x, State &dx, const double t ) {
    
        double N = x[ 0 ] + x[ 1 ] + x[ 2 ];
        dx[ 0 ] = ( (- 1) * x[ 0 ] * x[ 1 ] * beta ) / N ; // dS = - S * I * beta
        dx[ 1 ] = ( x[ 0 ] * x[ 1 ] * beta ) / N - lamda * x[ 1 ] ; // dI = 
        dx[ 2 ] = lamda * x[ 1 ];                            // dR = lambda * I
    }

void sisDynamic( const State &x, State &dx, const double t ) {
    
    double N = x[ 0 ] + x[ 1 ];
    
    dx[ 0 ] = ( lamda * x[ 1 ] - x[ 0 ] * x[ 1 ] * beta ) / N ; // dS = I * lamda    - S * I * beta
    dx[ 1 ] = ( x[ 0 ] * x[ 1 ] * beta - lamda * x[ 1 ] ) / N;  // dI = S * I * beta - lambda * I
}


class VisODE : public AppBasic {
public:
	void setup();
	void update();
	void draw();
    void toggleSim();
    void keyDown( KeyEvent event );
    void setInitialSIR( State &stateVec, double inS, double inI, double inR );
    void setInitialSIS( State &stateVec, double inS, double inI );
    
    
    State           x;
    State           xTest;
    State           dxdt;
    vector< double >    SOE;
    vector< double >    Diff;
    //        dynState = State( sirSys.totalStates() );
    
    StepCounter sc;
    
    float   controlPrecision = 100;
    bool    runSim      = false;
    float  initS       = 300.0;
    float  initI       = 1.0;
    float  initR       = 0.0;
    float   alpha      = 1.0;
    double   t       = 0.0;
    
    params::InterfaceGl		mParams;
    Color backgroundColor = Color( 1, 1, 1 );
    
    runge_kutta4< State >   rk;
    
};

void VisODE::setup()
{
    setWindowSize( 500, 700);
    gl::enableAlphaBlending();
    
//    x       = vector< double >(3);
//    dxdt    = vector< double >(3);
//    xTest   = vector< double >(3);
//    
//    setInitialSIR( x, initS, initI, initR );
//    setInitialSIR( xTest, initS, initI, initR);
    
    x       = vector< double >(2);
    dxdt    = vector< double >(2);
    xTest   = vector< double >(2);
    
    setInitialSIS( x, initS, initI );
    setInitialSIS( xTest, initS, initI);

    
    // -----------------------------------------------
    // ----- Parameter Interface -----
    mParams = params::InterfaceGl( "Parameters", Vec2i( 200, 200 ) );
    mParams.addParam( "Steps Size",         &sc.stepSize,  "min=0.0001 step=0.0001 precision=4" );
    mParams.addParam( "Control Precision", &controlPrecision, "min=0.0 step=1");
    mParams.addParam( "Steps per Run",   &sc.stepsPerRun,     "min=0.0 step=1" );
    mParams.addParam( "beta", &beta, "min=0.0 step=0.01");
    mParams.addParam( "alpha", &lamda, "min=0.0 step=0.001");
    mParams.addParam( "Initial S", &initS, "min=0.0");
    mParams.addParam( "Initial I", &initI, "min=0.0");
    mParams.addParam( "Initial R", &initR, "min=0.0");
}


void VisODE::update()
{
    if (runSim) {
//        integrate( sirDynamic , x ,     0.0 , (double)(sc.stepSize * sc.steps ) , (double)sc.stepSize );             // one Step
//        integrate( sirDynamic , xTest , 0.0 , (double)(sc.stepSize * sc.steps ) , ( (double)sc.stepSize/controlPrecision ) );    // Hundred Steps

        integrate( sisDynamic , x ,     0.0 , (double)(sc.stepSize * sc.stepsPerRun ) , (double)sc.stepSize );             // one Step
        integrate( sisDynamic , xTest , 0.0 , (double)(sc.stepSize * sc.stepsPerRun ) , ( (double)sc.stepSize/controlPrecision ) );    // Hundred Steps
        
        sc.runDone();   // calculetes the values for time&run counting
        
        SOE.push_back( x[ 1 ] );
        Diff.push_back( x[ 1 ] - xTest[ 1 ] );
        cout << "MaxVal Error: " << *max_element( Diff.begin(), Diff.end() ) << endl;
    }
    
}

void VisODE::draw()
{
	gl::clear( backgroundColor ); 
    
    string stepsTxt = "Passed Steps: ";
    string timeTxt = "Passed Virtual Time: ";
    stepsTxt += boost::lexical_cast<string>( sc.passedSteps );
    timeTxt += boost::lexical_cast<string>( sc.passedTime );
    
    // ----- Plot: Size of Epidemic -----
    drawPlot( SOE, Vec2i( 50, 400 ), 1.0, "SOE");
    drawPlot( Diff, Vec2i( 50, 600 ), 1.0, "Numerical precision difference");
    
    gl::drawString( "SIS Model Parameter Tester", Vec2i( getWindowWidth() / 2 + 50, 50 ), Color(0, 0, 0) );
    gl::drawString( stepsTxt, Vec2i( getWindowWidth() / 2 + 50, 95 ), Color(0, 0, 0) );
    gl::drawString( timeTxt, Vec2i( getWindowWidth() / 2 + 50, 110 ), Color(0, 0, 0) );
    
    params::InterfaceGl::draw();
}


void VisODE::keyDown( KeyEvent event )
{
    switch ( event.getChar() ) {
        case ' ':
            toggleSim();
            break;
        case 'r':
            SOE.clear();
            Diff.clear();
            setInitialSIS( x, initS, initI );
            setInitialSIS( xTest, initS, initI );
            sc.resetCounter();
            break;
        default:
            break;
    }
}


void VisODE::toggleSim()
{
    if (runSim){
        runSim = false;
    }
    else
        runSim = true;
}

void VisODE::setInitialSIR( State &stateVec, double inS, double inI, double inR )
{
    stateVec[ 0 ] = inS;
    stateVec[ 1 ] = inI;
    stateVec[ 2 ] = inR;
}

void VisODE::setInitialSIS( State &stateVec, double inS, double inI )
{
    stateVec[ 0 ] = inS;
    stateVec[ 1 ] = inI;
}



CINDER_APP_BASIC( VisODE, RendererGl )