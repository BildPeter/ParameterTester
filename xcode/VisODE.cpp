//
//  VisODE.cpp
//  ParameterTester
//
//  Created by BildPeter Visuals on 31.01.12.
//  Copyright (c) 2012 BildPeter Visuals. All rights reserved.
//

#include <boost/numeric/odeint.hpp>
#include <algorithm>

#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"
#include "cinder/ImageIO.h"
#include "cinder/gl/Texture.h"
#include "cinder/Perlin.h"
#include "cinder/Channel.h"
#include "cinder/Vector.h"
#include "cinder/Utilities.h"
#include "cinder/params/Params.h"

using namespace boost::numeric::odeint;

using namespace ci;
using namespace ci::app;
using namespace std;

/* --------------------------------------------------------------------------------------------------------
 AIM:
 
 TODO:
 
 - Rescale the axis
 
 -------------------------------------------------------------------------------------------------------- */

void drawPlot( vector< double > &inputVec, Color drawColor );

typedef vector< double >    State;
float   beta    = 20.0;
float   lamda   = 1.0;

void sirDynamic( const State &x, State &dx, const double t ) {
    
        double N = x[ 0 ] + x[ 1 ] + x[ 2 ];
        dx[ 0 ] = ( (- 1) * x[ 0 ] * x[ 1 ] * beta ) / N ; // dS = - S * I * beta
        dx[ 1 ] = ( x[ 0 ] * x[ 1 ] * beta ) / N - lamda * x[ 1 ] ; // dI = 
        dx[ 2 ] = lamda * x[ 1 ];                            // dR = lambda * I
    }




class VisODE : public AppBasic {
public:
	void setup();
	void update();
	void draw();
    void toggleSim();
    void keyDown( KeyEvent event );
    void setInitialSIR( State &stateVec, double inS, double inI, double inR );
    
    
    State           x;
    State           xTest;
    State           dxdt;
    vector< double >    SOE;
    vector< double >    Diff;
    //        dynState = State( sirSys.totalStates() );
    
    int     steps       = 1;
    float   stepSize    = 0.01;
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
    setWindowSize( 1300, 700);
    gl::enableAlphaBlending();
    
    x       = vector< double >(3);
    dxdt    = vector< double >(3);
    xTest   = vector< double >(3);
    
    setInitialSIR( x, initS, initI, initR );
    setInitialSIR( xTest, initS, initI, initR);
    
    // -----------------------------------------------
    // ----- Parameter Interface -----
    mParams = params::InterfaceGl( "Parameters", Vec2i( 200, 200 ) );
    mParams.addParam( "Steps Size",         &stepSize,  "min=0.0001 step=0.0001 precision=4" );
    mParams.addParam( "Control Precision", &controlPrecision, "min=0.0 step=1");
    mParams.addParam( "Steps per output",   &steps,     "min=0.0 step=1" );
    mParams.addParam( "beta", &beta, "min=0.0");
    mParams.addParam( "lamda", &lamda, "min=0.0");
    mParams.addParam( "Initial S", &initS, "min=0.0");
    mParams.addParam( "Initial I", &initI, "min=0.0");
    mParams.addParam( "Initial R", &initR, "min=0.0");
}


void VisODE::update()
{
    if (runSim) {
        integrate( sirDynamic , x ,     0.0 , (double)(stepSize * steps ) , (double)stepSize );             // one Step
        integrate( sirDynamic , xTest , 0.0 , (double)(stepSize * steps ) , ( (double)stepSize/controlPrecision ) );    // Hundred Steps            
        
        SOE.push_back( x[ 1 ] );
        Diff.push_back( x[ 1 ] - xTest[ 1 ] );
        cout << "MaxVal Error: " << *max_element( Diff.begin(), Diff.end() ) << endl;
    }
    
}

void VisODE::draw()
{
	gl::clear( backgroundColor ); 
    
    // ----- Plot: Size of Epidemic -----
    drawPlot( SOE,  Color( 0, 1, 0 ) );
    drawPlot( Diff, Color( 1, 0, 0 ) );
    
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
            setInitialSIR( x, initS, initI, initR );
            setInitialSIR( xTest, initS, initI, initR );
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


CINDER_APP_BASIC( VisODE, RendererGl )


void drawPlot( vector< double > &inputVec, Color drawColor ){
    int     timeLine = 0;
    float headLength = 6.0f;
    float headRadius = 3.0f;
    Vec2f   coordinate;
    
    glLineWidth( 1 );
    gl::color( Color( 0,0,0 ) );
    
    gl::pushMatrices();
    gl::translate( Vec2f( 50, 650 ) );
    gl::rotate( Vec3f( 180, 0, 0) );
    
    gl::drawVector( Vec3f( 0, 0, 0 ) , Vec3f( 1200, 0, 0 ),  headLength , headRadius );     // X-Axis
    gl::drawVector( Vec3f( 0, 0, 0 ) , Vec3f( 0, 300, 0 ), headLength , headRadius );       // Y-Axis
    
    glLineWidth( 2 );
    gl::color( drawColor );
    gl::scale( Vec3f( 1, 0.1, 1 ) );
    glBegin(GL_LINE_STRIP);    
    for (vector< double >::iterator iter = inputVec.begin(); iter != inputVec.end(); iter++){
        coordinate = Vec2f( timeLine, *iter );
        //        cout << "Output: " << *iter << endl;
        gl::vertex( coordinate );
        timeLine++;
    }
    glEnd();
    gl::popMatrices();
}






