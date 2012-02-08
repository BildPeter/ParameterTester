#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"
#include "cinder/ImageIO.h"
#include "cinder/gl/Texture.h"
#include "cinder/Perlin.h"
#include "cinder/Channel.h"
#include "cinder/Vector.h"
#include "cinder/Utilities.h"
#include "cinder/params/Params.h"

#include <netevo/netevo.h>

//#include <numeric/odeint/stepper/runge_kutta4.hpp>
//#include "peter/PeterDynamics.h"

using namespace netevo;
//using namespace boost::numeric;

using namespace ci;
using namespace ci::app;
using namespace std;

/* --------------------------------------------------------------------------------------------------------
 AIM:
 
 TODO:
 - Stepsize should be variable 
        -> Do_Step implementation
 - Include ODEint with Do_Step() instead of Netevo
 - Rescale the axis
 - Difference plot
        -> New plot
        -> additional precision
 
 -------------------------------------------------------------------------------------------------------- */
 

void initDegreeDistributed( System &sys, State &initial, int states, double degreeFactor );

//class SimulateOdeStep   : public SimulateOdeFixed{
//    public  :
//    SimulateOdeStep( double StepSize );   // Konstruktor
//    void simulateStep(System &sys, double time, State &inital );
//    
//    private :
//    odeint::runge_kutta4< State > mRk4Stepper;
//    double mStepSize;
//};
//
//SimulateOdeStep::SimulateOdeStep( double StepSize ) : SimulateOdeFixed( RK_4, StepSize ){
//    mStepSize = StepSize;
//}
//
//void SimulateOdeStep::simulateStep(System &sys, double time, State &inital ){
//    
//    if (!sys.validStateIDs()) { sys.refreshStateIDs(); }
//    
//    mRk4Stepper.do_step( Simulator(&sys), inital,  time, mStepSize ); 
//    
//}



class sirDynamics : public NodeDynamic {
    
public:
    string  getName()   { return "SIR"; }
    int     getStates() { return 3; }
    void    setDefaultParams( Node n, System &sys ){
        sys.nodeData( n ).dynamicParams.push_back( 20.0 );  // beta
        sys.nodeData( n ).dynamicParams.push_back( 1.0 ); // lamda
    }
    void fn( Node n, System &sys, const State &x, State &dx, const double t ) {
        
        int nID = sys.stateID( n );
        double N = x[ nID ] + x[ nID + 1 ] + x[ nID + 2 ];
        vector< double > &nParams = sys.nodeData( n ).dynamicParams;
        
        dx[  nID      ] = ( (- 1) * x[ nID ] * x[ (nID + 1)] * nParams[ 0 ] ) / N ; // dS = - S * I * beta
        dx[ (nID + 1) ] = ( x[ nID ] * x[ (nID + 1)] * nParams[ 0 ] ) / N - nParams[ 1 ] * x[ (nID + 1) ] ; // dI = 
        dx[ (nID + 2) ] = nParams[ 1 ] * x[ (nID + 1) ];                            // dR = lambda * I
    }
};

void drawPlot( vector< double > &inputVec );

class MarriageApp : public AppBasic {
  public:
	void setup();
	void update();
	void draw();
    void toggleSim();
    void keyDown( KeyEvent event );
    void setInitialSIR( State &stateVec, double inS, double inI, double inR );
    
    System              sirSys;
    sirDynamics         sirDyn;
    ChangeLog           nullLog;
    SimObserver         nullObserver;
    State           dynState;
    vector<double>  SOE;
//        dynState = State( sirSys.totalStates() );
    
    int     steps       = 1;
    float  stepSize    = 0.0001;
    bool    runSim      = false;
    float  initS       = 300.0;
    float  initI       = 1.0;
    float  initR       = 0.0;
    
    SimulateOdeFixed sirSimulator = SimulateOdeFixed(RK_4, stepSize );
    
    params::InterfaceGl		mParams;
    Color backgroundColor = Color( 1, 1, 1 );
    
};

void MarriageApp::setup()
{
    setWindowSize( 1300, 700);
    gl::enableAlphaBlending();

    // ----- Dynamik zuweisen und Graph kopieren (Reihenfolge beachten) -----
    sirSys.addNodeDynamic( &sirDyn );   
    sirSys.addNode( "SIR" );
    
    // ----- Anfangsbedingungen festlegen -----
    dynState = State( sirSys.totalStates() );
    setInitialSIR( dynState, initS, initI, initR );
    
    // -----------------------------------------------
    // ----- Parameter Interface -----
    mParams = params::InterfaceGl( "Parameters", Vec2i( 200, 200 ) );
    mParams.addParam( "Steps Size",         &stepSize,  "min=0.0001 step=0.0001 precision=4" );
    mParams.addParam( "Steps per output",   &steps,     "min=0.0 step=1" );
    mParams.addParam( "Initial S", &initS, "min=0.0");
    mParams.addParam( "Initial I", &initI, "min=0.0");
    mParams.addParam( "Initial R", &initR, "min=0.0");
    
}


void MarriageApp::update()
{
    if (runSim) {
//    SimObserverToStream coutObs( std::cout );
//    sirSimulator.simulate( sirSys, 0.01,  dynState,  coutObs,  nullLog );
    sirSimulator.simulate(    sirSys, steps * stepSize, dynState, nullObserver, nullLog );
    SOE.push_back( dynState[ 1 ] );
//    cout << "Nr: " << getElapsedFrames() << "\tSizeOfEpidemic: " << SOE.back() << endl;
    }
}




void MarriageApp::draw()
{
	gl::clear( backgroundColor ); 
    
    // ----- Plot: Size of Epidemic -----
    drawPlot( SOE );
    
    params::InterfaceGl::draw();
    }


void MarriageApp::keyDown( KeyEvent event )
{
    switch ( event.getChar() ) {
        case ' ':
            toggleSim();
            break;
        case 'r':
            SOE.clear();
            setInitialSIR( dynState, initS, initI, initR );
            break;
        default:
            break;
    }
}

void MarriageApp::toggleSim()
{
    if (runSim){
        runSim = false;
    }
    else
        runSim = true;
}

void MarriageApp::setInitialSIR( State &stateVec, double inS, double inI, double inR )
{
    stateVec[ 0 ] = inS;
    stateVec[ 1 ] = inI;
    stateVec[ 2 ] = inR;
}


CINDER_APP_BASIC( MarriageApp, RendererGl )

void initDegreeDistributed( System &sys, State &initial, int states, double degreeFactor ){
    int nodeNr = 0;     // Index des Knoten, um keinen Iterator zu benutzen
    
    for ( System::NodeIt n (sys); n != INVALID; ++n ){
        unsigned int degree = 0;
        for ( System::InArcIt a( sys, n ); a != INVALID; ++a) {
            degree++;
        }
        //        initial[ (nodeNr * states ) ] = degree * degreeFactor + 10;   //Anzahl des Grades 
        initial[ (nodeNr * states ) ] = 30;   //Anzahl des Grades 
        nodeNr++;
    }
    //    for_each( initial.begin(), initial.end(),printCout );
}

void drawPlot( vector< double > &inputVec ){
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
    gl::color( Color( 1, 0, 0 ) );
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






