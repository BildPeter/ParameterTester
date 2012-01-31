#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"

using namespace ci;
using namespace ci::app;
using namespace std;

class ParameterTesterApp : public AppBasic {
  public:
	void setup();
	void mouseDown( MouseEvent event );	
	void update();
	void draw();
};

void ParameterTesterApp::setup()
{
}

void ParameterTesterApp::mouseDown( MouseEvent event )
{
}

void ParameterTesterApp::update()
{
}

void ParameterTesterApp::draw()
{
	// clear out the window with black
	gl::clear( Color( 0, 0, 0 ) ); 
}


CINDER_APP_BASIC( ParameterTesterApp, RendererGl )
