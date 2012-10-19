//
//  peterVis.h
//  NE_Vis_Cinder
//
//  Created by BildPeter Visuals on 18.04.12.
//  Copyright (c) 2012 BildPeter Visuals. All rights reserved.
//

#include <boost/lexical_cast.hpp>
#include <algorithm>
#include <cmath>
#include "cinder/Vector.h"

using namespace ci;
using namespace ci::app;
using namespace std;


#ifndef NE_Vis_Cinder_peterVis_h
#define NE_Vis_Cinder_peterVis_h


void drawPlot( vector< double > &inputVec , Vec2f pos, int scale, string title){
    int   mHeight = 150;
    
    int   timeLine = 0;
    float headLength = 6.0f;
    float headRadius = 3.0f;
    Vec2f   coordinate;
    string  axisY;
    vector< double >::iterator  maxIter;
    double mScale = 1;
    
    axisY = "Y: ";
        
    if ( inputVec.size() != 0 ){
        maxIter = max_element( inputVec.begin(), inputVec.end() );
        double maxVal = (int)(*maxIter * 100);
        maxVal = maxVal /100;
//        double maxVal = fmod( *maxIter, 2);
        axisY += boost::lexical_cast< string >( maxVal );
        // --- scale the y-Axis
        if (  maxVal > 0  ){
            mScale = 10 * mHeight / maxVal;
        }
    }
    
    gl::pushMatrices();
    gl::translate( pos );
    gl::rotate( Vec3f( 180, 0, 0) );
    
    // --- Background
    gl::color( ColorA( 0.8, 0.8, 0.8, 0.7 ) );
    gl::drawSolidRect( Rectf( 380, 0, 0, mHeight) );

    // --- Axis
    glLineWidth( 1 );
    gl::color( Color (0,0,0) );
    
    gl::drawVector( Vec3f( 0, 0, 0 ) , Vec3f( 380, 0, 0 ),  headLength , headRadius );     // x-Axis
    gl::drawVector( Vec3f( 0, 0, 0 ) , Vec3f( 0, mHeight, 0 ), headLength , headRadius );     // Y-Axis
    
    // --- 
    glLineWidth( 2 );
    gl::color( Color( 1, 0, 0 ) );
    gl::scale( Vec2f( 0.7, 0.1 ) );
    glBegin(GL_LINE_STRIP);    
    for (vector< double >::iterator iter = inputVec.begin(); iter != inputVec.end(); iter++){
        coordinate = Vec2f( timeLine, (int)( (*iter) * mScale) );
        gl::vertex( coordinate );
        timeLine++;
    }
    glEnd();
    gl::popMatrices();

    // --- Title
    gl::drawString( title, (pos - Vec2i( - 200, 163) ), Color(0,0,0) );
    gl::drawString( axisY, pos - Vec2i( 0, 163 ), Color( 0, 0, 0) );

}

#endif
