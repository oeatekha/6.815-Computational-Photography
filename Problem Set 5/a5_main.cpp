/* -----------------------------------------------------------------
 * File:    a5_main.cpp
 * Author:  Michael Gharbi <gharbi@mit.edu>
 * Created: 2015-09-30
 * -----------------------------------------------------------------
 *
 *
 *
 * ---------------------------------------------------------------*/

#include "Image.h"
#include "basicImageManipulation.h"
#include "morphing.h"
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

void testNearestNeighbor() {
  // Test nearest neighbor
  const Image im("./Input/BostonRainbow-crop-400.png");
  float fact = 3.5;
  Image scaledNN = scaleNN(im, fact);
  scaledNN.write("./Output/testNN.png");
}

void testBilinearInterpolation() {
  // Test bilinear interpolation
  cout << endl << "linear interpolation" << endl;
  cout << "--------------------" << endl;
  Image test(2, 2, 1);
  test(0, 0, 0) = 0.0f;
  test(0, 1, 0) = 0.25f;
  test(1, 0, 0) = 0.5f;
  test(1, 1, 0) = 1.0f;

  float x = 0.25f, y = 0.0;
  cout << "interpolate at (" << x << ", " << y
       << "): " << interpolateLin(test, x, y, 0, false) << ", should be 0.125"
       << endl;

  x = 0.0;
  y = 0.25f;
  cout << "interpolate at (" << x << ", " << y
       << "): " << interpolateLin(test, x, y, 0, false) << ", should be 0.0625"
       << endl;

  x = 0.5;
  y = 0.5f;
  cout << "interpolate at (" << x << ", " << y
       << "): " << interpolateLin(test, x, y, 0, false) << ", should be 0.4375"
       << endl;
}

void testBilinearRescaling() {
  // Test bilinear
  const Image im("./Input/BostonRainbow-crop-400.png");
  float fact = 3.5;
  Image scaled = scaleLin(im, fact);
  scaled.write("./Output/testLin.png");
}

void testRotation() {
  const Image im("./Input/BostonRainbow-crop-400.png");

  float theta = 3.14159265358979323846 / 8; // M_PI

  Image rot = rotate(im, theta);
  rot.write("./Output/testRotate.png");
}

void testVectorOperations() {
  // Test vector lib
  Vec2f a(10.0f, 3.0f);
  Vec2f b(-4.1f, 2.7f);
  float f = 2.0f;
  cout << endl << "vector operations" << endl;
  cout << "-----------------" << endl;

  cout << "a(" << a.x << ", " << a.y << ")" << endl;
  cout << "b(" << b.x << ", " << b.y << ")" << endl;
  cout << "f=" << f << endl;

  // a+b
  Vec2f c = a + b;
  cout << "a+b: (" << c.x << ", " << c.y << ")" << endl;

  // a-b
  c = a - b;
  cout << "a-b: (" << c.x << ", " << c.y << ")" << endl;

  // a*f
  c = a * f;
  cout << "a*f: (" << c.x << ", " << c.y << ")" << endl;

  // a/f
  c = a / f;
  cout << "a/f: (" << c.x << ", " << c.y << ")" << endl;

  // <a,b>
  float s = dot(a, b);
  cout << "<a,b>=" << s << endl;

  // ||a||
  float l = length(a);
  cout << "||a||=" << l << endl;

  Vec2f p_a = perpendicular(a);
  cout << "<a, perpendicular(a)>=" << dot(a, p_a) << endl;
}

void testSegment() {
  // Test segment
  // are P,Q,e1 e1 correctly implemented?
  // are the u,v coordinates meaningful?
  // What should be u and v for P,Q ?
  // Come up with a few test cases !
    Vec2f P = { 0,0 };
    Vec2f Q = { 1,2 };
    Vec2f X = { 1,1 }; 
    Vec2f X2 = { -5,-5 };
    Segment PQ_VEC(P, Q);
    //Segment segBefore(Vec2f(81, 105), Vec2f(122, 112));

    float b = PQ_VEC.distance(X);
    float a = PQ_VEC.distance(X2);
    cout << "vec len is " << PQ_VEC.getLength() << endl;

    cout << "B is 1,1 " << b << endl;
    cout << "A is 5,5 " << a << endl;




}

void testWarpBy1() {
  // Test warpBy1 ----------------------------------
  Image bear("./Input/bear.png");

  // define before and after segments
  Segment segBefore(Vec2f(0, 0), Vec2f(10, 0)); //= Segment((0,0), (10,0))
  Segment segAfter(Vec2f(10, 10), Vec2f(30, 15)); //Segment((10, 10), (30, 15))


  // warp
  Image warp1 = warpBy1(bear, segBefore, segAfter);
  warp1.write("bear.png");
  // ------------------------------------------------
}

void testWarp() {
  // Make your own test !

    vector<Segment> segsBefore;
    segsBefore.push_back(Segment(Vec2f(56, 111), Vec2f(69, 90)));
    segsBefore.push_back(Segment(Vec2f(97, 98), Vec2f(119, 104)));


    vector<Segment> segsAfter;
    segsAfter.push_back(Segment(Vec2f(61, 103), Vec2f(78, 76)));
    segsAfter.push_back(Segment(Vec2f(93, 93), Vec2f(114, 101)));


    Image fredo("./Input/fredo2.png");
    Image werewolf("./Input/werewolf.png");
    Image output = warp(fredo, segsBefore, segsAfter);
    output.write("Checkwarp.png");
}

void testMorph() {
  // Test morph
  Image fredo("./Input/fredo2.png");
  Image werewolf("./Input/werewolf.png");

  // ... use the UI to specify segments
  vector<Segment> segsBefore;
  segsBefore.push_back(Segment(Vec2f(169, 255), Vec2f(230, 253)));
  segsBefore.push_back(Segment(Vec2f(344, 261), Vec2f(292, 255)));
  segsBefore.push_back(Segment(Vec2f(124, 371), Vec2f(247, 537)));
  segsBefore.push_back(Segment(Vec2f(293, 528), Vec2f(358, 377)));
  segsBefore.push_back(Segment(Vec2f(268, 241), Vec2f(274, 283)));
  segsBefore.push_back(Segment(Vec2f(274, 344), Vec2f(272, 469)));
  segsBefore.push_back(Segment(Vec2f(221, 407), Vec2f(312, 404)));
  segsBefore.push_back(Segment(Vec2f(120, 222), Vec2f(121, 351)));


  vector<Segment> segsAfter;
  segsAfter.push_back(Segment(Vec2f(169, 244), Vec2f(207, 248)));
  segsAfter.push_back(Segment(Vec2f(315, 240), Vec2f(271, 243)));
  segsAfter.push_back(Segment(Vec2f(150, 338), Vec2f(209, 395)));
  segsAfter.push_back(Segment(Vec2f(262, 392), Vec2f(334, 344)));
  segsAfter.push_back(Segment(Vec2f(242, 249), Vec2f(242, 299)));
  segsAfter.push_back(Segment(Vec2f(240, 336), Vec2f(238, 351)));
  segsAfter.push_back(Segment(Vec2f(200, 339), Vec2f(282, 341)));
  segsAfter.push_back(Segment(Vec2f(140, 201), Vec2f(145, 319)));

  // This should monsterify Fredo a little.
  vector<Image> out = morph(werewolf, fredo, segsBefore, segsAfter, 30);

  // This is how you can write an image sequence
  for (int i = 0; i < (int)out.size(); ++i) {
    ostringstream fname;
    fname << "class_morph_";
    fname << setfill('0') << setw(2);
    fname << i;
    fname << ".png";
    out[i].write(fname.str());
  }
}

// This is a way for you to test your functions.
// We will only grade the contents of morphing.cpp and
// basicImageManipulation.cpp
int main() {
  cout << "nothing done in a5_main.cpp, debug me !" << endl;
  Image stevie("./Input/BostonRainbow-crop-400.png");
  Image output(scaleLin(stevie, .5));
  output.write("./Output/testNN.png");
  //testVectorOperations();
  //testSegment();
  // testNearestNeighbor();
  // testBilinearInterpolation();
  // testBilinearRescaling();
  // testRotation();
  // testVectorOperations();
  // testSegment();
  // testWarpBy1();
  // testWarp();
   testMorph();
}
