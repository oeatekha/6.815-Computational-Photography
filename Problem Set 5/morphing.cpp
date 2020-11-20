/* -----------------------------------------------------------------
 * File:    morphing.cpp
 * Created: 2015-09-25
 * -----------------------------------------------------------------
 *
 *
 *
 * ---------------------------------------------------------------*/

#include "morphing.h"
#include <cassert>

using namespace std;

Vec2f operator+(const Vec2f &a, const Vec2f &b) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return the vector sum of a an b
    float x_val = a.x + b.x;
    float y_val = a.y + b.y;
    return Vec2f(x_val, y_val); 
}

Vec2f operator-(const Vec2f &a, const Vec2f &b) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a-b
    float x_val = b.x - a.x;
    float y_val = b.y - a.y;
    return Vec2f(x_val, y_val);
}

Vec2f operator*(const Vec2f &a, float f) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a*f
    float x_val = a.x * f;
    float y_val = a.y * f;

  return Vec2f(x_val, y_val); 
}

Vec2f operator/(const Vec2f &a, float f) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a/f
    float x_val = a.x / f;
    float y_val = a.y / f;

    return Vec2f(x_val, y_val);
}

float dot(const Vec2f &a, const Vec2f &b) {
  // --------- HANDOUT  PS05 ------------------------------
    float product = a.x * b.x + a.y * b.y;
    return product; // change me
}

float length(const Vec2f &a) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return the length of a.
    float euclidian_length = sqrt(pow(a.x, 2) + pow(a.y, 2));
    return euclidian_length; // change me
}

Vec2f perpendicular(const Vec2f &a) {
  // --------- HANDOUT  PS05 ------------------------------
  // Return a vector that is perpendicular to a.
  // Either direction is fine.
    Vec2f normal = { -a.y, a.x };   
    return normal;
}

// The Segment constructor takes in 2 points P(x1,y1) and Q(x2,y2) corresponding
// to the ends of a segment and initialize the local reference frame e1,e2.
Segment::Segment(Vec2f P_, Vec2f Q_) : P(P_), Q(Q_) {
  // // --------- HANDOUT  PS05 ------------------------------
  // // The initializer list above ": P(P_), Q(Q_)" already copies P_
  // // and Q_, so you don't have to do it in the body of the constructor.
  // You should:
  // * Initialize the local frame e1,e2 (see header file)
   // this->lPQ = lPQ;
    
    Vec2f PQ = P - Q; //??
    lPQ = length(PQ);
    e1 = PQ / lPQ;
    e2 = perpendicular(e1);
    //  e1 = (Q-P)/(||Q-P||) in the paper
    //  e2 = perpendicalur(Q-P)/(||Q-P||) in the paper
}

Vec2f Segment::XtoUV(Vec2f X) const {
  // --------- HANDOUT  PS05 ------------------------------
  // Compute the u,v coordinates of a point X with
  // respect to the local frame of the segment.
  // e2 ^
  //    |
  // v  +  * X
  //    | /
  //    |/
  //    *--+------>-----*
  //    P  u     e1     Q
  //                    u=1
  //
  // * Be careful with the different normalization for u and v
  //  Segment PX(P, X); //segment from p to u 
    //e1 PQ/lPQ is the unit vector u equals PQ/lPQ^2 while v is simply scalled by lPQ.
    //e1 just equals PX * e1/lPQ = PX * PQ/(lPQ^2)
    Segment PX(P, X); 
    Vec2f PX_v = P - X;
    float u = dot(PX_v, e1)/lPQ; //e1 projects that stuff and dividing by lPQ makes it less than 1.. dot products are commutative
    float v = dot(PX_v, e2); //v is not rescaled u is within the x local distance
    return Vec2f(u, v); // changeme
}

Vec2f Segment::UVtoX(Vec2f uv) const {
  // --------- HANDOUT  PS05 ------------------------------
  // compute the (x, y) position of a point given by the (u,v)
  // location relative to this segment.
  // * Be careful with the different normalization for u and v
  // X = P + u*PQ + v * perpindicularPQ / |PQ|
    Vec2f PQ = P - Q;
    Vec2f X = P + PQ * uv.x + (perpendicular(PQ) / lPQ) * uv.y; //order * correctly
    return Vec2f(X);
}

float Segment::distance(Vec2f X) const {
  // --------- HANDOUT  PS05 ------------------------------
  // Implement distance from a point X(x,y) to the segment. Remember the 3
  // cases from class.
    // top vertex
    // bottom vertex
    // inbetween line
    Vec2f uv = XtoUV(X); //uv vector.... u is the distance from P 
    Vec2f euclidian = { 0,0 };
    float output = 0;
    //cout << "u is " << uv.x << endl;
    //cout << "v is " << uv.y << endl;

    if (uv.x > 0 && uv.x < 1) { 
        
        //if inbetween P&Q
        output = abs(uv.y); //perpindicular distance of the segment PQ to the point in the e2 direction..
        return output;
    }
    if (uv.x > 1) { 

        //if past Q
        Segment QX(Q, X);
        return QX.lPQ;
        //euclidian = Q - X;
        //output = length(euclidian);
    }
    else { 

        //if near P
        Segment PX(P, X);
        return PX.lPQ;
        //euclidian = P - X;
        //output = length(euclidian);
    }
}

Image warpBy1(const Image &im, const Segment &segBefore,
              const Segment &segAfter) {
  // --------- HANDOUT  PS05 ------------------------------
  // Warp an entire image according to a pair of segments.
    segBefore;
    segAfter;
    int width = im.width();
    int height = im.height();
    int channels = im.channels();
    Vec2f X;
    Vec2f X_dest;
    Vec2f uv_corresponding;
    Image output(width, height, channels);
    float interpVal = 0;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            
               X = { float(i),float(j) }; //convert
               uv_corresponding = segAfter.XtoUV(X);
               X_dest = segBefore.UVtoX(uv_corresponding); // X location in the source image
                //use it for the new image...
                //bilinear??scaleLin

            for (int k = 0; k < channels; k++) {
               interpVal = interpolateLin(im, X_dest.x, X_dest.y, k, true);
               output(i, j, k) = interpVal;
            }
        }
    }


    return output;
}

float Segment::weight(Vec2f X, float a, float bs, float p) const {
  // --------- HANDOUT  PS05 ------------------------------
  // compute the weight of a segment to a point X(x,y) given the weight
  // parameters a,b, and p (see paper for details).

    float weight  = pow(pow(lPQ, p)/(a + distance(X)), 2);
    return weight; // changeme
}

Image warp(const Image &im, const vector<Segment> &src_segs,
           const vector<Segment> &dst_segs, float a, float b, float p) {
  // --------- HANDOUT  PS05 ------------------------------
  // Warp an image according to a vector of before and after segments using
  // segment weighting
    int lines = dst_segs.size(); //n segments to distort and move....
    int height = im.height();
    int width = im.width();
    int channels = im.channels();

    Image output(width, height, channels);

    Vec2f X_p; //prime
    Vec2f X_i; // source 
    Vec2f X_final;
    Vec2f DSUM;
    Vec2f uv_corresponding;
    float interpVal = 0;
    float dist = 0;
    float weight;
    output = im;
    Vec2f D_i = { 0, 0 };
    float weightsum = 0;

    
        // output = warpBy1(output, src_segs[n], dst_segs[n]);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            DSUM = { 0,0 };
            weightsum = 0;

            for (int n = 0; n < lines; n++) {
                X_p = { float(i),float(j) }; //convert
                uv_corresponding = dst_segs[n].XtoUV(X_p);
                X_i = src_segs[n].UVtoX(uv_corresponding); // X location in the destination
                D_i = X_p - X_i; // is this order correct
                weight = dst_segs[n].weight(X_i, a, b, p);
                DSUM = DSUM + (D_i * weight);
                weightsum = weight + weightsum;
            }
            X_final = X_p + DSUM / weightsum;
            //now interpolate
            for (int k = 0; k < channels; k++) {
                interpVal = interpolateLin(im, X_final.x, X_final.y, k, true);
                output(i, j, k) = interpVal;
            }
        }
    }

    return output;
}

vector<Image> morph(const Image &im_before, const Image &im_after,
                    const vector<Segment> &segs_before,
                    const vector<Segment> &segs_after, int N, float a, float b,
                    float p) {
  // --------- HANDOUT  PS05 ------------------------------
  // return a vector of N+2 images: the two inputs plus N images that morphs
  // between im_before and im_after for the corresponding segments. im_before
  // should be the first image, im_after the last.
    if (im_before.width() != im_after.width() || im_before.height() != im_after.height()) {
        throw MismatchedDimensionsException();
    }

    vector<Image> imSeq;
    imSeq.push_back(im_before);
    //n segments to distort and move....
    int height = im_before.height();
    int width = im_before.width();
    int channels = im_before.channels();
    Image output(width, height, channels);
   
    //vector<int> t;
    //t.push_back(0);
    cout << "N is " << N << endl;


    float tc = 1/float(N + 1); //t change every sequence...
    //cout << "tc is " << tc << endl;

    for (int n = 1; n < N+1; n++) {
        
        float t = tc * n; // each index has a larger t value starting from 0
        float delta = 1 - t;

        vector<Segment> InterpSeg;
        //cout << "t is " << t << endl;

        for (int l = 0; l < segs_before.size(); l++) {
            //interpolate the edges accrding to paper
            
            
            Vec2f seg1 = segs_before[l].getP()*delta + segs_after[l].getP()*t;
            Vec2f seg2 = segs_before[l].getQ()*delta + segs_after[l].getQ()*t;
            Segment newSeg(seg1, seg2);
            InterpSeg.push_back(newSeg);
        }

        Image output_0 = warp(im_before, segs_before, InterpSeg);
        Image output_1 = warp(im_after, segs_after, InterpSeg);

        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                for (int k = 0; k < channels; k++) {

                    output(i, j, k) = delta * output_0(i,j,k)  + t * output_1(i,j,k);
                    //you gotta linear interpolate again dont forger

                }
            
            }
        }
        imSeq.push_back(output);

    }






    imSeq.push_back(im_after);
    return imSeq;
}
