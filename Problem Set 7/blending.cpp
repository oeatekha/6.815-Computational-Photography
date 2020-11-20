#include "blending.h"
#include "matrix.h"

#include <ctime>

using namespace std;

Image blendingweight(int imwidth, int imheight) {
  // // --------- HANDOUT  PS07 ------------------------------
    Image output(imwidth, imheight, 1);
    float weight_x = 0;
    float weight_y = 0;
    float weight = 0;

    for (int j = 0; j < imheight; j++) {
        for (int i = 0; i < imwidth; i++) {
            weight_x = 1 - abs(float(i) - imwidth / 2) / (imwidth / 2); //get the absolute value of the x location from the center divided by the unit length
            weight_y = 1 - abs(float(j) - imheight / 2) / (imheight / 2);
            weight = weight_x * weight_y;
            output(i, j, 0) = weight;
        }
    }

    return output;
}

//  ****************************************************************************
//  * blending related functions re-written from previous assignments
//  ****************************************************************************

// instead of writing source in out, *add* the source to out based on the weight
// so out(x,y) = out(x, y) + weight * image
void applyhomographyBlend(const Image &source, const Image &weight, Image &out,
                          const Matrix &H, bool bilinear) {
  // --------- HANDOUT  PS07 ------------------------------
    Vec3f sourceVec; // is x' y' z'
    Vec3f outputVec; //multuplied by H  is x y z 
    Vec2f coordinate;

    int width = out.width();
    int height = out.height();
    int channels = out.channels();

    int s_width = source.width();
    int s_height = source.height();
    int s_channels = source.channels();
    Matrix invH = H.inverse();
    Image weight_Source = source * copychannels(weight, source.channels());


    // Loop through all images blah 
    //check if x or y is within the boundaries of the source image
    for (int j = 0; j < height; j++) {
        for (int i = 0;i < width; i++) {

            outputVec = { float(i), float(j), 1.0 };
            sourceVec = invH * outputVec; //multiply output vector by H
            coordinate = { sourceVec(0) / sourceVec(2), sourceVec(1) / sourceVec(2) };
            //cout << coordinate << endl;

            //gets x'/w' and y'/w'
            if ((coordinate(0) >= 0) && (coordinate(0) < float(s_width)) && (coordinate(1) >= 0) && (coordinate(1) < float(s_height))) {
                for (int k = 0; k < channels; k++) {

                    if (bilinear == true) {
                        out(i, j, k) = out(i, j, k) +  interpolateLin(weight_Source, coordinate(0), coordinate(1), k); //out(x, y) + weight * image
                    }

                    else {
                        float ys = round(coordinate(1));
                        float xs = round(coordinate(0));
                        out(i, j, k) = out(i, j, k) + weight_Source.smartAccessor(xs, ys, k, true);
                    }

                }

            }


        }
    }

}

Image stitchLinearBlending(const Image &im1, const Image &im2, const Image &we1,
                           const Image &we2, const Matrix &H) {
  // --------- HANDOUT  PS07 ------------------------------
  // stitch using image weights.
  // note there is no weight normalization.
    

    BoundingBox B1 = computeTransformedBBox(im1.width(), im1.height(), H);
    BoundingBox B2 = BoundingBox(0, im2.width(), 0, im2.height());
    BoundingBox BB = bboxUnion(B1, B2);
    Matrix T = makeTranslation(BB);
    int width = T(0, 2) + BB.x2;
    int height = T(1, 2) + BB.y2;

    Image output(width, height, im1.channels());
    Matrix HT = T * H;

    applyhomographyBlend(im1, we1, output, HT, true);
    applyhomographyBlend(im2, we2, output, T, true);
    

    return output;

}

/*****************************************************************************
 * blending functions Pset 08
 *****************************************************************************/

// low freq and high freq (2-scale decomposition)
vector<Image> scaledecomp(const Image &im, float sigma) {
  vector<Image> ims;
  ims.push_back(gaussianBlur_separable(im, sigma));
  ims.push_back(im - ims[0]);
  return ims;
}

// stitch using different blending models
// blend can be 0 (none), 1 (linear) or 2 (2-layer)
Image stitchBlending(const Image &im1, const Image &im2, const Matrix &H,
                     BlendType blend) {
  // --------- HANDOUT  PS07 ------------------------------
    //BlendType { BLEND_NONE, BLEND_LINEAR, BLEND_2LAYER }; 

    BoundingBox B1 = computeTransformedBBox(im1.width(), im1.height(), H);
    BoundingBox B2 = BoundingBox(0, im2.width(), 0, im2.height());
    BoundingBox BB = bboxUnion(B1, B2);
    Matrix T = makeTranslation(BB);
    int width = T(0, 2) + BB.x2;
    int height = T(1, 2) + BB.y2;
    Image we1 = blendingweight(im1.width(), im1.height());
    Image we2 = blendingweight(im2.width(), im2.height());

    Image output(width, height, im1.channels());
    Image output2(width, height, im1.channels());

    Image linearW(width, height, 1); //linear place holder for weights
    Matrix HT = T * H;
    // switch case; if or else too much
    Image lin_we1(im1.width(), im1.height(), im1.channels());
    Image lin_we2(im2.width(), im2.height(), im2.channels());
    lin_we1 = lin_we1 + 1; lin_we2 = lin_we2 + 1;

    switch (blend) {
    
        case BlendType(0): { //BLEND NONE

            applyHomography(im1, HT, output, true);
            applyHomography(im2, T, output, true);
            return output;

        }
        case BlendType(1): { // BLEND LINEAR

            output = stitchLinearBlending(im1, im2, we1, we2, H);
            
            linearW = stitchLinearBlending(lin_we1, lin_we2, we1, we2, H);
            cout << linearW.channels() << " channels" << endl;
            return output / (linearW + pow(10,-15));      //cant divide by zero    //Normalize output....
        }
        case BlendType(2): {

            vector <Image> lowhi1 = scaledecomp(im1, 2.0);
            vector <Image> lowhi2 = scaledecomp(im2, 2.0);

            Image lowF1 = lowhi1[0];
            Image lowF2 = lowhi2[0];
            Image hiF1 =  lowhi1[1];
            Image hiF2 =  lowhi2[1];     
            Image layer_lowf = stitchBlending(lowF1, lowF2, H, BlendType::BLEND_LINEAR);


            
            //Upper Layer aka Details
            Image tempW1(width, height, 1);
            Image tempW2(width, height, 1);
            //get the two images with the respective weights.... so you can compare which one to keep....
            applyHomography(we1, HT, tempW1, true);
            applyHomography(we2, T, tempW2, true);



            Image hiftemp1(width, height, im1.channels());
            Image hiftemp2(width, height, im1.channels());
            //get the two images with the respective weights.... so you can compare which one to keep....
            applyHomography(hiF1, HT, hiftemp1, true);
            applyHomography(hiF2, T, hiftemp2, true); //gives me the two hi frequecny transformations....
            Image hifW1(width, height, 1);
            Image hifW2(width, height, 1);


            for (int j = 0; j < layer_lowf.height(); j++) {
                for (int i = 0; i < layer_lowf.width(); i++) {
                    
                        if (tempW1(i, j) > tempW2(i, j)) {
                            for (int k = 0; k < layer_lowf.channels(); k++) {
                                output(i, j, k) = hiftemp1(i, j, k);
                                //hifW1(i, j, 0) = 1;
                                //hifW2(i, j, 0) = 0;
                            }

                        }
                        else {
                            for (int k = 0; k < layer_lowf.channels(); k++) {
                                output(i, j, k) = hiftemp2(i, j, k);
                                //hifW1(i, j, 0) = 0;
                                //hifW2(i, j, 0) = 1;
                            }
                        }

                    
                }
            }
            //cout << " navi says hey " << endl;
            //output = stitchLinearBlending(hiF1, hiF2, hifW1, hifW2, H); 
            

            // does not work cout << " navi says hey " << endl;
            return  layer_lowf + output;

            }

    }




    return output ;
}

// auto stitch
Image autostitch(const Image &im1, const Image &im2, BlendType blend,
                 float blurDescriptor, float radiusDescriptor) {
  // --------- HANDOUT  PS07 ------------------------------
    vector <Feature> feat1 = computeFeatures(im1, HarrisCorners(im1), blurDescriptor, radiusDescriptor);
    vector <Feature> feat2 = computeFeatures(im2, HarrisCorners(im2), blurDescriptor, radiusDescriptor);
    vector <FeatureCorrespondence> featcoro = findCorrespondences(feat1, feat2);
    Matrix H = RANSAC(featcoro);
    Image output = stitchBlending(im1, im2, H, blend);

    return output;
}

/************************************************************************
 * Tricks: mini planets.
 ************************************************************************/

Image pano2planet(const Image &pano, int newImSize, bool clamp) {
  // // --------- HANDOUT  PS07 ------------------------------
    Image pano360(newImSize, newImSize, pano.channels());

    int center = float(newImSize) / 2; //center of x and y coordinates
    float radius, x_sqr, y_sqr;
    float theta;
    float x_width, y_height;

    for (int j = 0; j < newImSize; j++) {
        for (int i = 0; i < newImSize; i++) {
            x_sqr = float(i) - center;
            y_sqr = float(j) - center; //get the distance from the center from the right, left, up, and down
            radius = sqrt(pow(x_sqr, 2) + pow(y_sqr, 2));
            theta = atan2(y_sqr, x_sqr); //matlab lol returns -pi to pi.....
            y_height = pano.height() * (1 - (radius / center)); 

            //atan2 is weird. only range values from -pi to pi accrdng to luke........
            // theta less than zero is just gonna be the similar to y height because its clockwise.. or the left side...
            if (theta < 0) {
                x_width = pano.width() * (theta / 3.14159265) * -0.5; //-0.5 to deal with the negative theta... and also to only go to half of the width.
            }
            else {
                x_width = pano.width()* (1 - (theta / 3.14159265)) * 0.5 + pano.width()/2; //start from mid of panorama width....
            }
            //cout << "hereh";
            //interpolate..
            for (int k = 0; k < pano.channels(); k++) {
                pano360(i, j, k) = interpolateLin(pano, x_width, y_height, k, true);
            }

        }
    }

    return pano360;

}

/************************************************************************
 * 6.865: Stitch N images into a panorama
 ************************************************************************/

// Pset07-865. Compute sequence of N-1 homographies going from Im_i to Im_{i+1}
// Implement me!
vector<Matrix> sequenceHs(const vector<Image> &ims, float blurDescriptor,
                          float radiusDescriptor) {
  // // --------- HANDOUT  PS07 ------------------------------
  return vector<Matrix>();
}

// stack homographies:
//   transform a list of (N-1) homographies that go from I_i to I_i+1
//   to a list of N homographies going from I_i to I_refIndex.
vector<Matrix> stackHomographies(const vector<Matrix> &Hs, int refIndex) {
  // // --------- HANDOUT  PS07 ------------------------------
  return vector<Matrix>();
}

// Pset07-865: compute bbox around N images given one main reference.
BoundingBox bboxN(const vector<Matrix> &Hs, const vector<Image> &ims) {
  // // --------- HANDOUT  PS07 ------------------------------
  return BoundingBox(0, 0, 0, 0);
}

// Pset07-865.
// Implement me!
Image autostitchN(const vector<Image> &ims, int refIndex, float blurDescriptor,
                  float radiusDescriptor) {
  // // --------- HANDOUT  PS07 ------------------------------
  return Image(1, 1, 1);
}

/******************************************************************************
 * Helper functions
 *****************************************************************************/

// copy a single-channeled image to several channels
Image copychannels(const Image &im, int nChannels) {
  // image must have one channel
  assert(im.channels() == 1);
  Image oim(im.width(), im.height(), nChannels);

  for (int i = 0; i < im.width(); i++) {
    for (int j = 0; j < im.height(); j++) {
      for (int c = 0; c < nChannels; c++) {
        oim(i, j, c) = im(i, j);
      }
    }
  }
  return oim;
}
