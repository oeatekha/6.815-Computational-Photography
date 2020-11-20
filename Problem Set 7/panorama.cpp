#include "panorama.h"
#include "matrix.h"
#include <cassert>
#include <ctime>
#include <unistd.h>

using namespace std;

Image computeTensor(const Image &im, float sigmaG, float factorSigma) {
  // // --------- HANDOUT  PS07 ------------------------------
  // Compute xx/xy/yy Tensor of an image. (stored in that order)
    int width = im.width();
    int height = im.height();
    int channels = im.channels();
    Image output(width, height, channels);
    Image tensors(width, height, channels);


    Image lumi = gaussianBlur_2D(lumiChromi(im)[0], sigmaG); //blurred luminance or black and white version of the image...
    Image luminanceX = gradientX(lumi);
    Image luminanceY = gradientY(lumi);
    float sigmaFinal = sigmaG * factorSigma;

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {

            tensors(i, j, 0) = pow(luminanceX(i,j,0), 2);
            tensors(i, j, 1) = luminanceX(i,j,0) * luminanceY(i,j,0);
            tensors(i, j, 2) = pow(luminanceY(i, j, 0), 2);
        }
    }
    output = gaussianBlur_separable(tensors, sigmaFinal);
    return output;
}

Image cornerResponse(const Image &im, float k, float sigmaG,
                     float factorSigma) {
  // // --------- HANDOUT  PS07 ------------------------------
  // Compute response = det(M) - k*[(trace(M)^2)] at every pixel location,
  // using the structure tensor of im.
    int width = im.width();
    int height = im.height();
    int channels = im.channels();
    Image R_output(width, height, 1);
    Image tensor = computeTensor(im, sigmaG, factorSigma);
    Matrix M = Matrix::Zero(2, 2);
    float Ixx, Iyy, Ixy;

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            
            //fill in Matrix
            Ixx = tensor(i, j, 0); Ixy = tensor(i, j, 1); Iyy = tensor(i, j, 2);
            M(0, 0) = Ixx; M(0, 1) = Ixy; M(1, 0) = Ixy; M(1, 1) = Iyy;
            float R = M.determinant() - k*pow(M.trace(), 2);
            R_output(i, j, 0) = R;

        }
    }

    //R_output = maximum_filter(R_output, 7);
    return R_output;
}

vector<Point> HarrisCorners(const Image& im, float k, float sigmaG,
    float factorSigma, float maxiDiam,
    float boundarySize) {
    // // --------- HANDOUT  PS07 ------------------------------
    // Compute Harris Corners by maximum filtering the cornerResponse map.
    // The corners are the local maxima.
    int width = im.width();
    int height = im.height();
    int channels = im.channels();
    float feature = 0;
    vector<Point> list;

    Image corners = cornerResponse(im, k, sigmaG, factorSigma);
    Image maxf = maximum_filter(corners, maxiDiam);



    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            feature = corners(i, j, 0);
            // Only positive features 
            // Cant be near the edge
            // check if it is the local maxima by comparing it to regions where everything is maxiumum...
            if (feature > 0 && feature == maxf(i,j,0)) {
                if (i > boundarySize && j > boundarySize && i < width - boundarySize && j < height - boundarySize) {
                    Point position(i, j);
                    list.push_back(position);
                }
            }

        }
    }

    return list;
}

Image descriptor(const Image &blurredIm, const Point &p,
                 float radiusDescriptor) {
  // --------- HANDOUT  PS07 ------------------------------
  // Extract a descriptor from blurredIm around point p, with a radius
  // 'radiusDescriptor'.
    int width = blurredIm.width();
    int height = blurredIm.height();
    int channels = blurredIm.channels();
    float y = 0; float x = 0;
    float radius = radiusDescriptor * 2 + 1;
    
    //get values
    Image descriptor(radius, radius, 1);

    for (int j = p.y - radiusDescriptor; j <= p.y + radiusDescriptor; j++) {
        x = 0;
        for (int i = p.x - radiusDescriptor; i <= p.x + radiusDescriptor; i++) {
            descriptor(x,y,0) = blurredIm(i, j, 0);
            x++;
        }
        y++;
    }

    float avg = descriptor.mean();
    float std = sqrt(descriptor.var());
    for (int j = 0; j < descriptor.height(); j++) {
        for (int i = 0; i < descriptor.width(); i++) {
            descriptor(i,j,0) = (descriptor(i, j, 0) - avg) / std;
        }
    }

    return descriptor;
}

vector<Feature> computeFeatures(const Image &im, const vector<Point> &cornersL,
                                float sigmaBlurDescriptor,
                                float radiusDescriptor) {
  // // --------- HANDOUT  PS07 ------------------------------
  // Pset07. obtain corner features from a list of corner points
    int width = im.width();
    int height = im.height();
    int channels = im.channels();

    Image blurredIm = gaussianBlur_2D(im, sigmaBlurDescriptor);
    vector <Feature> listFeature;
    
    for (int i = 0; i < cornersL.size(); i++) {
        Image desc = descriptor(blurredIm, cornersL[i], radiusDescriptor);
        Feature feat(cornersL[i], desc);
        listFeature.push_back(feat);
        
    }

  return listFeature;

}

float l2Features(const Feature &f1, const Feature &f2) {
  // // --------- HANDOUT  PS07 ------------------------------
  // Compute the squared Euclidean distance between the descriptors of f1, f2.
  Image point1 =  f1.desc(); Image point2 = f2.desc();
  float dist = 0;
  float total = 0;

  for (int j = 0; j < point1.height(); j++) {
    for (int i = 0; i < point1.width(); i++) {

      dist = point1(i,j,0) - point2(i,j,0);
      dist = pow(dist, 2);
      total = total + dist;

    }
  }

  return total; // return squared total...
}

vector<FeatureCorrespondence> findCorrespondences(const vector<Feature> &listFeatures1,
 const vector<Feature> &listFeatures2, float threshold) {
  // // --------- HANDOUT  PS07 ------------------------------
  // Find correspondences between listFeatures1 and listFeatures2 using the
  // second-best test.
    //return min distance....
    
    vector<FeatureCorrespondence> corr;

    float dist_temp;
    float min_dist;
    float snd_dist;
    float ratio = 0;
    Feature crrntfeature2 = listFeatures2[0];
    Feature most_similar = crrntfeature2;
    Feature feature1 = listFeatures1[0];
    FeatureCorrespondence coro = FeatureCorrespondence(feature1, crrntfeature2);
    //threshold = pow(threshold, 2); am i suppsoed to square this ???

    //cout << "where" << endl;
    for (int i = 0; i < listFeatures1.size(); i++) { //if i swith the feature 2 to 1 the correspondances are different
        
        feature1 = listFeatures1[i];
        min_dist = 47483647;
        snd_dist = 2147483647; 
  
        for (int j = 0; j < listFeatures2.size(); j++) {
            
            crrntfeature2 = listFeatures2[j];
            dist_temp = l2Features(feature1, crrntfeature2); //get new temporary distance..
            
            
            if (dist_temp < min_dist) {
                snd_dist = min_dist;// do this first or it wont work lol second current 
                min_dist = dist_temp;

                most_similar = crrntfeature2; //at this moment current feature is the most similar feature cause the squared distance is the lowest..
               // cout << "min dist is " << min_dist << endl;

            }
            else if(dist_temp < snd_dist){
                snd_dist = dist_temp; // if the new one is only lower than the second shortest distance replace it...
                //cout << "snd dist is " << snd_dist << endl;
            }
        }

        ratio = snd_dist / min_dist;      
        if (ratio > pow(threshold, 2)) {
            coro = FeatureCorrespondence(feature1, most_similar);
            corr.push_back(coro);
        }

    }

    return corr;

}

vector<bool> inliers(const Matrix &H,
                     const vector<FeatureCorrespondence> &listOfCorrespondences,
                     float epsilon) {
  // // --------- HANDOUT  PS07 ------------------------------
  // Pset07: Implement as part of RANSAC
  // return a vector of bools the same size as listOfCorrespondences
  // indicating whether each correspondance is an inlier according to the
  // homography H and threshold epsilon
    vector<bool> inlier;
    
    //vector<CorrespondencePair> rand4Sample;
    //randomSample = sampleFeatureCorrespondences(listOfCorrespondences);
    //rand4Sample = { randomSample[0].toCorrespondencePair(), randomSample[1].toCorrespondencePair(), randomSample[2].toCorrespondencePair(), randomSample[3].toCorrespondencePair() }; //4 random samples top four in random order..
    float product = 0;
    Vec3f pointI;
    float Euclidian;
    Vec3f delta;
    
    for (int i = 0; i < listOfCorrespondences.size(); i++) {

        pointI = H * listOfCorrespondences[i].toCorrespondencePair().point1; //so many dots lol
        // Normalize by dividing by z...
        pointI[0] = pointI[0] / pointI[2];
        pointI[1] = pointI[1] / pointI[2];
        
        delta = pointI - listOfCorrespondences[i].toCorrespondencePair().point2; //subtract adjusted vec
        //cout << delta << endl;

        Euclidian = sqrt(pow(delta[0], 2) + pow(delta[1], 2));
        if (Euclidian < epsilon) { //if it is leess than epi^2 then its cool
            inlier.push_back(true);           
        }
        else {
            inlier.push_back(false);
        }
        

    }

    return inlier;
}

Matrix RANSAC(const vector<FeatureCorrespondence> &listOfCorrespondences,
              int Niter, float epsilon) {
  // // --------- HANDOUT  PS07 ------------------------------
  // Put together the RANSAC algorithm.
    vector<FeatureCorrespondence> randomSample = listOfCorrespondences;
    vector<CorrespondencePair> correspondances;
    vector<bool> inlier;
    vector<bool> inlierMain;

    Matrix H;
    Matrix finalH;
    int total = 0;
    int max_total =  0;


    Matrix I = Matrix::Identity(3,3); //creates identity matrix to return if H is singular or if det(H) = 0..


    for (int i = 0; i < Niter; i++) {
        //cout << i << endl;
        //new random sample each time...
        total = 0;
        randomSample = sampleFeatureCorrespondences(listOfCorrespondences);
        correspondances = getListOfPairs(randomSample);       
        H = computeHomography(correspondances.data()); //data has to be used to retrive actual values
        
        if (H.determinant() == 0) {
            H = I;
        }

        inlier = inliers(H, listOfCorrespondences, epsilon);
        
        for (int j = 0; j < listOfCorrespondences.size(); j++) {
            if (inlier[j] == 1) {
                total = total + 1;
            }
            
        }

        //total = countBoolVec(inlier); //you dont have to use a for loop
        //cout << total << endl;

        if (total > max_total) {

            max_total = total; //gives the the inlier with the highest total...
            finalH = H;

        }
       // cout << H << endl;

    }
    
    return finalH;
}

Image autostitch(const Image &im1, const Image &im2, float blurDescriptor,
                 float radiusDescriptor) {
  // // --------- HANDOUT  PS07 ------------------------------
  // Now you have all the ingredients to make great panoramas without using a
  // primitive javascript UI !
    vector <Feature> feat1 = computeFeatures(im1, HarrisCorners(im1), blurDescriptor, radiusDescriptor);
    vector <Feature> feat2 = computeFeatures(im2, HarrisCorners(im2), blurDescriptor, radiusDescriptor);
    vector <FeatureCorrespondence> featcoro = findCorrespondences(feat1, feat2);
    Matrix H = RANSAC(featcoro);

    BoundingBox B = computeTransformedBBox(im1.width(), im1.height(), H);
    BoundingBox B2 = BoundingBox(0, im2.width(), 0, im2.height());
    BoundingBox BB = bboxUnion(B, B2);
    Matrix translation = makeTranslation(BB);

    //do homography then translate the image... so you can fill in the pixels
    int width = translation(0, 2) + BB.x2;
    int height = translation(1, 2) + BB.y2;
    Image output(width, height, im1.channels());

    // linear algebra is linear so we can multiply t and h ...
    Matrix HomTran = translation * H;
    applyHomography(im1, HomTran, output, true);
    applyHomography(im2, translation, output, true);
    


    return output;
}

// *****************************************************************************
//  * Helpful optional functions to implement
// ****************************************************************************

Image getBlurredLumi(const Image &im, float sigmaG) { return Image(1, 1, 1); }

int countBoolVec(const vector<bool> &ins) { return 0; }

// *****************************************************************************
//  * Do Not Modify Below This Point
// *****************************************************************************

// Pset07 RANsac helper. re-shuffle a list of correspondances
vector<FeatureCorrespondence> sampleFeatureCorrespondences(
    vector<FeatureCorrespondence> listOfCorrespondences) {
  random_shuffle(listOfCorrespondences.begin(), listOfCorrespondences.end());
  return listOfCorrespondences;
}

// Pset07 RANsac helper: go from 4 correspondances to a list of points [4][2][3]
// as used in Pset06.
// Note: The function uses the first 4 correspondences passed
vector<CorrespondencePair>
getListOfPairs(const vector<FeatureCorrespondence> &listOfCorrespondences) {
  assert(listOfCorrespondences.size() >= 4);
  vector<CorrespondencePair> out;
  for (int i = 0; i < 4; i++) {
    out.push_back(listOfCorrespondences[i].toCorrespondencePair());
  }
  return out;
}

// Corner visualization.
Image visualizeCorners(const Image &im, const vector<Point> &pts, int rad,
                       const vector<float> &color) {
  Image vim = im;
  for (int i = 0; i < (int)pts.size(); i++) {
    int px = pts[i].x;
    int py = pts[i].y;

    int minx = max(px - rad, 0);

    for (int delx = minx; delx < min(im.width(), px + rad + 1); delx++)
      for (int dely = max(py - rad, 0); dely < min(im.height(), py + rad + 1);
           dely++) {
        if (sqrt(pow(delx - px, 2) + pow(dely - py, 2)) <= rad) {
          for (int c = 0; c < im.channels(); c++) {
            vim(delx, dely, c) = color[c];
          }
        }
      }
  }
  return vim;
}

Image visualizeFeatures(const Image &im, const vector<Feature> &LF,
                        float radiusDescriptor) {
  // assumes desc are within image range
  Image vim = im;
  int rad = radiusDescriptor;

  for (int i = 0; i < (int)LF.size(); i++) {
    int px = LF[i].point().x;
    int py = LF[i].point().y;
    Image desc = LF[i].desc();

    for (int delx = px - rad; delx < px + rad + 1; delx++) {
      for (int dely = py - rad; dely < py + rad + 1; dely++) {
        vim(delx, dely, 0) = 0;
        vim(delx, dely, 1) = 0;
        vim(delx, dely, 2) = 0;

        if (desc(delx - (px - rad), dely - (py - rad)) > 0) {
          vim(delx, dely, 1) = 1;
        } else if (desc(delx - (px - rad), dely - (py - rad)) < 0) {
          vim(delx, dely, 0) = 1;
        }
      }
    }
  }
  return vim;
}

void drawLine(const Point &p1, const Point &p2, Image &im,
              const vector<float> &color) {
  float minx = min(p1.x, p2.x);
  float miny = min(p1.y, p2.y);
  float maxx = max(p1.x, p2.x);
  float maxy = max(p1.y, p2.y);

  int spaces = 1000;
  for (int i = 0; i < spaces; i++) {
    float x = minx + (maxx - minx) / spaces * (i + 1);
    float y = miny + (maxy - miny) / spaces * (i + 1);
    for (int c = 0; c < im.channels(); c++) {
      im(x, y, c) = color[c];
    }
  }
}

Image visualizePairs(const Image &im1, const Image &im2,
                     const vector<FeatureCorrespondence> &corr) {
  Image vim(im1.width() + im2.width(), im1.height(), im1.channels());

  // stack the images
  for (int j = 0; j < im1.height(); j++) {
    for (int c = 0; c < im1.channels(); c++) {
      for (int i = 0; i < im1.width(); i++) {
        vim(i, j, c) = im1(i, j, c);
      }
      for (int i = 0; i < im2.width(); i++) {
        vim(i + im1.width(), j, c) = im2(i, j, c);
      }
    }
  }

  // draw lines
  for (int i = 0; i < (int)corr.size(); i++) {
    Point p1 = corr[i].feature(0).point();
    Point p2 = corr[i].feature(1).point();
    p2.x = p2.x + im1.width();
    drawLine(p1, p2, vim);
  }
  return vim;
}

Image visualizePairsWithInliers(const Image &im1, const Image &im2,
                                const vector<FeatureCorrespondence> &corr,
                                const vector<bool> &ins) {
  Image vim(im1.width() + im2.width(), im1.height(), im1.channels());

  // stack the images
  for (int j = 0; j < im1.height(); j++) {
    for (int c = 0; c < im1.channels(); c++) {
      for (int i = 0; i < im1.width(); i++) {
        vim(i, j, c) = im1(i, j, c);
      }
      for (int i = 0; i < im2.width(); i++) {
        vim(i + im1.width(), j, c) = im2(i, j, c);
      }
    }
  }

  // draw lines
  vector<float> red(3, 0);
  vector<float> green(3, 0);
  red[0] = 1.0f;
  green[1] = 1.0f;

  for (int i = 0; i < (int)corr.size(); i++) {
    Point p1 = corr[i].feature(0).point();
    Point p2 = corr[i].feature(1).point();
    p2.x = p2.x + im1.width();
    if (ins[i]) {
      drawLine(p1, p2, vim, green);
    } else {
      drawLine(p1, p2, vim, red);
    }
  }
  return vim;
}

// Inliers:  Detected corners are in green, reprojected ones are in red
// Outliers: Detected corners are in yellow, reprojected ones are in blue
vector<Image> visualizeReprojection(const Image &im1, const Image &im2,
                                    const Matrix &H,
                                    const vector<FeatureCorrespondence> &corr,
                                    const vector<bool> &ins) {
  // Initialize colors
  vector<float> red(3, 0);
  vector<float> green(3, 0);
  vector<float> blue(3, 0);
  vector<float> yellow(3, 0);
  red[0] = 1.0f;
  green[1] = 1.0f;
  blue[2] = 1.0f;
  yellow[0] = 1.0f;
  yellow[1] = 1.0f;

  vector<Point> detectedPts1In;
  vector<Point> projectedPts1In;
  vector<Point> detectedPts1Out;
  vector<Point> projectedPts1Out;

  vector<Point> detectedPts2In;
  vector<Point> projectedPts2In;
  vector<Point> detectedPts2Out;
  vector<Point> projectedPts2Out;

  for (int i = 0; i < (int)corr.size(); i++) {
    Point pt1 = corr[i].feature(0).point();
    Point pt2 = corr[i].feature(1).point();
    Matrix P1 = pt1.toHomogenousCoords();
    Matrix P2 = pt2.toHomogenousCoords();
    Matrix P2_proj = H * P1;
    Matrix P1_proj = H.inverse() * P2;
    Point reproj1 = Point(P1_proj(0) / P1_proj(2), P1_proj(1) / P1_proj(2));
    Point reproj2 = Point(P2_proj(0) / P2_proj(2), P2_proj(1) / P2_proj(2));
    if (ins[i]) { // Inlier
      detectedPts1In.push_back(pt1);
      projectedPts1In.push_back(reproj1);
      detectedPts2In.push_back(pt2);
      projectedPts2In.push_back(reproj2);
    } else { // Outlier
      detectedPts1Out.push_back(pt1);
      projectedPts1Out.push_back(reproj1);
      detectedPts2Out.push_back(pt2);
      projectedPts2Out.push_back(reproj2);
    }
  }

  vector<Image> output;
  Image vim1(im1);
  Image vim2(im2);
  vim1 = visualizeCorners(im1, detectedPts1In, 2, green);
  vim1 = visualizeCorners(vim1, projectedPts1In, 1, red);
  vim1 = visualizeCorners(vim1, detectedPts1Out, 2, yellow);
  vim1 = visualizeCorners(vim1, projectedPts1Out, 1, blue);

  vim2 = visualizeCorners(im2, detectedPts2In, 2, green);
  vim2 = visualizeCorners(vim2, projectedPts2In, 1, red);
  vim2 = visualizeCorners(vim2, detectedPts2Out, 2, yellow);
  vim2 = visualizeCorners(vim2, projectedPts2Out, 1, blue);

  output.push_back(vim1);
  output.push_back(vim2);
  return output;
}

/***********************************************************************
 * Point and Feature Definitions *
 **********************************************************************/
Point::Point(int xp, int yp) : x(xp), y(yp) {}

Point::Point() : x(0.0f), y(0.0f) {}

void Point::print() const { printf("(%d, %d)\n", x, y); }

Vec3f Point::toHomogenousCoords() const { return Vec3f(x, y, 1.0f); }

// Feature Constructors
Feature::Feature(const Point &ptp, const Image &descp) : pt(ptp), dsc(descp) {
  pt = ptp;
  dsc = descp;
}

// getter functions
const Point &Feature::point() const { return pt; }
const Image &Feature::desc() const { return dsc; }

// printer
void Feature::print() const {
  printf("Feature:");
  point().print();
  for (int j = 0; j < dsc.height(); j++) {
    for (int i = 0; i < dsc.width(); i++) {
      printf("%+07.2f ", dsc(i, j));
    }
    printf("\n");
  }
}

// FeatureCorrespondence Constructors
FeatureCorrespondence::FeatureCorrespondence(const Feature &f0p,
                                             const Feature &f1p)
    : f0(f0p), f1(f1p) {}

vector<Feature> FeatureCorrespondence::features() const {
  return vector<Feature>{f0, f1};
}

const Feature &FeatureCorrespondence::feature(int i) const {
  assert(i >= 0 && i < 2);
  if (i == 0)
    return f0;
  else
    return f1;
}

// printer
void FeatureCorrespondence::print() const {
  printf("FeatureCorrespondence:");
  f0.print();
  f1.print();
}

CorrespondencePair FeatureCorrespondence::toCorrespondencePair() const {
  return CorrespondencePair((float)f0.point().x, (float)f0.point().y, 1,
                            (float)f1.point().x, (float)f1.point().y, 1);
}
