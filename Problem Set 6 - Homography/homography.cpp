#include "homography.h"
#include "matrix.h"

using namespace std;

void applyHomography(const Image& source, const Matrix& H, Image& out,
    bool bilinear) {
    // --------- HANDOUT  PS06 ------------------------------
    // Transform image source using the homography H, and composite in onto out.
    // if bilinear == true, using bilinear interpolation. Use nearest neighbor
    // otherwise.

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


    // Loop through all images blah 
    //check if x or y is within the boundaries of the source image
    for (int j = 0; j < height; j++) {
        for (int i= 0;i< width; i++) {

            outputVec = { float(i), float(j), 1.0 };
            sourceVec = invH * outputVec; //multiply output vector by H
            coordinate = {sourceVec(0)/sourceVec(2), sourceVec(1)/sourceVec(2)};
            //cout << coordinate << endl;

            //gets x'/w' and y'/w'
            if ((coordinate(0) >= 0) && (coordinate(0) < float(s_width)) && (coordinate(1) >= 0) && (coordinate(1) < float(s_height))) {
                for (int k = 0; k < channels; k++) {
                
                        if (bilinear == true) {
                            out(i, j, k) = interpolateLin(source, coordinate(0), coordinate(1), k);
                        }

                        else {
                            float ys = round(coordinate(1));
                            float xs = round(coordinate(0));
                            out(i, j, k) = source.smartAccessor(xs, ys, k, true);
                        }

                    }

                }


        }
    }

}

void subProcess(Matrix& mat, int space, CorrespondencePair correspondance_Line) {

    // fill in a through h....
    //axi + byi + c +0d +0e + 0f -xp*xi g - xp*yi h = 0 move to one side and fill in that way..
    int up1 = space + 1;

    float x1 = correspondance_Line.point1(0);
    float y1 = correspondance_Line.point1(1);

    float xp = correspondance_Line.point2(0); // x'
    float yp = correspondance_Line.point2(1); // y'

    mat(space, 0) = x1;
    mat(space, 1) = y1;
    mat(space, 2) = 1;
    mat(space, 3) = 0;
    mat(space, 4) = 0;
    mat(space, 5) = 0;
    mat(space, 6) = -x1 * xp; //
    mat(space, 7) = -y1 * xp;
    mat(space, 8) = -xp;
    

    mat(up1, 0) = 0;
    mat(up1, 1) = 0;
    mat(up1, 2) = 0;
    mat(up1, 3) = x1;
    mat(up1, 4) = y1;
    mat(up1, 5) = 1;
    mat(up1, 6) = -x1 * yp;
    mat(up1, 7) = -y1 * yp;
    mat(up1, 8) = -yp;

}

Matrix computeHomography(const CorrespondencePair correspondences[4]) {
  // --------- HANDOUT  PS06 ------------------------------
  // Compute a homography from 4 point correspondences.

    Matrix A = Matrix::Zero(9, 9);

    for (int i = 0; i < 4; i++) {
        subProcess(A, 2*i, correspondences[i]); // fill in each row by 2 so multiply the index by 2.. to ge two at once
    }

    A(8, 8) = 1;

    cout << "wait" << endl;
    // A is filled
    // B is ??? 0
    Matrix B = Matrix::Zero(9, 1); //said its (1,8) transposed...???
    B(8) = 1;

    //Matrix Bz = B.transposeInPlace();
    Matrix x = A.inverse() * B; //has to be in column form
   // for (int i = 0; i < x.size(); i++) {
    //    x(i) = x(i)/x()
   // }
    //for loop??
    cout << A << endl;
    cout << endl;

    Matrix Homography(3, 3);
    Homography(0, 0) = x(0);
    Homography(0, 1) = x(1);
    Homography(0, 2) = x(2);
    Homography(1, 0) = x(3);
    Homography(1, 1) = x(4);
    Homography(1, 2) = x(5);
    Homography(2, 0) = x(6);
    Homography(2, 1) = x(7);
    Homography(2, 2) = 1; //assumed it is one 

    return Homography;
}


BoundingBox computeTransformedBBox(int imwidth, int imheight, Matrix H) {
  // --------- HANDOUT  PS06 ------------------------------
  // Predict the bounding boxes that encompasses all the transformed
  // coordinates for pixels frow and Image with size (imwidth, imheight)
    
    Vec3f x1y1 = { 0,0,1 };
    Vec3f x2y1 = { float(imwidth) - 1, 0 , 1};
    Vec3f x1y2 = { 0, float(imheight) - 1, 1 };
    Vec3f x2y2 = { float(imwidth) - 1, float(imheight) - 1, 1 };

    Matrix tlc = H * x1y1;
    Matrix trc = H * x2y1;
    Matrix blc = H * x1y2;
    Matrix brc = H * x2y2;

    float xmin = 0;
    float xmax = 0;
    float ymax = 0;
    float ymin = 0;


    vector <Vec3f> points;
    Matrix x_vals = Matrix::Zero(4, 1);
    Matrix y_vals = Matrix::Zero(4, 1);

    points.push_back(tlc); points.push_back(trc); points.push_back(blc); points.push_back(brc);
    // Divid by w'
    for (int i = 0; i < 4; i++) {
        points[i](0) = points[i](0) / points[i](2);
        points[i](1) = points[i](1) / points[i](2);
        points[i](2) = 1;

        x_vals(i) = points[i](0);
        y_vals(i) = points[i](1);

    }
    xmin = x_vals.minCoeff();
    ymin = y_vals.minCoeff();
    xmax = x_vals.maxCoeff();
    ymax = y_vals.maxCoeff();


    return BoundingBox(xmin, xmax, ymin, ymax);
}

BoundingBox bboxUnion(BoundingBox B1, BoundingBox B2) {
    // --------- HANDOUT  PS06 ------------------------------
    // Compute the bounding box that tightly bounds the union of B1 an B2.

    float xWidth_max = max(B1.x2, B2.x2);
    float xWidth_min = min(B1.x1, B2.x1);
    float yWidth_max = max(B1.y2, B2.x2);
    float yWidth_min = min(B1.y1, B2.y1);


  return BoundingBox(xWidth_min, xWidth_max, yWidth_min, yWidth_max);
}

Matrix makeTranslation(BoundingBox B) {
  // --------- HANDOUT  PS06 ------------------------------
  // Compute a translation matrix (as a homography matrix) that translates the
  // top-left corner of B to (0,0).
    float tx = -B.x1;
    float ty = -B.y1;
    
    Matrix translation = Matrix::Identity(3, 3);
    translation(0, 2) = tx;
    translation(1, 2) = ty;  
    return translation;
}

Image stitch(const Image &im1, const Image &im2,
             const CorrespondencePair correspondences[4]) {
  // --------- HANDOUT  PS06 ------------------------------
  // Transform im1 to align with im2 according to the set of correspondences.
  // make sure the union of the bounding boxes for im2 and transformed_im1 is
  // translated properly (use makeTranslation)
    Matrix H = computeHomography(correspondences);
    Matrix invH = H.inverse();
    BoundingBox B = computeTransformedBBox(im1.width(), im1.height(), H);
    BoundingBox B2 = BoundingBox(0, im2.width() - 1, 0, im2.height()-1);
    BoundingBox BB = bboxUnion(B, B2);
    Matrix translation = makeTranslation(BB);
    
    //do homography then translate the image... so you can fill in the pixels
    int width = translation(0, 2) + BB.x2;
    int height = translation(1, 2) + BB.y2;
    Image output(width, height, im1.channels());

    // linear algebra is linear so we can multiply t and h ...
    Matrix HomTran = translation * H;
    applyHomography(im2, translation, output, true);
    applyHomography(im1, HomTran, output, true);


    return output;
}

// debug-useful
Image drawBoundingBox(const Image &im, BoundingBox bbox) {
  // --------- HANDOUT  PS06 ------------------------------
  /*
    ________________________________________
   / Draw me a bounding box!                \
   |                                        |
   | "I jumped to my                        |
   | feet, completely thunderstruck.space     |
   | blinked my eyes hard.spacelooked         |
   | carefully all around me. Andspacesaw a   |
   | most extraordinary small person, who   |
   | stood there examining me with great    |
   | seriousness."                          |
   \              Antoine de Saint-Exupery  /
    ----------------------------------------
           \   ^__^
            \  (oo)\_______
               (__)\       )\/\
                   ||----w |
                   ||     ||
  */
    Image output = im;
    for (int j = bbox.y1; j < bbox.y2 + 1; j++) {
        for (int i = bbox.x1; i < bbox.x2 + 1; i++) {
            for (int k = 0; k < im.channels(); k++) {
                output(i, j, k) = 0;
            }
            output(i, j, 0) = .4;
        }
    }
    //output.create_line(bbox.x1, bbox.x2, bbox.y1, bbox.y2, .5);
    return output;
    
    //return im;
}

void applyHomographyFast(const Image &source, const Matrix &H, Image &out,
                         bool bilinear) {
  // --------- HANDOUT  PS06 ------------------------------
  // Same as apply but change only the pixels of out that are within the
  // predicted bounding box (when H maps source to its new position).
}
