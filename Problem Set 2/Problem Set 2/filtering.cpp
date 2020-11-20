/* -----------------------------------------------------------------
 * File:    filtering.cpp
 * Created: 2015-09-22
 * -----------------------------------------------------------------
 *
 * Image convolution and filtering
 *
 * ---------------------------------------------------------------*/

#include "filtering.h"
#include <cassert>
#include <cmath>

using namespace std;

Image boxBlur(const Image &im, int k, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Convolve an image with a box filter of size k by k
  // It is safe to asssume k is odd.
    int center = (k - 1) / 2;
    float total = 0;
    Image blured(im.width(), im.height(), im.channels());
    for (int i = 0; i < im.width(); i++) {
        for (int j = 0; j < im.height(); j++) {
            
            for (int m = 0; m < im.channels(); m++) {
                total = 0;
                //im(i, j, m);
                for (int x = -center; x <= center; x++) {
                    for (int y = -center; y <= center; y++) {

                        total = total + im.smartAccessor(i + x, j + y, m, clamp);
                        
                    }
                }

               // cout << endl << total / (pow(center * 2 + 1, 2));
                blured(i, j, m) = total / (pow(center * 2 + 1, 2)); //divide by num of elements in k square

            }
        }
    }
    blured.write("bihblur.png");
    return blured; // change this
}

Image Filter::convolve(const Image &im, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Write a convolution function for the filter class
    int w_center = (width - 1) / 2;
    int h_center = (height - 1) / 2;

    float total = 0;
    float temp_val = 0;
    Image output(im.width(), im.height(), im.channels());
    for (int i = 0; i < im.width(); i++) {
        for (int j = 0; j < im.height(); j++) {
            for (int k = 0; k < im.channels(); k++) {
                total = 0;

                for (int x = 0; x < width; x++) {
                    for (int y = 0; y < height; y++) {

                        temp_val = im.smartAccessor(i + x -w_center, j + y - h_center, k, clamp);
                        total = total + temp_val * (*this)(width - 1 - x, height - 1 - y); //take the current value then multiply it by the x,y coordin of kernal
                        //cout << "check" << endl;
                                                                                           //this got messy doing it in reverse bec I wrote mine frontwards
                        // so you have to flip and adjust for backwards
                        // okay remeber to subtrack 1 and to flip kernels. width is starting from 0...
                    }
                }

                output(i, j, k) = total; //divide by num of elements in k square

            }
        }
    }

  return output; // change this
}

Image boxBlur_filterClass(const Image &im, int k, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Reimplement the box filter using the filter class.
  // check that your results match those in the previous function "boxBlur"
    vector<float> kernel; //create a row major kernel 
    Image output(im.width(), im.height(), im.channels());

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) { 

            kernel.push_back(1 / pow(k, 2)); // all values are the same for this

        }
    }

    Filter f(kernel, k, k); //creates kernel with 10 10 that can access the x,y position in row major order.
    output = f.convolve(im, clamp);
    //output.write("bihbih.png");
    return output; // change this
}

Image gradientMagnitude(const Image &im, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Uses a Sobel kernel to compute the horizontal and vertical components
  // of the gradient of an image and returns the gradient magnitude.
    vector<float> Sobel_H = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    vector<float> Sobel_V = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    Image output(im.width(), im.height(), im.channels());
    Image hsk_filtered(im.width(), im.height(), im.channels());
    Image vsk_filtered(im.width(), im.height(), im.channels());

    Filter h(Sobel_H, 3, 3);
    Filter v(Sobel_V, 3, 3);

    hsk_filtered = h.convolve(im, clamp);
    vsk_filtered = v.convolve(im, clamp);

    float didx = 0;
    float didy = 0;
    float temp = 0;

    for (int i = 0; i < im.width(); i++) {
        for (int j = 0; j < im.height(); j++) {
            for (int k = 0; k < im.channels(); k++) {
               
                didx = hsk_filtered(i, j, k);
                didy = vsk_filtered(i, j, k);
                temp = sqrt(pow(didx, 2) + pow(didy, 2));
                output(i, j, k) = temp;
            }
        }
    }

    output.write("gradient_magnitude.png");
    return output; // change this
}

vector<float> gauss1DFilterValues(float sigma, float truncate) {
  // --------- HANDOUT  PS02 ------------------------------
  // Create a vector containing the normalized values in a 1D Gaussian filter
  // Truncate the gaussian at truncate*sigma.
    float normilazation = 1 / sqrt(2 * atan(1)*4 * sigma * sigma); //atan(1)*4 == pi
    float temp_total = 0;
    int limit = ceil(sigma * truncate); // returns int
    float k_val = 0;
    vector<float> kernel;

    for (int i = -limit; i <= limit; i++) {      
        temp_total = temp_total + exp((-i*i)/(2*sigma*sigma));
    }

    for (int i = -limit; i <= limit; i++) {

        k_val = exp((-i * i) / (2 * sigma * sigma)) / temp_total;
        kernel.push_back(k_val);

    }

    return kernel;
}

Image gaussianBlur_horizontal(const Image &im, float sigma, float truncate,
                              bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Gaussian blur across the rows of an image
    Image output(im.width(), im.height(), im.channels());
    vector<float> kernel = gauss1DFilterValues(sigma, truncate);
    Filter k(kernel, kernel.size(), 1);

    output = k.convolve(im, clamp);
    //output.write("gaussian.png");
    return output;
}

vector<float> gauss2DFilterValues(float sigma, float truncate) {
  // --------- HANDOUT  PS02 ------------------------------
  // Create a vector containing the normalized values in a 2D Gaussian
  // filter. Truncate the gaussian at truncate*sigma.

    float temp_total = 0;
    int limit = ceil(sigma * truncate); // returns int
    float k_val;
    vector<float> kernel;

    // Get normalization...
    for (int i = -limit; i <= limit; i++) {
        for (int j = -limit; j <= limit; j++) {
            temp_total = temp_total + exp(-(i*i + j*j)/ (2 * sigma*sigma)); //x^2 + y^2
        }
    }

    // do the same thingy...
    for (int i = -limit; i <= limit; i++) {
        for (int j = -limit; j <= limit; j++) {

            k_val = exp(-(i * i + j * j) / (2 * sigma * sigma)) / temp_total; //x^2 + y^2 /sigma^2 /Total
            kernel.push_back(k_val);
        }
    }

  return kernel;
}

Image gaussianBlur_2D(const Image &im, float sigma, float truncate,
                      bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Blur an image with a full 2D rotationally symmetric Gaussian kernel
    Image output(im.width(), im.height(), im.channels());
    vector<float> kernel = gauss2DFilterValues(sigma, truncate);
    float kernel_l = 2 * sigma * truncate + 1; //kernel length... like x and y...
    Filter g(kernel, kernel_l, kernel_l); //kernel size gives the total amount of elements you need row..

    output = g.convolve(im, clamp);
    //output.write("gaussian_2d.png");
    return output;
}

Image gaussianBlur_separable(const Image &im, float sigma, float truncate,
                             bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Use principles of seperabiltity to blur an image using 2 1D Gaussian
  // Filters
    Image output(im.width(), im.height(), im.channels());
    Image horizontal(im.width(), im.height(), im.channels());
    vector<float> kernel = gauss1DFilterValues(sigma, truncate);
    float kernel_l = 2 * sigma * truncate + 1;

    Filter h_conv(kernel, kernel_l, 1);
    Filter v_conv(kernel, 1, kernel_l);

    horizontal = h_conv.convolve(im, clamp);
    output = v_conv.convolve(horizontal, clamp);
    
    //output.write("gaussian_seperated.png");
    return output;
}

Image unsharpMask(const Image &im, float sigma, float truncate, float strength,
                  bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Sharpen an image
    Image output(im.width(), im.height(), im.channels());
    Image gaussian(im.width(), im.height(), im.channels());
    gaussian = gaussianBlur_separable(im, sigma);
    output = im + (strength *(im - gaussian));
    //output.write("sharpedn.png");
    return output;

}

Image bilateral(const Image &im, float sigmaRange, float sigmaDomain,
                float truncateDomain, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // Denoise an image using the bilateral filter

    //compute rgb distance...
    //compute xyz difference 
    float point_distance = 0;
    float val = 0;
    float tempWeight = 0;
    int limit = ceil(sigmaDomain * truncateDomain); // returns int
    float total;
    float weight;
    float difSqrd;
    float x_c, y_c, z_c;
    float x_n, y_n, z_n;
    float denominator = 2.0 * sigmaRange * sigmaRange;

    Image output(im.width(), im.height(), im.channels());  
    vector<float> kernel = gauss1DFilterValues(sigmaDomain, truncateDomain);; //make kernel...



    for (int i = 0; i < im.width(); i++) {
        for (int j = 0; j < im.height(); j++) {
            for (int k = 0; k < im.channels(); k++) {


                total = 0;
                val = 0;
                //easier to conceptualize if we go from -limit unlike in gaus...func
                for (int x = -limit; x <= limit; x++) {
                    for (int y = -limit; y <= limit; y++) {
                        
                        z_c = im(i, j, 2);
                        z_n = im.smartAccessor(i + x, j + y, 2, clamp);
                        y_c = im(i, j, 1);
                        y_n = im.smartAccessor(i + x, j + y, 1, clamp);
                        x_c = im(i, j, 0);  
                        x_n = im.smartAccessor(i + x ,j + y, 0, clamp);                        
                        
                        
                        difSqrd = pow(x_n - x_c, 2) + pow(y_n - y_c, 2) + pow(z_n - z_c, 2); //x^2 + y^2 + z^2
                        tempWeight =  kernel[limit + x] * kernel[limit + y] * exp(-difSqrd/denominator); //(G(x) * G(Y) * G(rgb_dist/sigma_r^2)...
                        total = total + tempWeight; //temp = weight total weights...
                        val = val +  tempWeight * im.smartAccessor(i + x, j + y, k, clamp); //total weights multiplied by the different im kernel locations

                        //idea is that the kernel is basically 1d symmetrical so we can just use the 1d gaus
                        //add x bc we start at -limit...so its flipped so go from there and..move x++
                        //so pull the kernel[index +x] for th especific G(x-x').. or y..
                        // so you have to flip and adjust for backwards
                        // okay remeber to subtrack 1 and to flip kernels. width is starting from 0...
                    }
                }

                output(i, j, k) = val / total; 
                // sum/sumweights                                                      
                //sqrt( x^2 + y^2 +z^2) should be (x-x2)^2  
                // I dunno add to gaussian?
            }
        }
    }

    output.write("bilateral.png");
    return output;
}

Image bilaYUV(const Image &im, float sigmaRange, float sigmaY, float sigmaUV,
              float truncateDomain, bool clamp) {
  // --------- HANDOUT  PS02 ------------------------------
  // 6.865 only
  // Bilaterial Filter an image seperatly for
  // the Y and UV components of an image
  return im;
}

/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// Create an image of 0's with a value of 1 in the middle. This function
// can be used to test that you have properly set the kernel values in your
// Filter object. Make sure to set k to be larger than the size of your kernel
Image impulseImg(int k) {
  // initlize a kxkx1 image of all 0's
  Image impulse(k, k, 1);

  // set the center pixel to have intensity 1
  int center = floor(k / 2);
  impulse(center, center, 0) = 1.0f;

  return impulse;
}

// ------------- FILTER CLASS -----------------------
Filter::Filter(const vector<float> &fData, int fWidth, int fHeight)
    : kernel(fData), width(fWidth), height(fHeight) {
  assert(fWidth * fHeight == (int)fData.size());
}

Filter::Filter(int fWidth, int fHeight)
    : kernel(std::vector<float>(fWidth * fHeight, 0)), width(fWidth),
      height(fHeight) {}

Filter::~Filter() {}

const float &Filter::operator()(int x, int y) const {
  if (x < 0 || x >= width)
    throw OutOfBoundsException();
  if (y < 0 || y >= height)
    throw OutOfBoundsException();

  return kernel[x + y * width];
}

float &Filter::operator()(int x, int y) {
  if (x < 0 || x >= width)
    throw OutOfBoundsException();
  if (y < 0 || y >= height)
    throw OutOfBoundsException();

  return kernel[x + y * width];
}
// --------- END FILTER CLASS -----------------------
 