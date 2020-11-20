// hdr.cpp
// Assignment 5

#include "hdr.h"
#include "filtering.h"
#include <algorithm>
#include <math.h>

using namespace std;

/**************************************************************
 //                       HDR MERGING                        //
 *************************************************************/

Image computeWeight(const Image &im, float epsilonMini, float epsilonMaxi) {
  // --------- HANDOUT  PS04 ------------------------------
  // Generate a weight image that indicates which pixels are good to use in
  // HDR, i.e. weight=1 when the pixel value is in [epsilonMini, epsilonMaxi].
  // The weight is per pixel, per channel.

  // Generate weights that are usable in the HDR combination. We want to ignore values that are either too bright or too noisy. 
  // These values clip and the values that are too dark are removed from the total. So we use binary weights

    int height = im.height();
    int width = im.width();
    int channels = im.channels();
    Image output(width, height, channels);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < channels; k++) {
                float crnt_Val = im(i, j, k);
                if (crnt_Val < epsilonMini || crnt_Val > epsilonMaxi) {
                    //If out the range either way ignore...
                    output(i, j, k) = 0;
                }
                else {
                    output(i, j, k) = 1;
                }
            }
        }
    }


    return output;
}

float computeFactor(const Image &im1, const Image &w1, const Image &im2,
                    const Image &w2) {
  // --------- HANDOUT  PS04 ------------------------------
  // Compute the multiplication factor between a pair of images. This
  // gives us the relative exposure between im1 and im2. It is computed as
  // the median of im2/(im1+eps) for some small eps, taking into account
  // pixels that are valid in both images.
  // gamma_code?? 
  // How do I get a median for an image..... I kinda get the ratio part but besides that..
  // wouldn't ratios be pretty similar often?
    int height = im1.height();
    int width = im1.width();
    int channels = im1.channels();
    vector <float> ratio_vals;
    //vector <float> sorted_vals;
    float eps = pow(10, -10); // small number
    float ratio = 0;


    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < channels; k++) {
                if (w1(i, j, k) == 1 && w2(i, j, k) == 1) {

                    float crnt_Val1 = im1(i, j, k);
                    float crnt_Val2 = im2(i, j, k); // im 2 value
                    ratio = crnt_Val2 / (crnt_Val1 + eps); 
                    ratio_vals.push_back(ratio); //push back the list....
                    //cout << endl << ratio;

                }
                     
            }
        }
    }
    sort(ratio_vals.begin(), ratio_vals.end()); //sort just reorders that mess
    //size of set../2 floor
    float k = ratio_vals[floor(ratio_vals.size()/2)];
    return k;
}

Image makeHDR(vector<Image> &imSeq, float epsilonMini, float epsilonMaxi) {
  // --------- HANDOUT  PS04 ------------------------------
  // Merge images to make a single hdr image
  // For each image in the sequence, compute the weight map (special cases
  // for the first and last images).
  // Compute the exposure factor for each consecutive pair of image.
  // Write the valid pixel to your hdr output, taking care of rescaling them
  // properly using the factor.

    //Dont forget to convert to gamma lol
    // its already inverted stupid...

    int height = imSeq[0].height();
    int width = imSeq[0].width();
    int channels = imSeq[0].channels();
    Image output(width, height, channels);
    vector <float> ratio;
    vector<Image> wSeq; //should be one 
    float total;
    float total_weights;
    float crrnt_ki = 0;
    float hdr_val = 0;

    //get Weights first val, in between, and last
    wSeq.push_back(computeWeight(imSeq[0], epsilonMini, 1)); //the light values are allowed in this context while epsilonMini is 0.001 or whatever
    for (int n = 1; n < imSeq.size()-1; n++) {
        wSeq.push_back(computeWeight(imSeq[n], epsilonMini, epsilonMaxi)); // get the weight for each image
    }
    wSeq.push_back(computeWeight(imSeq[(imSeq.size() - 1)], 0, epsilonMaxi)); //the dark values are allowed in this context while epsilonMaxi is .99 or whatever
    
    // get ratio ki for values
    ratio.push_back(1); //1 to 1 ratio with image 0......
    for (int n = 1; n < imSeq.size(); n++) {
        crrnt_ki = ratio[n-1] * (computeFactor(imSeq[n-1], wSeq[n-1], imSeq[n], wSeq[n])); // im1 then im2 ratio prev multiplied by adjacent..
        ratio.push_back(crrnt_ki);
    } 

    //Now perform equation...
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < channels; k++) { //this loop is honestly unecessary but for some reason I could not run my code without it so?? i dunno
                total = 0;
                total_weights = 0;

                for (int n = 0; n < imSeq.size(); n++) {

                    total = total + wSeq[n](i, j, k) * 1/ratio[n] * imSeq[n](i, j, k);
                    total_weights += wSeq[n](i, j, k); //tw = tw + w

                }
                
                hdr_val = 1 / total_weights * total; // 1/sumweights * sum(wi*1/ki*Ii(x,y);
                // We got to make sure the total weight isn't equal to 0. b/c we do 1/0...... and get a giant val so we check...
                if (total_weights == 0) { //parsing right?
                    output(i, j, k) = imSeq[0](i, j, k);
                }
                else {
                    output(i, j, k) = hdr_val;
                    //cout << hdr_val << endl;
                }

            }
        }
    }

    
    return output;
}

/**************************************************************
 //                      TONE MAPPING                        //
 *************************************************************/

Image toneMap(const Image &im, float targetBase, float detailAmp, bool useBila,
              float sigmaRange) {
  // --------- HANDOUT  PS04 ------------------------------
  // tone map an hdr image
  // - Split the image into its luminance-chrominance components.
  // - Work in the log10 domain for the luminance
  // -
    int height = im.height();
    int width = im.width();
    int channels = im.channels();
    float sigma = 0;
    vector <Image> LumChrom(lumiChromi(im));

    Image lumiLog(log10Image(LumChrom[0]));
    Image filtered(width, height, lumiLog.channels()); // error cause 3 channels//filtered image
    Image logOutput(width, height, lumiLog.channels());
    cout << endl << "kfactor";

    if (width > height) {  sigma = width / 50; }   else {  sigma = height / 50; }


    if (useBila == true) {
        
         filtered = bilateral(lumiLog, sigmaRange, sigma);
         cout << endl << "bilateral";
    }   
    else {
         
         filtered = gaussianBlur_separable(lumiLog, sigma);
         cout << endl << "gaussian";
    }
    

    Image detail = lumiLog - filtered;  
    float logLarge = lumiLog.max() - lumiLog.min(); //on log scale now so it should be really massive tbh
    float k_factor = log10(100) / logLarge;
    float lumimax = lumiLog.max();
    cout << endl << "details made";

    // I THOUGHT I COULD DO THIS : (
    // output = detailAmp*detail + k*(lumiLog - lumiLog.max());
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < lumiLog.channels(); k++) {
                logOutput(i, j, k) = detailAmp * detail(i, j, k) + k_factor*(lumiLog(i,j,k)-lumimax);
            }
        }
    }
    cout << endl <<" done with it";

    vector <Image> outputs;
    outputs.push_back(exp10Image(logOutput)); //convert lumi to exp10
    outputs.push_back(LumChrom[1]); //add chrominance back to the vector
    return lumiChromi2rgb(outputs);

}

/*********************************************************************
 *                       Tone mapping helpers                        *
 *********************************************************************/

// image --> log10Image
Image log10Image(const Image &im) {
  // --------- HANDOUT  PS04 ------------------------------
  // Taking a linear image im, transform to log10 scale.
  // To avoid infinity issues, make any 0-valued pixel be equal the the
  // minimum non-zero value. See image_minnonzero(im).
    int height = im.height();
    int width = im.width();
    int channels = im.channels();
    Image output(width, height, channels);
    float nonZero = image_minnonzero(im);
    // Not for this part
    // float largeRange = im.max() - image_minnonzero(im); // do i do in log range?
    // float k = log10(100) / largeRange;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < channels; k++) {
                if (im(i, j, k) == 0) {
                    output(i, j, k) = log10(nonZero);
                }
                else {
                    output(i, j, k) = log10(im(i, j, k));
                }
            }
        }
    }

    return output;
}

// Image --> 10^Image
Image exp10Image(const Image &im) {
  // --------- HANDOUT  PS04 ------------------------------
  // take an image in log10 domain and transform it back to linear domain.
  // see pow(a, b)
    int height = im.height();
    int width = im.width();
    int channels = im.channels();
    Image output(width, height, channels);

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            for (int k = 0; k < channels; k++) {
                output(i,j,k) = pow(10, im(i, j, k)); //base 10^y = x ............log10(x) =  y
            }
        }
    }

    return output;

}

// min non-zero pixel value of image
float image_minnonzero(const Image &im) {
  // --------- HANDOUT  PS04 ------------------------------
  // return the smallest value in the image that is non-zeros (across all
  // channels too)
    float minimum = 2147483647; // Big number or whatever

    for (int j = 0; j < im.height(); j++) {
        for (int i = 0; i < im.width(); i++) {
            for (int k = 0; k < im.channels(); k++) {
                if (im(i, j, k) != 0 && im(i, j, k) < minimum) { //ignores when it is equal to 0 anything but 0 can be checked though
                    minimum = im(i, j, k);
                }
                
            }

        }
    }
  return minimum;
}
