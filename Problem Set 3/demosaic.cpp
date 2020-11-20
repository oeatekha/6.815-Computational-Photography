/* --------------------------------------------------------------------------
 * File:    demosaic.cpp
 * Created: 2015-10-01
 * --------------------------------------------------------------------------
 *
 *
 *
 * ------------------------------------------------------------------------*/

#include "demosaic.h"
#include <cmath>

using namespace std;


Image basicGreen(const Image &raw, int offset) {
  // --------- HANDOUT  PS03 ------------------------------
  // Takes as input a raw image and returns a single-channel
  // 2D image corresponding to the green channel using simple interpolation
    //G 0 G 0 
    //0 G 0 G
    //G 0 G 0
    //0 G 0 G

    Image output(raw.width(), raw.height(), 1); //its still one channel of 3..
    int size = raw.number_of_elements();
    int innerboud_x = raw.width() - 1;
    int innerboud_y = raw.height() - 1;
    int offset_thisrow = offset;

    
    for (int j = 1; j < raw.height()-1; j++) {
        offset_thisrow = 1 - offset_thisrow;
            for (int i = 1; i < raw.width()-1; i++) {
                 //Lukas is clutch
                                                
                if (i % 2 == offset_thisrow) { //if first row of mosaic then we select it...                     
                        output(i, j, 0) = raw(i, j, 0);                                              
                }
                else {
                    output(i, j, 0) = ((raw(i, j - 1, 0) + raw(i, j + 1, 0) + raw(i + 1, j, 0) + raw(i - 1, j, 0)) / 4);
                }
                             
            }
           
        }
    
    return output;
}
// fix
Image basicRorB(const Image &raw, int offsetX, int offsetY) {
  // --------- HANDOUT  PS03 ------------------------------
  //  Takes as input a raw image and returns a single-channel
  //  2D image corresponding to the red or blue channel using simple
  //  interpolation
    int w = raw.width();
    int h = raw.height();
    int c = raw.channels();
    Image output(w, h, 1); //one channel r or b demosaic
    int offset_thisrow_x = offsetX;
    int offset_thisrow_y = offsetY;
    // So we have basically all the values we know then we have to interpolate the edges....using the one to the right or left..
    // Where are we? Is it a normal pixel? How is it averaged?
    //some thing is wrong with this....
    for (int j = offsetY ; j < h - 2; j = j + 2) {
        for (int i = offsetX; i < w - 2; i = i + 2) {
            
            output(i, j, 0) = raw(i, j, 0);  
            output(i + 1, j, 0) = (raw(i, j, 0) + raw(i + 2, j, 0)) / 2; // average of the anchor and 2+ the anchor to right
            output(i, j + 1, 0) = (raw(i, j, 0) + raw(i, j + 2, 0)) / 2; //go to down from the anchor and average it
            output(i + 1, j + 1, 0) = ((raw(i, j, 0) + raw(i, j + 2, 0)) + raw(i + 2, j, 0) + raw(i + 2, j + 2, 0)) / 4; //average of diagnolas

        }
    }
    
    return output;
} 

Image basicDemosaic(const Image &raw, int offsetGreen, int offsetRedX,
                    int offsetRedY, int offsetBlueX, int offsetBlueY) {
  // --------- HANDOUT  PS03 ------------------------------
  // takes as input a raw image and returns an rgb image
  // using simple interpolation to demosaic each of the channels
    vector <Image> rgb;
    rgb.push_back(basicRorB(raw, offsetRedX, offsetRedY)); //R
    rgb.push_back(basicGreen(raw, offsetGreen)); //G
    rgb.push_back(basicRorB(raw, offsetBlueX, offsetBlueY)); //B
    
    int width = raw.width();
    int height = raw.height();
    Image output(width, height, 3);

    
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            for (int k = 0; k < 3; k++) { //cycle through the array of the different image channels

                output(i, j, k) = rgb[k](i, j, 0);

            }
        }
    }
    
    
    return output;
}

float greenEdge(const Image& raw, int x, int y) {

    float output;
    float dy = abs(raw(x, y + 1, 0) - raw(x, y - 1, 0)); //up down scope //.subtract cause umm it is look for the updown diff
    float dx = abs(raw(x + 1 , y, 0) - raw(x - 1, y, 0)); //left right scope

    if (dx > dy) {
        output = (raw(x, y + 1, 0) + raw(x, y - 1, 0)) / 2.0; // if the difference between the up down is greater than the xy distance is the edge because they are more similar.. it is more gradual... which indicates that is the edge we should use
    }
    else {
        output = (raw(x + 1, y, 0) + raw(x - 1, y, 0)) / 2.0;
    }

    return output; 

}

Image edgeBasedGreen(const Image &raw, int offset) {
  // --------- HANDOUT  PS03 ------------------------------
  // Takes a raw image and outputs a single-channel
  // image corresponding to the green channel taking into account edges
    Image output(raw.width(), raw.height(), 1); //its still one channel of 3..
    int size = raw.number_of_elements();
    int innerboud_x = raw.width() - 1;
    int innerboud_y = raw.height() - 1;
    int offset_thisrow = offset;


    for (int j = 1; j < raw.height() - 1; j++) {
        offset_thisrow = 1 - offset_thisrow;
        for (int i = 1; i < raw.width() - 1; i++) {
             //Lukas is clutch
            if (i % 2 == offset_thisrow) { //if first row of mosaic then we select it...                     
                output(i, j, 0) = raw(i, j, 0);
            }
            else {
                output(i, j, 0) = greenEdge(raw, i, j);
            }

        }

    }

    return output;
}

Image edgeBasedGreenDemosaic(const Image &raw, int offsetGreen, int offsetRedX,
                             int offsetRedY, int offsetBlueX, int offsetBlueY) {
  // --------- HANDOUT  PS03 ------------------------------
  // Takes as input a raw image and returns an rgb image
  // using edge-based green demosaicing for the green channel and
  // simple interpolation to demosaic the red and blue channels
    vector <Image> rgb;
    rgb.push_back(basicRorB(raw, offsetRedX, offsetRedY)); //R
    rgb.push_back(edgeBasedGreen(raw, offsetGreen)); //G
    rgb.push_back(basicRorB(raw, offsetBlueX, offsetBlueY)); //B

    int width = raw.width();
    int height = raw.height();
    Image output(width, height, 3);


    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            for (int k = 0; k < 3; k++) { //cycle through the array of the different image channels

                output(i, j, k) = rgb[k](i, j, 0);

            }
        }
    }


    return output;
}

Image greenBasedRorB(const Image &raw, Image &green, int offsetX, int offsetY) {
    // --------- HANDOUT  PS03 ------------------------------
    // Takes as input a raw image and returns a single-channel
    // 2D image corresponding to the red or blue channel using green based
    // interpolation

    int w = raw.width();
    int h = raw.height();
    int c = raw.channels();
    Image output(w, h, 1); //one channel r or b demosaic
    int offset_thisrow_x = offsetX;
    int offset_thisrow_y = offsetY;
    // So we have basically all the values we know then we have to interpolate the edges....using the one to the right or left..
    // Where are we? Is it a normal pixel? How is it averaged?
    //some thing is wrong with this....
    for (int j = offsetY; j < h - 2; j = j + 2) {
        for (int i = offsetX; i < w - 2; i = i + 2) {

            output(i, j, 0) = raw(i, j, 0);
            output(i + 1, j, 0) = (raw(i, j, 0) + raw(i + 2, j, 0) - green(i,j,0) - green(i + 2, j, 0)) / 2 + green(i +1,j,0); // Do i subtract interpolate than add it back average of the anchor and 2+ the anchor to right
            output(i, j + 1, 0) = (raw(i, j, 0) + raw(i, j + 2, 0) - green(i, j, 0) - green(i, 2 + j, 0)) / 2 + green(i,j + 1,0); //go to down from the anchor and average it
            output(i + 1, j + 1, 0) = (raw(i, j, 0) + raw(i, j + 2, 0) + raw(i + 2, j, 0) + raw(i + 2, j + 2, 0) - green(i, j, 0) - green(i + 2, j, 0) - green(i, j +2, 0) -green(i + 2, j + 2, 0)) / 4 + green(i + 1,j +1,0); //average of diagnolas

        }
    }

    return output;
}

Image improvedDemosaic(const Image &raw, int offsetGreen, int offsetRedX,
                       int offsetRedY, int offsetBlueX, int offsetBlueY) {
    // // --------- HANDOUT  PS03 ------------------------------
    // Takes as input a raw image and returns an rgb image
    // using edge-based green demosaicing for the green channel and
    // simple green based demosaicing of the red and blue channels

    vector <Image> rgb;
    int width = raw.width();
    int height = raw.height();
    Image output(width, height, 3);
    Image green = edgeBasedGreen(raw, offsetGreen);


    rgb.push_back(greenBasedRorB(raw, green, offsetRedX, offsetRedY)); //R
    rgb.push_back(green); //G
    rgb.push_back(greenBasedRorB(raw, green, offsetBlueX, offsetBlueY)); //B

    


    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            for (int k = 0; k < 3; k++) { //cycle through the array of the different image channels

                output(i, j, k) = rgb[k](i, j, 0);

            }
        }
    }


    return output;
}
