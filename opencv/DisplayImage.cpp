#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <stdlib.h>
#include <stdio.h>

using namespace cv;

/// Global variables

Mat src, src_gray;
Mat dst, detected_edges;
Mat dest_gray;

int edgeThresh = 1;
int lowThreshold = 40;
int const max_lowThreshold = 50;
int ratio = 3;
int kernel_size = 3;
char window_name[] = "Edge Map";

#define MAX_BRIGHTNESS 255

/**
 * @function CannyThreshold
 * @brief Trackbar callback - Canny thresholds input with a ratio 1:3
 */
void CannyThreshold(int, void*)
{
  /// Reduce noise with a kernel 5x5
  blur( src_gray, detected_edges, Size(5,5) );

  /// Canny detector
  Canny( detected_edges, detected_edges, lowThreshold, 50, kernel_size );

  /// Using Canny's output as a mask, we display our result
  dst = Scalar::all(0);

  src_gray.copyTo( dst, detected_edges);

  cvtColor(dst, dest_gray, CV_GRAY2BGR);

  for (int i = 0; i < dest_gray.rows; i++) {
    for (int j = 0; j < dest_gray.cols; j++) {

      if (dest_gray.at<cv::Vec3b>(i,j)[0] > 0) {
        dest_gray.at<cv::Vec3b>(i,j)[0] = MAX_BRIGHTNESS;
      }

      if (dest_gray.at<cv::Vec3b>(i,j)[1] > 0) {
        dest_gray.at<cv::Vec3b>(i,j)[1] = MAX_BRIGHTNESS;
      }

     if (dest_gray.at<cv::Vec3b>(i,j)[2] > 0) {
        dest_gray.at<cv::Vec3b>(i,j)[2] = MAX_BRIGHTNESS;
      }
    }
  }

  imshow( window_name, dest_gray );
 }


/** @function main */
int main( int argc, char** argv )
{
  /// Load an image
  src = imread( argv[1] );

  if( !src.data )
  { return -1; }

  /// Create a matrix of the same type and size as src (for dst)
  dst.create( src.size(), src.type() );

  /// Convert the image to grayscale
  cvtColor( src, src_gray, CV_BGR2GRAY );

  /// Create a window
  namedWindow( window_name, CV_WINDOW_AUTOSIZE );

  /// Create a Trackbar for user to enter threshold
  createTrackbar( "Min Threshold:", window_name, &lowThreshold, max_lowThreshold, CannyThreshold );

  /// Show the image
  CannyThreshold(0, 0);

  /// Wait until user exit program by pressing a key
  waitKey(0);

  imwrite("opencv_out.bmp", dest_gray);


  return 0;
  }