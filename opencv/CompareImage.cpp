#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <stdlib.h>
#include <stdio.h>

using namespace cv;

/// Global variables

Mat src_ref, src_comp;

/** @function main */
int main( int argc, char** argv )
{

  if (argc != 3) {
    perror("Number of arguments invalid: Try ./<exec> <image_ref> <image_comp>");
    exit(-1);
  }

  /// Load an image
  src_ref = imread( argv[1] );
  src_comp = imread( argv[2] );

  if( !src_ref.data )
  { return -1; }

  if( !src_comp.data )
  { return -1; }


  int rows = src_ref.rows;
  int cols = src_ref.cols;

  double total_pixels = rows*1.0F * cols*1.0F * 3.0;
  double diff_pixels = 0.0F;

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {

      if (src_ref.at<cv::Vec3b>(i,j)[0] !=  src_comp.at<cv::Vec3b>(i,j)[0]) {
	//printf("[%d][%d] Ref %d != Comp %d \n",i, j, src_ref.at<cv::Vec3b>(i,j)[0], src_comp.at<cv::Vec3b>(i,j)[0]);
        diff_pixels += 1.0F;
      }

      if (src_ref.at<cv::Vec3b>(i,j)[1] !=  src_comp.at<cv::Vec3b>(i,j)[1]) {
	//printf("[%d][%d] Ref %d != Comp %d \n",i, j, src_ref.at<cv::Vec3b>(i,j)[1], src_comp.at<cv::Vec3b>(i,j)[1]);
        diff_pixels += 1.0F;
      }

      if (src_ref.at<cv::Vec3b>(i,j)[2] !=  src_comp.at<cv::Vec3b>(i,j)[2]) {
	//printf("[%d][%d] Ref %d != Comp %d \n",i, j, src_ref.at<cv::Vec3b>(i,j)[2], src_comp.at<cv::Vec3b>(i,j)[2]);
        diff_pixels += 1.0F;
      }
    }
  }

  double diff_percentage = diff_pixels / total_pixels;

  printf("The difference percentage between images is %.2lf \n", diff_percentage * 100.0F);

  return 0;
}
