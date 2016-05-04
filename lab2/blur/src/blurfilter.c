/*
  File: blurfilter.c

  Implementation of blurfilter function.
    
 */
#include <stdio.h>
#include <stdlib.h>
#include "blurfilter.h"


pixel* pix(pixel* image, const int xx, const int yy, const int xsize)
{
  register int off = xsize*yy + xx;

#ifdef DBG
  if(off >= MAX_PIXELS) {
    fprintf(stderr, "\n Terribly wrong: %d %d %d\n",xx,yy,xsize);
  }
#endif
  return (image + off);
}

void blurfilter_x(const int ystart, const int yend,  const int xsize, const int ysize, pixel* src, pixel* dst, const int radius, const double *w){
  int y;
  for (y=ystart; y<yend; y++) {
    int x;
    for (x=0; x<xsize; x++) {
      double r = w[0] * pix(src, x, y, xsize)->r;
      double g = w[0] * pix(src, x, y, xsize)->g;
      double b = w[0] * pix(src, x, y, xsize)->b;
      double n = w[0];

      int wi;
      for (wi=1; wi <= radius; wi++) {
  	double wc = w[wi];
  	int x2 = x - wi;
  	if(x2 >= 0) {
  	  r += wc * pix(src, x2, y, xsize)->r;
  	  g += wc * pix(src, x2, y, xsize)->g;
  	  b += wc * pix(src, x2, y, xsize)->b;
  	  n += wc;
  	}
  	x2 = x + wi;
  	if(x2 < xsize) {
  	  r += wc * pix(src, x2, y, xsize)->r;
  	  g += wc * pix(src, x2, y, xsize)->g;
  	  b += wc * pix(src, x2, y, xsize)->b;
  	  n += wc;
  	}
      }
      pix(dst,x,y, xsize)->r = r/n;
      pix(dst,x,y, xsize)->g = g/n;
      pix(dst,x,y, xsize)->b = b/n;
    }
  }
}

void blurfilter_y(const int ystart, const int yend,  const int xsize, const int ysize, pixel* src, pixel* dst, const int radius, const double *w){
  int y;
  for (y=ystart; y<yend; y++) {
    int x;
    for (x=0; x<xsize; x++) {
      double r = w[0] * pix(src, x, y, xsize)->r;
      double g = w[0] * pix(src, x, y, xsize)->g;
      double b = w[0] * pix(src, x, y, xsize)->b;
      double n = w[0];

      int wi;
      for ( wi=1; wi <= radius; wi++) {
  	double wc = w[wi];
  	
	int y2 = y - wi;
  	if(y2 >= 0) {
  	  r += wc * pix(src, x, y2, xsize)->r;
  	  g += wc * pix(src, x, y2, xsize)->g;
  	  b += wc * pix(src, x, y2, xsize)->b;
  	  n += wc;
  	}
  	y2 = y + wi;
  	if(y2 < ysize) {
  	  r += wc * pix(src, x, y2, xsize)->r;
  	  g += wc * pix(src, x, y2, xsize)->g;
  	  b += wc * pix(src, x, y2, xsize)->b;
  	  n += wc;
  	}
      }
      pix(dst,x,y, xsize)->r = r/n;
      pix(dst,x,y, xsize)->g = g/n;
      pix(dst,x,y, xsize)->b = b/n;
    }
  }
}
