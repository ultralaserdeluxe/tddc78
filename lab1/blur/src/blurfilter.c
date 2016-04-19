/*
  File: blurfilter.c

  Implementation of blurfilter function.
    
 */
#include <stdio.h>
#include "blurfilter.h"
#include "ppmio.h"


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

void blurfilter(const int xsize, const int ysize, pixel* src, pixel* dst, const int radius, const double *w){
  int x,y,x2,y2, wx, wy;
  double r,g,b,n, wc;

  for (y=0; y<ysize; y++) {
    for (x=0; x<xsize; x++) {
      r = w[0] * pix(src, x, y, xsize)->r;
      g = w[0] * pix(src, x, y, xsize)->g;
      b = w[0] * pix(src, x, y, xsize)->b;
      n = w[0];
      for ( wx=1; wx <= radius; wx++) {
	for ( wy=1; wy <= radius; wy++) {
	  wc = w[wx] * w[wy];

	  x2 = x - wx;
	  y2 = y - wy;

	  if(x2 >= 0 || y2 >= 0) {
	    r += wc * pix(src, x2, y2, xsize)->r;
	    g += wc * pix(src, x2, y2, xsize)->g;
	    b += wc * pix(src, x2, y2, xsize)->b;
	    n += wc;
	  }
	  x2 = x + wx;
	  y2 = y + wy;
	  if(x2 < xsize || y2 < ysize) {
	    r += wc * pix(src, x2, y2, xsize)->r;
	    g += wc * pix(src, x2, y2, xsize)->g;
	    b += wc * pix(src, x2, y2, xsize)->b;
	    n += wc;
	  }
	}
      }
      pix(dst,x,y, xsize)->r = r/n;
      pix(dst,x,y, xsize)->g = g/n;
      pix(dst,x,y, xsize)->b = b/n;
    }
  }
}



