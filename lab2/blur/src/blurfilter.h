/*
  File: blurfilter.h

  Declaration of pixel structure and blurfilter function.
    
 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_

/* NOTE: This structure must not be padded! */
typedef struct _pixel {
    unsigned char r,g,b;
} pixel;

void blurfilter_x(const int ystart, const int yend, const int xsize, const int ysize, pixel* src, pixel* dst, const int radius, const double *w);

void blurfilter_y(const int ystart, const int yend, const int xsize, const int ysize, pixel* src, pixel* dst, const int radius, const double *w);

#endif

