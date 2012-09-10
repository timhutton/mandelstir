/*  Copyright 2012 Tim Hutton

    This file is part of MandelStir.

    MandelStir is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MandelStir is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MandelStir. If not, see <http://www.gnu.org/licenses/>.         */

// http://code.google.com/p/mandelstir
// Contact: tim.hutton@gmail.com
// http://www.sq3.org.uk

// stdlib:
#include <math.h>

// STL:
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
using namespace std;

// ---------------------------------------------------------------------------------

// return the counter-clockwise angle [0,2pi]
double atanCCW(double y,double x)
{
    double a = atan2(y,x);
    if(a<0) a += 2.0 * 3.14159265358979;
    return a;
}

// ---------------------------------------------------------------------------------

// raise a complex number x+yi to the p'th power (where p is a real number) in-place
void ComplexPower(double &x,double &y,double p)
{
    // http://mathworld.wolfram.com/ComplexNumber.html
    double r2 = x*x + y*y;
    double rp = pow(r2,p/2.0);
    double theta = atanCCW(y,x);
    x = rp * cos(p*theta);
    y = rp * sin(p*theta);
}

// ---------------------------------------------------------------------------------

// draw a number on the image
void DrawNumber(unsigned char* image,int W,int H,double f);

// ---------------------------------------------------------------------------------

int main()
{
    const double aspect = 4.0/3.0;
    const double scale = 1.5;
    const int WIDTH = int(512*aspect*scale);
    const int HEIGHT = int(512*scale);
    unsigned char * image = new unsigned char[HEIGHT*WIDTH*4]; // row,column BGRA
    float * float_image = new float[HEIGHT*WIDTH*3]; // row,column RGB
    for(int i=0;i<HEIGHT*WIDTH;i++)
        image[i*4+3] = 255; // don't want transparency
    
    double xmin,xmax,ymin,ymax; // area to explore

    int max_its = 100;

    enum { Mandelbrot, Julia }                      algorithm =     Mandelbrot;
    enum { All, CircleRadius2, Checkerboard, Set }  source_points = Set;
    enum { Stirring, Linear }                       interpolation = Linear;

    const double sub = 0.1; // subsample the pixels for better resolution

    double speed;

    int max_frames;

    double zx,zy,cx,cy;
    double tempx,t;
    int r,g,b;

    if(algorithm==Mandelbrot)
    {
       xmin=-2;
       xmax=2;
       ymin=-1.5;
       ymax=1.5;
       speed = 0.005;
       max_frames = 4000; // or as high as you like
    }
    else
    {
        cx=-0.11; cy=0.6557;   // spirally blob
        //cx=-0.8;  cy=0.15;   // whiskery dragon
        //cx = -0.743643887037151l; cy = 0.131825904205330; // separated whorls
        //cx = 0.0; cy = -0.636; // romanesco broccoli
        xmin=-1.8;
        xmax=1.8;
        ymin=-1.5;
        ymax=1.5;
        speed = 0.005;
        max_frames = int(3.0 / speed); // just run a few times
    }

    for(int frame=0;frame<max_frames;frame++)
    {
        cout << "Frame " << frame+1 << " of " << max_frames << "\n";

        double target = frame * speed;

        // plot the static background image
        for(int y=0;y<HEIGHT;y++)
        {
            for(int x=0;x<WIDTH;x++)
            {
                zx = xmin + (xmax-xmin)*(x/double(WIDTH));
                zy = ymin + (ymax-ymin)*(y/double(HEIGHT));
                if(algorithm == Mandelbrot)
                {
                    cx = zx;
                    cy = zy;
                }
                r = 0; g = 0; b = 128; // dark blue inside
                for(int it=0;it<max_its;it++)
                {
                    tempx = zx*zx - zy*zy + cx;
                    zy = 2.0*zx*zy + cy;
                    zx = tempx;
                    if(zx*zx+zy*zy>4.0)
                    {
                        r = g = b = 0; // black outside
                        break;
                    }
                }
                image[y*WIDTH*4+x*4+0] = b; 
                image[y*WIDTH*4+x*4+1] = g;
                image[y*WIDTH*4+x*4+2] = r;
                float_image[y*WIDTH*3+x*3+0] = float_image[y*WIDTH*3+x*3+1] = float_image[y*WIDTH*3+x*3+2] = 0.0f;
            }
        }
        // show the iteration count (a floating point number!)
        DrawNumber(image,WIDTH,HEIGHT,target);
        int itarget = int(floor(target));
        t = target - itarget;
        // plot the moving points
        for( double subx = xmin; subx < xmax; subx += sub*(xmax-xmin)/double(WIDTH) )
        {
            for(double suby = ymin; suby < ymax; suby += sub*(ymax-ymin)/double(HEIGHT) )
            {
                switch(source_points)
                {
                    case All: default:
                        break;
                    case CircleRadius2:
                        if(subx*subx+suby*suby>4.0) continue;
                        break;
                    case Checkerboard:
                        {
                            if( ((int(WIDTH*(subx-xmin)/(xmax-xmin))/30)%2) ^ ((int(HEIGHT*(suby-ymin)/(ymax-ymin))/30)%2) )
                                continue;
                        }
                        break;
                    case Set:
                        {
                            // does this pixel escape?
                            zx = subx;
                            zy = suby;
                            if(algorithm == Mandelbrot)
                            {
                                cx = zx;
                                cy = zy;
                            }
                            bool escaped = false;
                            for(int it=0;it<max_its;it++)
                            {
                                tempx = zx*zx - zy*zy + cx;
                                zy = 2.0*zx*zy + cy;
                                zx = tempx;
                                if(zx*zx+zy*zy>4.0)
                                {
                                    escaped = true;
                                    break;
                                }
                            }
                            if(escaped) continue;
                        }
                        break;
                }
                // find the transformed location
                zx = subx;
                zy = suby;
                if(algorithm == Mandelbrot)
                {
                    cx = zx;
                    cy = zy;
                }
                switch(interpolation)
                {
                    case Linear:
                    {
                        double ox=zx,oy=zy; // remember the previous iteration
                        for(int i=0;i<itarget+1;i++)
                        {
                            ox = zx;
                            oy = zy;
                            tempx = zx*zx - zy*zy + cx;
                            zy = 2.0*zx*zy + cy;
                            zx = tempx;
                        }
                        // use linear interpolation on the last jump
                        zx = ox + (zx-ox)*t;
                        zy = oy + (zy-oy)*t;
                    }
                    break;
                    case Stirring:
                    {
                        // perform the integer stages
                        for(int i=0;i<itarget;i++)
                        {
                            tempx = zx*zx - zy*zy + cx;
                            zy = 2.0*zx*zy + cy;
                            zx = tempx;
                        }
                        // perform the leftover fractional iterations
                        ComplexPower(zx,zy,1.0+t);
                        zx += cx * t;
                        zy += cy * t;
                    }
                    break;
                }
                // brighten the new location (if on the image)
                if(zx>xmin && zx<xmax && zy>ymin && zy<ymax)
                {
                    int ix = int(  WIDTH * (zx-xmin) / (xmax-xmin) );
                    int iy = int( HEIGHT * (zy-ymin) / (ymax-ymin) );
                    float_image[iy*WIDTH*3+ix*3+0] += 1.0f;
                    float_image[iy*WIDTH*3+ix*3+1] += 1.0f;
                    float_image[iy*WIDTH*3+ix*3+2] += 1.0f;
                }
            }
        }

        // composite the float image and the background one
        for(int y=0;y<HEIGHT;y++)
        {
            for(int x=0;x<WIDTH;x++)
            {
                float f = float_image[y*WIDTH*3+x*3+0];
                if(f>0.0f)
                {
                    f = log(log(f))*120;
                    if(f<0) continue;
                    int ib = int(image[y*WIDTH*4+x*4+0] + f);
                    image[y*WIDTH*4+x*4+0] = min(255,max(0,ib));
                    int ig = int(image[y*WIDTH*4+x*4+1] + f);
                    image[y*WIDTH*4+x*4+1] = min(255,max(0,ig));
                    int ir = int(image[y*WIDTH*4+x*4+2] + f);
                    image[y*WIDTH*4+x*4+2] = min(255,max(0,ir));
                }
            }
        }

        // simple TGA output
        {
            ostringstream filename;
            filename << "linear/frame_" << frame << ".tga";
            ofstream o(filename.str().c_str(), ios::out | ios::binary);
            for(int c=0;c<12;c++) o.put("002000000000"[c]-'0');
            o.put((WIDTH & 0x00FF));
            o.put((WIDTH & 0xFF00) / 256);
            o.put((HEIGHT & 0x00FF));
            o.put((HEIGHT & 0xFF00) / 256);
            o.put(32);
            o.put(0);
            o.write((const char*)image,HEIGHT*WIDTH*4);
        }
    }

    delete []image;
    
    return EXIT_SUCCESS;
}

// ---------------------------------------------------------------------------------

void DrawChar(unsigned char* image,int W,int H,int &sx,char c)
{
    const char *segments[11] = {"012456","25","02346","02356","1235","01356","013456","025","0123456","012356","7"}; // 0123456789.
    const int coords[8][4] = {{1,0,3,1},{0,1,1,4},{3,1,4,4},{1,4,3,5},{0,5,1,8},{3,5,4,8},{1,8,3,9},{1,4,2,5}}; // {xmin,ymin,xmax,ymax}
    const int scale = 3;
    const int sy = scale*9*2;

    const char *segs;
    if(c>='0' && c<='9') segs = segments[c-'0'];
    else segs = segments[10];
    for(int iSeg=0;iSeg<strlen(segs);iSeg++)
    {
        int seg = segs[iSeg]-'0';
        for(int x=sx+coords[seg][0]*scale;x<=sx+coords[seg][2]*scale;x++)
            for(int y=sy-coords[seg][3]*scale;y<=sy-coords[seg][1]*scale;y++)
                image[y*W*4+x*4+0] = image[y*W*4+x*4+1] = image[y*W*4+x*4+2] = 255;
    }
    sx += scale*6;
}

// ---------------------------------------------------------------------------------

void DrawNumber(unsigned char* image,int W,int H,double f)
{
    ostringstream oss;
    oss << f;
    int sx=0;
    for(int iC=0;iC<oss.str().length();iC++)
        DrawChar(image,W,H,sx,oss.str()[iC]);
}

// ---------------------------------------------------------------------------------
