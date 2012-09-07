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
    const int WIDTH = 614;
    const int HEIGHT = 512;
    unsigned char * image = new unsigned char[HEIGHT*WIDTH*4]; // row,column BGRA
    for(int i=0;i<HEIGHT*WIDTH;i++)
        image[i*4+3] = 255; // don't want transparency
    
    double xmin,xmax,ymin,ymax; // area to explore

    const int max_its = 300;

    const bool MSet = true;
    const bool move_set = true;

    double speed;

    int max_frames;

    double zx,zy,cx,cy;
    double tempx,t;
    int r,g,b;

    if(MSet)
    {
       xmin=-2;
       xmax=1;
       ymin=-1.5;
       ymax=1;
       speed = 0.01;
       max_frames = 4000; // or as high as you like
    }
    else
    {
        cx=-0.11; cy=0.6557;   // spirally blob
        //cx=-0.8;  cy=0.15;   // whiskery dragon
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

        // plot the static background image
        for(int y=0;y<HEIGHT;y++)
        {
            for(int x=0;x<WIDTH;x++)
            {
                zx = xmin + (xmax-xmin)*(x/double(WIDTH));
                zy = ymin + (ymax-ymin)*(y/double(HEIGHT));
                if(MSet) // (else J-Set)
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
            }
        }
        // show the iteration count (a floating point number!)
        double target = frame * speed;
        DrawNumber(image,WIDTH,HEIGHT,target);
        int itarget = int(floor(target));
        // plot the moving points
        const double sub = 0.25; // use a higher density of points
        for( double subx = xmin; subx < xmax; subx += sub*(xmax-xmin)/double(WIDTH) )
        {
            for(double suby = ymin; suby < ymax; suby += sub*(ymax-ymin)/double(HEIGHT) )
            {
                if(!move_set)
                {
                    // only move a checkerboard
                    int cb = ((int(WIDTH*(subx-xmin)/(xmax-xmin))/10)%2) ^ ((int(HEIGHT*(suby-ymin)/(ymax-ymin))/10)%2);
                    if(cb) continue;
                }
                else
                {
                    // does this pixel escape?
                    zx = subx;
                    zy = suby;
                    if(MSet) // (else J-Set)
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
                // find the transformed location
                zx = subx;
                zy = suby;
                if(MSet)
                {
                    cx = zx;
                    cy = zy;
                }
                // perform the integer stages
                for(int i=0;i<itarget;i++)
                {
                    tempx = zx*zx - zy*zy + cx;
                    zy = 2.0*zx*zy + cy;
                    zx = tempx;
                }
                // perform the leftover fractional iterations
                t = target - itarget;
                double pt = 1.0 + t;
                ComplexPower(zx,zy,pt);
                zx += cx * t;
                zy += cy * t;
                // brighten the new location (if on the image)
                if(zx>xmin && zx<xmax && zy>ymin && zy<ymax)
                {
                    int ix = int(WIDTH * ( (zx-xmin) / (xmax-xmin) ));
                    int iy = int(HEIGHT * ( (zy-ymin) / (ymax-ymin) ));
                    int ir = int(image[iy*WIDTH*4+ix*4+2]) + 9;
                    image[iy*WIDTH*4+ix*4+2] = min(255,max(0,ir));
                    int ig = int(image[iy*WIDTH*4+ix*4+1]) + 9;
                    image[iy*WIDTH*4+ix*4+1] = min(255,max(0,ig));
                    int ib = int(image[iy*WIDTH*4+ix*4+0]) + 9;
                    image[iy*WIDTH*4+ix*4+0] = min(255,max(0,ib));
                }
            }
        }

        // simple TGA output
        {
            ostringstream filename;
            filename << "frame_" << frame << ".tga";
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
    oss << setprecision(3) << f;
    int sx=0;
    for(int iC=0;iC<oss.str().length();iC++)
        DrawChar(image,W,H,sx,oss.str()[iC]);
}

// ---------------------------------------------------------------------------------
