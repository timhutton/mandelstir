/*  Copyright 2020 Tim Hutton

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

// stdlib
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

struct Point {
    double x, y;
};

struct MandelbrotPoint {
    Point c;
    Point z;
    Point prev_z;

    void update()
    {
        prev_z.x = z.x;
        prev_z.y = z.y;
        z.x = prev_z.x * prev_z.x - prev_z.y * prev_z.y + c.x;
        z.y = 2 * prev_z.x * prev_z.y + c.y;
    }

    Point getIntermediate(double u) const
    {
        Point p;
        p.x = prev_z.x + u * (z.x - prev_z.x);
        p.y = prev_z.y + u * (z.y - prev_z.y);
        return p;
    }
};

struct Rect {
    double xmin, xmax, ymin, ymax;
};

void writeTGA(const vector<unsigned char>& char_bgra_image, int width, int height, const string& filename)
{
    ofstream o(filename, ios::out | ios::binary);
    for(int c = 0; c < 12; c++) o.put("002000000000"[c] - '0');
    o.put((width & 0x00FF));
    o.put((width & 0xFF00) / 256);
    o.put((height & 0x00FF));
    o.put((height & 0xFF00) / 256);
    o.put(32);
    o.put(0);
    o.write(reinterpret_cast<const char*>(char_bgra_image.data()), height * width * 4);
}

void accumulatePoint(double x, double y, vector<float>& float_rgb_image, const Rect& image_range, int width, int height)
{
    const int ix = static_cast<int>( width * (x - image_range.xmin) / (image_range.xmax - image_range.xmin));
    const int iy = static_cast<int>(height * (y - image_range.ymin) / (image_range.ymax - image_range.ymin));
    if(ix >= 0 && ix < width && iy >= 0 && iy < height)
    {
        const size_t iPixel = iy * width + ix;
        float_rgb_image[iPixel * 3 + 0] += 1;
        float_rgb_image[iPixel * 3 + 1] += 1;
        float_rgb_image[iPixel * 3 + 2] += 1;
    }
}

void writePointsToFloatImage(const vector<MandelbrotPoint>& points, double u, vector<float>& float_rgb_image, const Rect& image_range, int width, int height)
{
    fill(float_rgb_image.begin(), float_rgb_image.end(), 0.0f);
    for(const MandelbrotPoint& pt : points)
    {
        const Point p = pt.getIntermediate(u);
        if(p.x * p.x + p.y * p.y < 4)
        {
            // accumulate in both halves since we only track points starting with y > 0
            accumulatePoint(p.x,  p.y, float_rgb_image, image_range, width, height);
            accumulatePoint(p.x, -p.y, float_rgb_image, image_range, width, height);
        }
    }
}

unsigned char scaleFloatToChar(float val)
{
    if(val > 0)
    {
        val = log(log(val)) * 110;
    }
    return static_cast<unsigned char>(min(255.0f, max(0.0f, val)));
}

void writeFloatImageToCharImage(const vector<float>& float_rgb_image, vector<unsigned char>& char_bgra_image, int width, int height)
{
    for(int y = 0; y < height; y++)
    {
        for(int x = 0; x < width; x++)
        {
            const size_t iPixel = y * width + x;
            const float r = float_rgb_image[iPixel * 3 + 0];
            const float g = float_rgb_image[iPixel * 3 + 1];
            const float b = float_rgb_image[iPixel * 3 + 2];
            char_bgra_image[iPixel * 4 + 0] = scaleFloatToChar(b);
            char_bgra_image[iPixel * 4 + 1] = scaleFloatToChar(g);
            char_bgra_image[iPixel * 4 + 2] = scaleFloatToChar(r);
            char_bgra_image[iPixel * 4 + 3] = 255;
        }
    }
}

void initializeTrackingPoints(vector<MandelbrotPoint>& points, const Rect& sample_range, size_t samples_per_unit_length)
{
    const double d = 1 / static_cast<double>(samples_per_unit_length);
    for(double x = sample_range.xmin; x < sample_range.xmax; x += d)
    {
        for(double y = sample_range.ymin; y < sample_range.ymax; y += d)
        {
            if(x * x + y * y > 4)
            {
                continue; // skip points outside the r=2 disk
            }
            MandelbrotPoint pt;
            pt.c.x = pt.z.x = x;
            pt.c.y = pt.z.y = y;
            pt.prev_z.x = pt.prev_z.y = 0;
            if(points.size() < points.capacity())
            {
                points.push_back(pt);
            }
        }
    }
}

void updatePoints(vector<MandelbrotPoint>& points)
{
    for(MandelbrotPoint& pt : points)
    {
        pt.update();
    }
}

string getFilename(size_t iFrame)
{
    ostringstream oss;
    oss << "frame_" << setw(6) << setfill('0') << iFrame << ".tga";
    return oss.str();
}

void drawCharacterOnCharImage(vector<unsigned char>& image, int W, int H, int &sx, char c)
{
    const char *segments[11] = {"012456","25","02346","02356","1235","01356","013456","025","0123456","012356","7"}; // 0123456789.
    const int coords[8][4] = {{1,0,3,1},{0,1,1,4},{3,1,4,4},{1,4,3,5},{0,5,1,8},{3,5,4,8},{1,8,3,9},{1,4,2,5}}; // {xmin,ymin,xmax,ymax}
    const int scale = 3;
    const int sy = H - scale*9*2;

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

void drawNumberOnCharImage(vector<unsigned char>& image, int width, int height, double f)
{
    ostringstream oss;
    oss << f;
    int sx = width / 20;
    for(int iC = 0; iC < oss.str().length(); iC++)
    {
        drawCharacterOnCharImage(image, width, height, sx, oss.str()[iC]);
    }
}

void run()
{
    // define the region that we sample with tracking points (probably don't want to change this)
    constexpr Rect sample_range = { -2, 2, 0, 2 }; // (we duplicate the upper half to save compute)

    // sample the region at as high a resolution as we can afford
    constexpr size_t samples_per_unit_length = 8000; // samples per unit length on the complex plane
    constexpr size_t MAX_POINTS = static_cast<size_t>(samples_per_unit_length * samples_per_unit_length
                                            * (sample_range.ymax - sample_range.ymin) * (sample_range.xmax - sample_range.xmin));
    vector<MandelbrotPoint> points;
    points.reserve(MAX_POINTS);

    // define the region that we want to see in the final image (can change this freely)
    constexpr Rect image_range = { -1.8, 0.9, -0.8, 1.2};

    // make a destination image to write into
    constexpr int HEIGHT = 1024;
    constexpr int WIDTH = static_cast<int>(HEIGHT * (image_range.xmax - image_range.xmin) / (image_range.ymax - image_range.ymin));
    cout << "Writing into " << WIDTH << " x " << HEIGHT << " image\n";
    vector<float> float_rgb_image(HEIGHT*WIDTH*3); // row,column RGB - we accumulate point counts into this one
    vector<unsigned char> char_bgra_image(HEIGHT*WIDTH*4); // row,column BGRA - this contains our final image

    initializeTrackingPoints(points, sample_range, samples_per_unit_length);

    const size_t gigabytes = MAX_POINTS * sizeof(MandelbrotPoint) / (1024 * 1024 * 1024);
    const double million_pts = points.size() / 1e6;
    cout << "Tracking " << million_pts << " million points. Memory: " << gigabytes << "GB" << std::endl;

    const int n_iterations = 100;
    const int n_sub_iterations = 20;
    size_t iFrame = 0;
    for(size_t iIteration = 0; iIteration <= n_iterations; iIteration++)
    {
        updatePoints(points);
        for(size_t iSubIteration = 0; iSubIteration < n_sub_iterations; iSubIteration++)
        {
            const double u = iSubIteration / static_cast<double>(n_sub_iterations);
            writePointsToFloatImage(points, u, float_rgb_image, image_range, WIDTH, HEIGHT);
            writeFloatImageToCharImage(float_rgb_image, char_bgra_image, WIDTH, HEIGHT);
            drawNumberOnCharImage(char_bgra_image, WIDTH, HEIGHT, iIteration + u);
            const string filename = getFilename(iFrame++);
            writeTGA(char_bgra_image, WIDTH, HEIGHT, filename);
            cout << "Iteration " << iIteration + u << " - wrote " << filename << "\n";
        }
    }
}

int main()
{
    try
    {
        run();
    }
    catch(const exception& e)
    {
        cout << "Caught an exception: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    catch(...)
    {
        cout << "Caught an unknown exception\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
