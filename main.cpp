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
    Point(double x, double y) : x(x), y(y) {}
};

bool inR2Disk(double x, double y)
{
    return x * x + y * y < 4;
}

bool inBulb1(double x, double y)
{
    const double cx = -1;
    const double cy = 0;
    const double r = 0.25;
    return (x - cx) * (x - cx) + (y - cy) * (y - cy) < r * r;
}

bool inBulb2(double x, double y)
{
    // this is only approximate
    const double cx = -0.1253;
    const double cy = 0.745;
    const double r = 0.094;
    return (x - cx) * (x - cx) + (y - cy) * (y - cy) < r * r;
}

bool inCardioid(double x, double y)
{
    const double tx = 0.25 - x;
    const double ty = y;
    const double theta = atan2(ty, tx);
    const double cr = pow(cos(theta / 2), 2);
    const double r2 = tx * tx + ty * ty;
    return r2 < cr * cr;
}

struct MandelbrotPoint {
    Point c;
    Point z;
    Point prev_z;
    bool is_in_set;

    MandelbrotPoint(double x, double y)
        : c(x, y)
        , z(x, y)
        , prev_z(0, 0)
    {
        is_in_set = isInSet();
    }

    void update()
    {
        prev_z.x = z.x;
        prev_z.y = z.y;
        z.x = prev_z.x * prev_z.x - prev_z.y * prev_z.y + c.x;
        z.y = 2 * prev_z.x * prev_z.y + c.y;
    }

    Point getIntermediate(double u) const
    {
        return Point(prev_z.x + u * (z.x - prev_z.x), prev_z.y + u * (z.y - prev_z.y));
    }

private:

    bool isInSet(size_t n_iterations = 100) const
    {
        double x = c.x;
        double y = c.y;
        for(size_t i = 0; i < n_iterations; i++)
        {
            double tempx = x * x - y * y + c.x;
            y = 2 * x * y + c.y;
            x = tempx;
            if(!inR2Disk(x,y))
            {
                return false;
            }
        }
        return true;
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

void accumulatePoint(double x, double y, bool is_in_set, vector<float>& float_rgb_image, const Rect& image_range, int width, int height)
{
    const int ix = static_cast<int>( width * (x - image_range.xmin) / (image_range.xmax - image_range.xmin));
    const int iy = static_cast<int>(height * (y - image_range.ymin) / (image_range.ymax - image_range.ymin));
    if(ix >= 0 && ix < width && iy >= 0 && iy < height)
    {
        const size_t iPixel = iy * width + ix;
        const size_t iFloat = iPixel * 3 + (is_in_set ? 1 : 0); // green: in set, red: out of set
        float_rgb_image[iFloat] += 1;
    }
}

void writePointsToFloatImage(const vector<MandelbrotPoint>& points,
                             double u,
                             vector<float>& float_rgb_image,
                             const Rect& image_range,
                             int width,
                             int height,
                             bool mirror)
{
    fill(float_rgb_image.begin(), float_rgb_image.end(), 0.0f);
    for(const MandelbrotPoint& pt : points)
    {
        const Point p = pt.getIntermediate(u);
        if(inR2Disk(p.x, p.y))
        {
            accumulatePoint(p.x,  p.y, pt.is_in_set, float_rgb_image, image_range, width, height);
            if(mirror)
            {
                // accumulate in both halves since we only track points starting with y > 0
                accumulatePoint(p.x, -p.y, pt.is_in_set, float_rgb_image, image_range, width, height);
            }
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
            char_bgra_image[iPixel * 4 + 1] = scaleFloatToChar(g);
            char_bgra_image[iPixel * 4 + 2] = scaleFloatToChar(r);
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
            if(inBulb2(x, y))
            {
                MandelbrotPoint pt(x, y);
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

void writeBackgroundImage(vector<unsigned char>& char_bgra_background_image, int width, int height, const Rect& image_range)
{
    for(int y = 0; y < height; y++)
    {
        for(int x = 0; x < width; x++)
        {
            MandelbrotPoint pt(image_range.xmin + x * (image_range.xmax - image_range.xmin) / width,
                               image_range.ymin + y * (image_range.ymax - image_range.ymin) / height);
            const size_t iPixel = y * width + x;
            char_bgra_background_image[iPixel * 4 + 0] = pt.is_in_set ? 150 : 0; // Mandelbrot set in blue pixel
            char_bgra_background_image[iPixel * 4 + 1] = 0;
            char_bgra_background_image[iPixel * 4 + 2] = 0;
            char_bgra_background_image[iPixel * 4 + 3] = 255; // alpha=255 for opaque
        }
    }
}

void copyBackgroundImageToCharImage(const vector<unsigned char>& char_bgra_background_image, vector<unsigned char>& char_bgra_image)
{
    char_bgra_image.assign(char_bgra_background_image.begin(), char_bgra_background_image.end());
}

void run()
{
    // define the region that we sample with tracking points (probably don't want to change this)
    constexpr Rect sample_range = { -2, 2, 0, 2 }; // (we duplicate the upper half to save compute)
    constexpr bool mirror = false;

    // sample the region at as high a resolution as we can afford
    constexpr size_t samples_per_unit_length = 8000; // samples per unit length on the complex plane
    constexpr size_t MAX_POINTS = static_cast<size_t>(samples_per_unit_length * samples_per_unit_length
                                            * (sample_range.ymax - sample_range.ymin) * (sample_range.xmax - sample_range.xmin));
    const size_t gigabytes = MAX_POINTS * sizeof(MandelbrotPoint) / (1024 * 1024 * 1024);
    vector<MandelbrotPoint> points;
    points.reserve(MAX_POINTS);

    // define the region that we want to see in the final image (can change this freely)
    constexpr Rect image_range = { -1.4, 0.7, -0.5, 1.1};

    // make a destination image to write into
    constexpr int HEIGHT = 1024;
    constexpr int WIDTH = static_cast<int>(HEIGHT * (image_range.xmax - image_range.xmin) / (image_range.ymax - image_range.ymin));
    vector<float> float_rgb_image(HEIGHT*WIDTH*3); // row,column RGB - we accumulate point counts into this one
    vector<unsigned char> char_bgra_background_image(HEIGHT*WIDTH*4); // row,column BGRA - this contains a static background image
    vector<unsigned char> char_bgra_image(HEIGHT*WIDTH*4); // row,column BGRA - this contains our final image

    cout << "Writing the background image..." << endl;
    writeBackgroundImage(char_bgra_background_image, WIDTH, HEIGHT, image_range);
    cout << "Initializing the tracking points..." << endl;
    initializeTrackingPoints(points, sample_range, samples_per_unit_length);

    const double million_pts = points.size() / 1e6;
    cout << "Tracking " << million_pts << " million points. (" << gigabytes << " GB)" << endl;
    cout << "Writing into " << WIDTH << " x " << HEIGHT << " image\n";

    const int n_iterations = 1000;
    int iFrame = 0;
    for(int iIteration = 0; iIteration <= n_iterations; iIteration++)
    {
        const int n_sub_iterations = max(1, 50 - iIteration);
        updatePoints(points);
        for(int iSubIteration = 0; iSubIteration < n_sub_iterations; iSubIteration++)
        {
            const double u = iSubIteration / static_cast<double>(n_sub_iterations);
            cout << "Iteration " << iIteration + u << " ";
            writePointsToFloatImage(points, u, float_rgb_image, image_range, WIDTH, HEIGHT, mirror);
            copyBackgroundImageToCharImage(char_bgra_background_image, char_bgra_image);
            writeFloatImageToCharImage(float_rgb_image, char_bgra_image, WIDTH, HEIGHT);
            drawNumberOnCharImage(char_bgra_image, WIDTH, HEIGHT, iIteration + u);
            const string filename = getFilename(iFrame++);
            writeTGA(char_bgra_image, WIDTH, HEIGHT, filename);
            cout << "Wrote " << filename << endl;
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
        cout << "Caught an exception: " << e.what() << endl;
        return EXIT_FAILURE;
    }
    catch(...)
    {
        cout << "Caught an unknown exception" << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
