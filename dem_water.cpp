/* ========================================================================
    File: @(#)dem_water.cpp
   ------------------------------------------------------------------------
    dem_water DEM waterbody flattening
    Copyright (C) 2013 Christoph Hormann <chris_hormann@gmx.de>
   ------------------------------------------------------------------------

    This file is part of dem_water

    glaciers_gen is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    dem_water is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with dem_water.  If not, see <http://www.gnu.org/licenses/>.

    Version history:

      0.1: initial public version, September 2013

   ========================================================================
 */

const char PROGRAM_TITLE[] = "DEM waterbody flattening 0.1";

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <string>

#include <gdal_priv.h>

#define cimg_use_magick 1
#define cimg_use_tiff 1
#define cimg_display 0

#include "CImg.h"

using namespace cimg_library;

const double Earth_Radius=6378.137;

const int xo[9] = { 0,-1, 0, 1,1,1,0,-1,-1 };
const int yo[9] = { 0,-1,-1,-1,0,1,1, 1, 0 };

struct Point {
  int x;
  int y;
  int z;
  Point(const int X = 0, const int Y = 0, const int Z = 0): x(X), y(Y), z(Z) {}
  bool operator<(const Point& pt) const {
    if( x < pt.x )
      return true;
    if( pt.x < x)
      return false;
    return y < pt.y;
  }

  bool operator==(const Point& pt) const {
    if (( x == pt.x ) && ( y == pt.y ))
      return true;
    return false;
  }

  bool operator!=(const Point& pt) const {
    if ( x != pt.x )
      return true;
    if ( y != pt.y )
      return true;
    return false;
  }
};

struct Waterbody {
	int min_x;
	int max_x;
	int min_y;
	int max_y;

	int pcnt;
	double el_sum;
	int pcnt_ol;
	double el_sum_ol;
	float el;
	float el_avg;
	float el_max;
	float el_min;
	float el_max_ol;
	float el_min_ol;
	float max_diff;
	Point max_diff_pnt;
};

struct Triangle {
	Point p1;
	Point p2;
	Point p3;
  Triangle(const Point &P1, const Point &P2, const Point &P3): p1(P1), p2(P2), p3(P3) {}
};


void load_image(std::string file, CImg<float> &img, int px, int py)
{
	GDALDataset  *poDataset;

	poDataset = (GDALDataset *) GDALOpen( file.c_str(), GA_ReadOnly );
	if( poDataset == NULL )
	{
		std::fprintf(stderr,"  opening file %s failed.\n\n", file.c_str());
		std::exit(1);
	}

	GDALRasterBand  *poBand;
	poBand = poDataset->GetRasterBand( 1 );
	float *pafScanline;
	int nXSize = poBand->GetXSize();
	int nYSize = poBand->GetYSize();

	CImg<short> img_raw = CImg<short>(nXSize, nYSize, 1, 1);

	poBand->RasterIO( GF_Read, 0, 0, nXSize, nYSize, 
										img_raw.data(), nXSize, nYSize, GDT_Int16, 
										0, 0 );

	int sz = 3600;
	if (nXSize < 3600) sz = 1200;

	for (int iy=0; iy <= sz; iy++)
	{
		for (int ix=0; ix <= sz; ix++)
		{
			img(px+ix,py+iy) = img_raw(ix,iy);
		}
	}
}

int main(int argc,char **argv)
{
	std::fprintf(stderr,"%s\n", PROGRAM_TITLE);
	std::fprintf(stderr,"-------------------------------------------------------\n");
	std::fprintf(stderr,"Copyright (C) 2013 Christoph Hormann\n");
	std::fprintf(stderr,"This program comes with ABSOLUTELY NO WARRANTY;\n");
	std::fprintf(stderr,"This is free software, and you are welcome to redistribute\n");
	std::fprintf(stderr,"it under certain conditions; see COPYING for details.\n");

  cimg_usage("Usage : dem_water [options]");

  // --- Read command line parameters ---

	// Files
  std::string file_m = cimg_option("-m","","mask image");
  std::string file_ms = cimg_option("-ms","","supersampled mask image");
  std::string file_b = cimg_option("-b","","bathymetry input image");
  std::string file_o = cimg_option("-o","","output image");
  std::string file_p = cimg_option("-p","","output POV-Ray file");
  std::string file_ce = cimg_option("-ce","","output coastline error file");
  std::string file_e = cimg_option("-e","","output error file");
	std::string file_i = cimg_option("-i","","DEM input image");

	std::string file_is1 = cimg_option("-is1","","SRTM input file 1");
	std::string file_is2 = cimg_option("-is2","","SRTM input file 2");
	std::string file_is3 = cimg_option("-is3","","SRTM input file 3");
	std::string file_is4 = cimg_option("-is4","","SRTM input file 4");
	std::string file_is5 = cimg_option("-is5","","SRTM input file 5");
	std::string file_is6 = cimg_option("-is6","","SRTM input file 6");
	std::string file_is7 = cimg_option("-is7","","SRTM input file 7");
	std::string file_is8 = cimg_option("-is8","","SRTM input file 8");
	std::string file_is9 = cimg_option("-is9","","SRTM input file 9");

	const int SRTMRes = cimg_option("-sr",1,"SRTM input resolution (1 or 3)");
	const float Rad = cimg_option("-r",2.1,"taper radius");
	const int Zero = cimg_option("-z",0,"height value for ocean areas");
	const int LandMin = cimg_option("-lm",-32767,"minimum height value for land areas");
	const int Wall = cimg_option("-w",0,"wall height around water areas");
	const int Carving = cimg_option("-c",0,"carving depth");
	const int Thr = cimg_option("-t",50,"elevation difference threshold");
	const int ErrorThr = cimg_option("-et",80,"elevation difference error threshold");
	const float Lon = cimg_option("-lon",0.0,"tile longitude");
	const float Lat = cimg_option("-lat",0.0,"tile latitude");
  const bool true3D = cimg_option("-3",false,"generate true 3d mesh");
  const bool out_cimg = cimg_option("-oc",false,"use cimg for output");
	const float Noise = cimg_option("-n",0.0,"noise to add to height values");

	const bool Debug = cimg_option("-debug",false,"generate debug ontput");

  const bool helpflag = cimg_option("-h",false,"Display this help");
  if (helpflag) std::exit(0);

	GDALAllRegister();

  // ---------------------------

	CImg<float> img_dem;

	int iWidth;
	int iHeight;

	int sWidth;
	int sHeight;

	int xStart;
	int yStart;

	int xEnd;
	int yEnd;

	if (file_i.size() != 0)
	{
		std::fprintf(stderr,"Loading DEM image...\n");

		img_dem = CImg<float>(file_i.c_str());
		iWidth = img_dem.width();
		iHeight = img_dem.height();
		sWidth = img_dem.width()-2;
		sHeight = img_dem.height()-2;
		xStart = 1;
		yStart = 1;
		xEnd = img_dem.width()-2;
		yEnd = img_dem.height()-2;
	}
	else
	{
		std::fprintf(stderr,"Loading DEM images...\n");

		iWidth = 3*3600/SRTMRes + 1;
		iHeight = 3*3600/SRTMRes + 1;
		sWidth = 3600/SRTMRes + 1;
		sHeight = 3600/SRTMRes + 1;
		xStart = 3600/SRTMRes;
		yStart = 3600/SRTMRes;
		xEnd = 2*3600/SRTMRes + 1;
		yEnd = 2*3600/SRTMRes + 1;

		img_dem = CImg<float>(iWidth,iHeight,1,1);

		cimg_forXY(img_dem,x,y)
		{
			img_dem(x,y) = -32767;
		}
		
		if (file_is1.size() != 0)
		{
			std::fprintf(stderr,"Loading SRTM input image 1...\n");
			load_image(file_is1, img_dem, 0, 0);
		}
		if (file_is2.size() != 0)
		{
			std::fprintf(stderr,"Loading SRTM input image 2...\n");
			load_image(file_is2, img_dem, 3600/SRTMRes, 0);
		}
		if (file_is3.size() != 0)
		{
			std::fprintf(stderr,"Loading SRTM input image 3...\n");
			load_image(file_is3, img_dem, 7200/SRTMRes, 0);
		}
		if (file_is4.size() != 0)
		{
			std::fprintf(stderr,"Loading SRTM input image 4...\n");
			load_image(file_is4, img_dem, 0, 3600/SRTMRes);
		}
		if (file_is6.size() != 0)
		{
			std::fprintf(stderr,"Loading SRTM input image 6...\n");
			load_image(file_is6, img_dem, 7200/SRTMRes, 3600/SRTMRes);
		}
		if (file_is7.size() != 0)
		{
			std::fprintf(stderr,"Loading SRTM input image 7...\n");
			load_image(file_is7, img_dem, 0, 7200/SRTMRes);
		}
		if (file_is8.size() != 0)
		{
			std::fprintf(stderr,"Loading SRTM input image 8...\n");
			load_image(file_is8, img_dem, 3600/SRTMRes, 7200/SRTMRes);
		}
		if (file_is9.size() != 0)
		{
			std::fprintf(stderr,"Loading SRTM input image 9...\n");
			load_image(file_is9, img_dem, 7200/SRTMRes, 7200/SRTMRes);
		}
		
		if (file_is5.size() != 0)
		{
			std::fprintf(stderr,"Loading SRTM input image 5...\n");
			load_image(file_is5, img_dem, 3600/SRTMRes, 3600/SRTMRes);
		}
	}

	CImg<short> img_bathy;

	if (file_b.size() != 0)
	{
		std::fprintf(stderr,"Loading bathymetry input image...\n");

		img_bathy = CImg<short>(file_b.c_str());
	}

	//img_dem.save("debug.tif");

	std::fprintf(stderr,"Loading mask image...\n");

	GDALDataset  *poDataset;

	poDataset = (GDALDataset *) GDALOpen( file_m.c_str(), GA_ReadOnly );
	if( poDataset == NULL )
	{
		std::fprintf(stderr,"  opening file %s failed.\n\n", file_m.c_str());
		return 1;
	}

	double adfGeoTransform[6], adfInvGeoTransform[6];

	if( poDataset->GetGeoTransform( adfGeoTransform ) != CE_None )
	{
		std::fprintf(stderr,"  error reading coordinates from file %s.\n\n", file_m.c_str());
		return 1;
	}

	GDALInvGeoTransform( adfGeoTransform, adfInvGeoTransform );

	GDALRasterBand  *poBand;
	poBand = poDataset->GetRasterBand( 1 );
	int nXSize = poBand->GetXSize();
	int nYSize = poBand->GetYSize();

	CImg<int> img_mask = CImg<int>(nXSize, nYSize, 1, 1);
	CImg<short> img_f = CImg<short>(nXSize, nYSize, 1, 1);

	poBand->RasterIO( GF_Read, 0, 0, nXSize, nYSize, 
										img_mask.data(), nXSize, nYSize, GDT_Int32, 
										0, 0 );

	CImg<unsigned char> img_ms;

	if (file_ms.size() != 0)
	{
		std::fprintf(stderr,"Loading supersampled mask image...\n");

		img_ms = CImg<unsigned char>(file_ms.c_str());
	}

	std::fprintf(stderr," dem image: %d x %d pixel\n", img_dem.width(), img_dem.height());
	std::fprintf(stderr," mask image: %d x %d pixel\n", img_mask.width(), img_mask.height());

	std::fprintf(stderr,"Analyzing waterbodies...\n");

	std::map<int,Waterbody> wbs;

	CImg<short> img_oe = CImg<short>(sWidth, sHeight, 1, 1);

	cimg_forXY(img_mask,px,py)
	{
		img_f(px,py) = 0;
		if ((px > 0) && (py > 0) && (px < img_mask.width()-1) && (py < img_mask.height()-1))
		if (img_mask(px,py) != 0)
		{
			std::map<int,Waterbody>::iterator it = wbs.find(img_mask(px,py));
			if (it == wbs.end())
			{
				Waterbody wb;
				wb.min_x = px;
				wb.max_x = px;
				wb.min_y = py;
				wb.max_y = py;
				wb.max_diff = -32767;
				if (img_mask(px,py) < 0) // Ocean
				{
					wb.el = 0;
					wb.el_avg = 0;
					wb.el_min = 0;
					wb.el_max = 0;
					wb.el_min_ol = 0;
					wb.el_max_ol = 0;
				}
				else
				{
					wb.pcnt = 1;
					wb.el_sum = img_dem(px,py);
					wb.el_min = img_dem(px,py);
					wb.el_max = img_dem(px,py);
					if ((img_mask(px+1,py) != img_mask(px,py))  || (img_mask(px-1,py) != img_mask(px,py)) || 
							(img_mask(px,py+1) != img_mask(px,py))  || (img_mask(px,py-1) != img_mask(px,py)))
					{
						wb.pcnt_ol = 1;
						wb.el_sum_ol = img_dem(px,py);
						wb.el_min_ol = img_dem(px,py);
						wb.el_max_ol = img_dem(px,py);
					}
				}
				wbs.insert(std::pair<int,Waterbody>(img_mask(px,py), wb));
			}
			else
			{
				(*it).second.min_x = std::min((*it).second.min_x, px);
				(*it).second.max_x = std::max((*it).second.max_x, px);
				(*it).second.min_y = std::min((*it).second.min_y, py);
				(*it).second.max_y = std::max((*it).second.max_y, py);
				if (img_mask(px,py) > 0) 
				{
					(*it).second.pcnt += 1;
					(*it).second.el_sum += img_dem(px,py);
					(*it).second.el_min = std::min((*it).second.el_min, img_dem(px,py));
					(*it).second.el_max = std::max((*it).second.el_max, img_dem(px,py));
					if ((img_mask(px+1,py) != img_mask(px,py))  || (img_mask(px-1,py) != img_mask(px,py)) || 
							(img_mask(px,py+1) != img_mask(px,py))  || (img_mask(px,py-1) != img_mask(px,py)))
					{
						(*it).second.pcnt_ol += 1;
						(*it).second.el_sum_ol += img_dem(px,py);
						(*it).second.el_min_ol = std::min((*it).second.el_min_ol, img_dem(px,py));
						(*it).second.el_max_ol = std::max((*it).second.el_max_ol, img_dem(px,py));
					}
				}
			}
		}
	}

	std::fprintf(stderr,"  found %d waterbodies.\n", wbs.size());

	CImg<float> img_dem_orig = img_dem;

	{
		std::fprintf(stderr,"Flattening waterbodies...\n");

		int cnt = 0;
		int cnto = 0;

		cimg_forXY(img_mask,px,py)
		{
			if (img_mask(px,py) < 0)
			{
				img_dem(px,py) = Zero;
				cnto++;
			}
			else if (img_mask(px,py) != 0)
			{
				std::map<int,Waterbody>::iterator it = wbs.find(img_mask(px,py));
				if (it != wbs.end())
				{
					float el = std::min((*it).second.el_min_ol, float((*it).second.el_sum/(*it).second.pcnt));
					if (el < ((*it).second.el_sum_ol/(*it).second.pcnt_ol-Thr))
						// outline average too low
						el = (*it).second.el_sum_ol/(*it).second.pcnt_ol - Thr*0.5;
					else
					{
						// heurestric to determine approximate lake level
						float d = ((*it).second.el_sum/(*it).second.pcnt) - el;
						if (d > 0.0) // outline minimum lower than area average
						{
							el += std::min(0.5*Thr, 0.75*d);
						}
					}
					(*it).second.el = el;
					img_f(px,py) = img_dem(px,py) - el;
					float f = 1.0;
					if (file_ms.size() != 0)
						if ((px >= xStart) && (py >= yStart) && (px <= xEnd) && (py <= yEnd))
							f = float(img_ms(px-xStart,py-yStart))/255.0;
					img_dem(px,py) += f*((el-Carving) - img_dem(px,py));
					cnt++;
				}
			}
		}

		if (Debug)
			img_dem.save("debug-f.tif");

		std::fprintf(stderr,"  flattened %d+%d pixel.\n", cnt, cnto);

		std::fprintf(stderr,"Computing differences...\n");

		bool Coastline_Error_Found = false;

		for (int iy=0; iy < yEnd-yStart; iy++)
		{
			for (int ix=0; ix < xEnd-xStart; ix++)
			{
				img_oe(ix,iy) = 0;

				if (img_mask(xStart+ix,yStart+iy) > 0)
				{
					std::map<int,Waterbody>::iterator it = wbs.find(img_mask(xStart+ix,yStart+iy));
					if (it != wbs.end())
					{
						if (std::abs(img_dem_orig(xStart+ix,yStart+iy)-img_dem(xStart+ix,yStart+iy)) > (*it).second.max_diff)
						{
							(*it).second.max_diff = std::abs(img_dem_orig(xStart+ix,yStart+iy)-img_dem(xStart+ix,yStart+iy));
							(*it).second.max_diff_pnt = Point(ix,iy);
						}
					}
				}
				else if (img_mask(xStart+ix,yStart+iy) < 0)
				{
					if (img_dem_orig(xStart+ix,yStart+iy) > Zero+ErrorThr)
					{
						img_oe(ix,iy) = img_dem_orig(xStart+ix,yStart+iy)-Zero;
						Coastline_Error_Found = true;
					}
				}
			}
		}

		// thin out dense error points leaving only the largest
		for (int iy=0; iy < yEnd-yStart; iy++)
		{
			for (int ix=0; ix < xEnd-xStart; ix++)
			{
				if (img_oe(ix,iy) > 0)
				{
					int d = 12;
					bool Found = false;
					for (int yn=iy-d; yn <=iy+d; yn++)
						for (int xn=ix-d; xn <=ix+d; xn++)
							if (!Found)
								if (xn >= 0)
									if (yn >= 0)
										if (xn < xEnd-xStart)
											if (yn < xEnd-xStart)
												if (img_oe(xn,yn) > img_oe(ix,iy))
												{
													img_oe(ix,iy) = 0;
													Found = true;
													break;
												}
				}
			}
		}
		
		if (file_ce.size() != 0)
		{
			if (Coastline_Error_Found)
			{
				std::fprintf(stderr,"  writing coastline error list...\n");

				std::ofstream OStream ( file_ce.c_str() );
				OStream.precision(8);
				if ( OStream.is_open() )
				{
					OStream << "{\n";
					OStream << "\"type\": \"FeatureCollection\",\n";
					OStream << "\"features\": [\n";

					for (int iy=0; iy < yEnd-yStart; iy++)
					{
						for (int ix=0; ix < xEnd-xStart; ix++)
						{
							if (img_oe(ix,iy) > 0)
							{
								OStream << "{\n";
								OStream << "  \"type\":\"Feature\",\n";
								OStream << "  \"properties\":{ " ;
								OStream << "\"error\": " << img_oe(ix,iy);
								OStream << "}, \n" ;
								OStream << "  \"geometry\":{ ";
								OStream << "\"type\": \"Point\", \"coordinates\": [" << Lon + (float(ix)-0.5)/sWidth << ", " << Lat + 1.0 - (float(iy)-0.5)/sHeight << "] ";
								OStream << "} \n" ;
								OStream << "}, \n" ;
							}
						}
					}
					
					OStream << "\n]\n";
					OStream << "}\n";
				}
			}
		}
	}

	if ((file_e.size() != 0) && (wbs.size() != 0))
	{
		std::ofstream OStream ( file_e.c_str() );
		if ( OStream.is_open() )
		{
			std::fprintf(stderr,"Generating error list...\n");

			OStream << "{\n";
			OStream << "\"type\": \"FeatureCollection\",\n";
			OStream << "\"features\": [\n";

			for (std::map<int,Waterbody>::iterator it = wbs.begin(); it != wbs.end(); ++it)
			{
				if ((*it).first < 0) continue;
				if ((*it).second.max_diff >= ErrorThr)
				{
					OStream << "{\n";
					OStream << "  \"type\":\"Feature\",\n";
					OStream << "  \"properties\":{ " ;
					OStream << "\"id\": " << (*it).first;
					OStream << ", \"error\": " << (*it).second.max_diff;
					OStream << "}, \n" ;
					OStream << "  \"geometry\":{ ";
					OStream << "\"type\": \"Point\", \"coordinates\": [" << Lon + (float((*it).second.max_diff_pnt.x)-0.5)/sWidth << ", " << Lat + 1.0 - (float((*it).second.max_diff_pnt.y)-0.5)/sHeight << "] ";
					OStream << "} \n" ;
					OStream << "}, \n" ;
				}
			}

			OStream << "\n]\n";
			OStream << "}\n";
		}
	}

	cimg_forXY(img_mask,px,py)
	{
		if (Wall > 0)
		if (img_mask(px,py) == 0)
		{
			float el_max = -32767.0;
			for (int i = 1; i < 9; i++)
			{
				int xn = px + xo[i];
				int yn = py + yo[i];
				if (xn >= 0)
					if (yn >= 0)
						if (xn < img_mask.width())
							if (yn < img_mask.height())
								if (img_mask(xn,yn) > 0)
								{
									std::map<int,Waterbody>::iterator it = wbs.find(img_mask(xn,yn));
									if (it != wbs.end())
									{
										el_max = std::max(el_max, img_dem(xn,yn));
									}
								}
			}
			if (el_max > -32000)
			{
				float f = 1.0;
				if (file_ms.size() != 0)
					if ((px >= xStart) && (py >= yStart) && (px <= xEnd) && (py <= yEnd))
						f = 1.0-float(img_ms(px,py))/255.0;
				img_dem(px,py) += f*(std::max(el_max, img_dem(px,py) + Wall) - img_dem(px,py));
				//img_dem(px,py) = std::max(el_max, img_dem(px,py) + Wall);
			}
		}
	}

	if (Debug)
		img_dem.save("debug-w.tif");

	std::fprintf(stderr,"Tapering water areas...\n");
	{
		CImg<unsigned char> morph_mask(Rad*2 + 1, Rad*2 + 1, 1, 1);

		cimg_forXY(morph_mask,px,py)
		{
			int dx = px-(morph_mask.width()/2);
			int dy = py-(morph_mask.height()/2);
			if (std::sqrt(dx*dx+dy*dy) <= Rad)
				morph_mask(px,py) = 255;
			else
				morph_mask(px,py) = 0;
		}

		std::fprintf(stderr,"  step 1...\n");

		CImg<float> img_el(img_mask.width(), img_mask.height(), 1, 1);
		CImg<float> img_eld(img_mask.width(), img_mask.height(), 1, 1);
		CImg<unsigned char> img_m(img_mask.width(), img_mask.height(), 1, 1);

		cimg_forXY(img_mask,px,py)
		{
			if (img_mask(px,py) != 0)
			{
				img_m(px,py) = 255;
				img_el(px,py) = img_dem(px,py);
				img_eld(px,py) = img_dem_orig(px,py)-img_dem(px,py);
			}
			else
			{
				img_m(px,py) = 0;
				img_el(px,py) = -32767.0;
				img_eld(px,py) = 0.0;
			}
		}
		img_el.dilate(morph_mask);
		img_eld.dilate(morph_mask);

		std::fprintf(stderr,"  step 2...\n");

		CImg<float> img_dist = img_m.get_distance(255);

		if (file_ms.size() != 0)
		{
			cimg_forXY(img_mask,px,py)
			{
				if (img_ms(px,py) != 0)
					img_dist(px,py) = 1.0 - float(img_ms(px,py))/255.0;
			}
		}

		std::fprintf(stderr,"  step 3...\n");

		cimg_forXY(img_dem,px,py)
		{
			if ((img_dist(px,py) <= Rad) && (img_el(px,py) > -32000))
			if (img_m(px,py) != 255)
			{
				img_dem(px,py) = std::max(img_el(px,py)+Wall, float(img_dem(px,py)-(1.0 - img_dist(px,py)/Rad)*img_eld(px,py)));
			}
		}

		if (Debug)
			img_dem.save("debug-t.tif");
	}

	if (file_b.size() != 0)
	{
		std::fprintf(stderr,"Merging bathymetry...\n");

		CImg<unsigned char> img_m(img_mask.width(), img_mask.height(), 1, 1);

		cimg_forXY(img_mask,px,py)
		{
			if (img_mask(px,py) >= 0) img_m(px,py) = 0;
			else img_m(px,py) = 255;
		}

		CImg<float> img_dist = img_m.get_distance(0);

		int cntb = 0;

		for (int iy=0; iy < yEnd-yStart; iy++)
		{
			for (int ix=0; ix < xEnd-xStart; ix++)
			{
				if (file_b.size() != 0)
				{
					if (img_mask(xStart+ix,yStart+iy) < 0)
					{
						img_dem(xStart+ix,yStart+iy) = std::min(int(img_bathy(ix,iy)), Zero);
						cntb++;
					}
					else
						img_dem(xStart+ix,yStart+iy) = std::max(img_dem(xStart+ix,yStart+iy), float(LandMin));
				}
				else img_dem(xStart+ix,yStart+iy) = std::max(img_dem(xStart+ix,yStart+iy), float(LandMin));
			}
		}
		std::fprintf(stderr,"  %d bathymetry pixel.\n", cntb);
	}

	if ((file_p.size() != 0) && (wbs.size() != 0))
	{
		std::ofstream OStream ( file_p.c_str() );
		if ( OStream.is_open() )
		{
			std::fprintf(stderr,"Generating POV-Ray mesh...\n");

			OStream.precision(10);

			std::vector<Triangle> mesh;

			for (std::map<int,Waterbody>::iterator it = wbs.begin(); it != wbs.end(); ++it)
			{
				if ((*it).first < 0) continue;
				// generate mesh
				for (int py = (*it).second.min_y-1; py <= (*it).second.max_y+1; py++)
				for (int px = (*it).second.min_x-1; px <= (*it).second.max_x+1; px++)
					if ((px >= xStart) && (py >= yStart) && (px <= xEnd) && (py <= yEnd))
					{
						if ((img_mask(px,py) == (*it).first) || 
								(img_mask(px+1,py) == (*it).first) || (img_mask(px-1,py) == (*it).first) ||
								(img_mask(px,py+1) == (*it).first) || (img_mask(px,py-1) == (*it).first))
						{
							if ((img_mask(px+1,py+1) == (*it).first) ||
									(img_mask(px+2,py+1) == (*it).first) || (img_mask(px,py+1) == (*it).first) ||
									(img_mask(px+1,py+2) == (*it).first) || (img_mask(px+1,py) == (*it).first))
							{
								if ((img_mask(px,py+1) == (*it).first) ||
										(img_mask(px+1,py+1) == (*it).first) || (img_mask(px-1,py+1) == (*it).first) ||
										(img_mask(px,py+2) == (*it).first) || (img_mask(px,py) == (*it).first))
								{
									mesh.push_back(Triangle(Point(px,py,(*it).second.el), Point(px,py+1,(*it).second.el), Point(px+1,py+1,(*it).second.el)));
								}
								if ((img_mask(px+1,py) == (*it).first) ||
										(img_mask(px+2,py) == (*it).first) || (img_mask(px,py) == (*it).first) ||
										(img_mask(px+1,py+1) == (*it).first) || (img_mask(px+1,py-1) == (*it).first))
								{
									mesh.push_back(Triangle(Point(px,py,(*it).second.el), Point(px+1,py+1,(*it).second.el), Point(px+1,py,(*it).second.el)));
								}
							}
							if (img_mask(px+1,py+1) != (*it).first)
							{
								if ((img_mask(px,py+1) == (*it).first) ||
										(img_mask(px+1,py+1) == (*it).first) || (img_mask(px-1,py+1) == (*it).first) ||
										(img_mask(px,py+2) == (*it).first) || (img_mask(px,py) == (*it).first))
									if ((img_mask(px+1,py) == (*it).first) ||
											(img_mask(px+2,py) == (*it).first) || (img_mask(px,py) == (*it).first) ||
											(img_mask(px+1,py+1) == (*it).first) || (img_mask(px+1,py-1) == (*it).first))
									mesh.push_back(Triangle(Point(px,py,(*it).second.el), Point(px,py+1,(*it).second.el), Point(px+1,py,(*it).second.el)));
							}

							if (img_mask(px-1,py-1) != (*it).first)
							{
								if ((img_mask(px,py-1) == (*it).first) ||
										(img_mask(px+1,py-1) == (*it).first) || (img_mask(px-1,py-1) == (*it).first) ||
										(img_mask(px,py) == (*it).first) || (img_mask(px,py-2) == (*it).first))
									if ((img_mask(px-1,py) == (*it).first) ||
											(img_mask(px,py) == (*it).first) || (img_mask(px-2,py) == (*it).first) ||
											(img_mask(px-1,py+1) == (*it).first) || (img_mask(px-1,py-1) == (*it).first))
										mesh.push_back(Triangle(Point(px,py,(*it).second.el), Point(px,py-1,(*it).second.el), Point(px-1,py,(*it).second.el)));

							}
						}
						
					}
			}

			std::fprintf(stderr,"  %d triangles.\n", mesh.size());

			std::vector<Point> vertices;
			std::map<Point,int> vertex_ids;
			int vi = 0;

			// generate vertex list
			for (std::vector<Triangle>::iterator itm = mesh.begin(); itm != mesh.end(); ++itm)
			{
				std::map<Point,int>::iterator itv = vertex_ids.find(Point((*itm).p1.x, (*itm).p1.y));
				if (itv == vertex_ids.end())
				{
					vertex_ids.insert(std::pair<Point,int>(Point((*itm).p1.x, (*itm).p1.y), vi));
					vertices.push_back((*itm).p1);
					vi++;
				}
				itv = vertex_ids.find(Point((*itm).p2.x, (*itm).p2.y));
				if (itv == vertex_ids.end())
				{
					vertex_ids.insert(std::pair<Point,int>(Point((*itm).p2.x, (*itm).p2.y), vi));
					vertices.push_back((*itm).p2);
					vi++;
				}
				itv = vertex_ids.find(Point((*itm).p3.x, (*itm).p3.y));
				if (itv == vertex_ids.end())
				{
					vertex_ids.insert(std::pair<Point,int>(Point((*itm).p3.x, (*itm).p3.y), vi));
					vertices.push_back((*itm).p3);
					vi++;
				}
			}

			std::fprintf(stderr,"  %d vertices.\n", vertices.size());

			// write out mesh2
			std::fprintf(stderr,"  writing mesh2...\n");

			OStream << "// POV-Ray mesh2 object for waterbodies" << std::endl;
			OStream << "// written by " << PROGRAM_TITLE << std::endl;
			OStream << std::endl;

			if (vertices.size() > 0)
			{
				OStream << "mesh2 {" << std::endl;
				OStream << "vertex_vectors {" << std::endl;
				OStream << vertices.size();
				for (std::vector<Point>::iterator it = vertices.begin(); it != vertices.end(); ++it)
				{
					double x, y, z;
					double lon = adfGeoTransform[0] + adfGeoTransform[1] * (*it).x + adfGeoTransform[2] * (*it).y;
					double lat = adfGeoTransform[3] + adfGeoTransform[4] * (*it).x + adfGeoTransform[5] * (*it).y;
					if (true3D)
					{
						double lon_r = lon*cimg::PI/180.0;
						double lat_r = lat*cimg::PI/180.0;
						x = -(Earth_Radius+(*it).z*0.001)*cos(lat_r) * sin(lon_r);
						y = (Earth_Radius+(*it).z*0.001)*cos(lat_r) * cos(lon_r);
						z = (Earth_Radius+(*it).z*0.001)*sin(lat_r);
					}
					else
					{
						x = lon;
						y = lat;
						z = (*it).z;
					}

					OStream << "," << std::endl;
					OStream << "<" << x << "," << y << "," << z << ">";
				}

				OStream << std::endl;
				OStream << "}" << std::endl;
			
				if (true3D)
				{
					OStream.precision(5);
					// normals on sphere: same as coordinates
					OStream << "normal_vectors {" << std::endl;
					OStream << vertices.size();
					for (std::vector<Point>::iterator it = vertices.begin(); it != vertices.end(); ++it)
					{
						double x, y, z;
						double lon = adfGeoTransform[0] + adfGeoTransform[1] * (*it).x + adfGeoTransform[2] * (*it).y;
						double lat = adfGeoTransform[3] + adfGeoTransform[4] * (*it).x + adfGeoTransform[5] * (*it).y;
						double lon_r = lon*cimg::PI/180.0;
						double lat_r = lat*cimg::PI/180.0;
						x = -cos(lat_r) * sin(lon_r);
						y = cos(lat_r) * cos(lon_r);
						z = sin(lat_r);

						OStream << "," << std::endl;
						OStream << "<" << x << "," << y << "," << z << ">";
					}
					OStream << std::endl;
					OStream << "}" << std::endl;
				}
			
				OStream <<"  face_indices {" << std::endl;
				OStream << mesh.size();
				for (std::vector<Triangle>::iterator itm = mesh.begin(); itm != mesh.end(); ++itm)
				{
					OStream << "," << std::endl;
					OStream << "<" << vertex_ids[Point((*itm).p1.x, (*itm).p1.y)] << "," << vertex_ids[Point((*itm).p2.x, (*itm).p2.y)] << "," << vertex_ids[Point((*itm).p3.x, (*itm).p3.y)] << ">";
				}
				OStream << std::endl;
				OStream << "}" << std::endl;
			
				OStream << "}" << std::endl;
			}

			std::fprintf(stderr,"  mesh written to %s\n", file_p.c_str());
		}
	}

	if (file_o.size() != 0)
	{
		std::fprintf(stderr,"Writing results...\n");

		if (Noise > 0) img_dem.noise(Noise);

		if (out_cimg)
		{
			img_dem.crop(xStart, yStart, 0, 0, xEnd, yEnd, 0, 0);
			img_dem.save(file_o.c_str());			
		}
		else
		{
			std::ofstream OStream ( file_o.c_str(), std::ios_base::out | std::ios_base::binary );
			if ( OStream.is_open() )
			{
				for (int iy=0; iy < yEnd-yStart; iy++)
				{
					for (int ix=0; ix < xEnd-xStart; ix++)
					{
						union v {
							short	s;
							char c[2];
						};
	
						union v val;

						char buf;
						val.s = img_dem(xStart+ix,yStart+iy);
						buf = val.c[0];
						val.c[0] = val.c[1];
						val.c[1] = buf;

						OStream.write (val.c,2);
					}
				}
			}
		}

		std::fprintf(stderr,"Results written to file %s\n", file_o.c_str());
	}
}
