#ifndef DATA_H_ABHYUI
#define DATA_H_ABHYUI

#include <QtOpenGL>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <map>

#include "myfont.h"
#include "MfixData.h"


template <typename T1 , typename T2>
T2 Convert(const T1 & t1)
{
  T2 t2;
  
  std::stringstream ss;
  ss << t1;
  ss >> t2;
  
  return t2;
}


struct Particle
{
  double x,y,z,d,u,v,w;
};

struct Position
{
    double x,y,z;
};

struct Data
{  
  
  static MfixData mfix;
  static MyFont myfont;
  static bool bBusy;
  static bool bFake_diameter;
  static bool bImage;
  
  static float xmin;
  static float xmax;
  static float ymin;
  static float ymax;
  
  static float scale;
  static float xtrans;
  static float ytrans;
  static int   xrot;
  static int   yrot;

  static int sphere_resolution;

  static double time;
  
  static bool bPaused;
  static bool b3D_bubbles;
  static bool bAdd_to_vpos;
  static bool bTrace;

  static double bAdd_time;
  static double vel_max;
  
  static int gl_width;
  static int gl_height;
  
  static int npixels_x;
  static int npixels_y;
    
  static int color_method;
  
  static int cur_ts;
  
  static std::string run_name;
  
  static std::vector<Particle> vP;
  static std::vector<Particle> vPos;
  static std::vector<double>   vTimes;
  static std::vector<int>      vColors_position;
  static std::vector<int>      vColors_diameter;

      
  static void init();  
  static void ReadNext_XML_file();
  static void GetTimes();

  static float mfix_xlength;
  static float mfix_ylength;
  static float mfix_zlength;
  static bool  mfix_bCartesian;
  static int   mfix_imax2;
  static int   mfix_jmax2;
  static int   mfix_kmax2;

};

#endif

