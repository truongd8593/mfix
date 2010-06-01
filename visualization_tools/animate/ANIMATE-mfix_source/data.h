#ifndef DATA_H_ABHYUI
#define DATA_H_ABHYUI

#include <QtOpenGL>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <map>

#include "MfixData.h"
#include "myfont.h"
//#include "some_font.h"

template <typename T1 , typename T2>
T2 Convert(const T1 & t1)
{
  T2 t2;
  
  std::stringstream ss;
  ss << t1;
  ss >> t2;
  
  return t2;
}




struct WindowVariable
{
  int                ts;
  int                mfix_variable_index;
  
  bool               bMustCalculateColors;
  bool               bColorsCalculated;
  
  bool               bUsed;
  bool               bMustRead;
};


struct MfixVariable
{
  std::vector<float> values;
  std::vector<int>   colors;
  bool               bUsed;
  bool               bMustCalculateColor;
};



struct PlotWindowInfo
{  
  bool bIJ_plot;
  bool bIK_plot; // not needed unless JK plot becones an option
  bool bContourFilled;
  bool bContourLines;
  bool bPlotGasVectors;
  bool bPlotSolidVectors;
  bool bPlotGrid;
  bool bFlipX;
  bool bPlot;
  
  int  win_var;

  int  win_i;
  int  win_j;
  int  win_k;
  
  int  win_imin;
  int  win_jmin;
  int  win_kmin;

  int  win_imax;
  int  win_jmax;
  int  win_kmax;
  
  
  int  win_m;
  int  win_ng;
  int  win_ns;
  int  win_scalar;
  int  win_rr;
  
  float contour_line_min;
  float contour_line_max;
  int   n_contours;
  
  PlotWindowInfo();
};

struct VariableRange
{
  float vmin , vmax;
};

struct MyTable
{
  float              time;
  std::vector<float> values;
};

struct Data
{
  static std::vector<MyTable> vTable;
  
  static bool bTable;
  static bool bSlice_only;
  static bool bBusy;
  static bool bScaleChanged;

  
  static MyFont myfont;
  static bool bPaused;
  
  static int gl_width;
  static int gl_height;
  
  static int npixels_x;
  static int npixels_y;
    
  static MfixData mfix;

  static std::string sCWD;
  
  // information for each plot window
  
  static std::vector<PlotWindowInfo> vPlotWindowInfo;
  static std::vector<VariableRange>  vVariableRange;
  static std::vector<WindowVariable> vWindowVariable;

  
  static int cur_ts;
  
  static int cur_win;
  
  static int num_plots;
  
  static int step_code;
  static int zoom_state;
  static int zoom_factor;
  
  static int jstep;
  
  static float xzoom_i;
  static float yzoom_i;
  static float xzoom_j;
  static float yzoom_j;
  static float xzoom_k;
  static float yzoom_k;
  
  static float xshift_i;
  static float yshift_i;
  static float xshift_j;
  static float yshift_j;
  static float xshift_k;
  static float yshift_k;
  
  static float table_last_time;

  static std::map<int,WindowVariable> map_windowVariables;
  static std::map<int,MfixVariable>   map_mfixVariables;
  
  static void init();  
  static void SetVariableRanges();


 
  
  
  static int  CalculateColor(int win , int ijk , float vmin , float vmax);
  static int  CalculateColor_v2(int ijk , float vmin , float vmax);
  static void Plot_IJ_slice(int win);
  static void Set_IJ_ortho(int win);
  static void Set_Viewport(int win);
  static void Set_Viewport_title(int win , int & x , int & y);
  static void Set_Viewport_time(int & x , int & y);
  static void DrawTitle(int win);
  static void DrawTime();
  static void Draw_IJ_vectors(int win , std::vector<float> & uvel , std::vector<float> & vvel);
  
  static int  FindVariableIndex(const char * var_name);
  static void Update_vWindowVariable();
  static void ReadVariables();

  
  
  
};

#endif

