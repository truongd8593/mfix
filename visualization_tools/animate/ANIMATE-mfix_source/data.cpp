
#include <QtOpenGL>

#include <iostream>
#include <algorithm>
#include <string>
#include <cmath>
#include <set>

#include "data.h"

bool Data::bSlice_only = false;
bool Data::bPaused = true;
bool Data::bTable  = false;
bool Data::bBusy   = false;
bool Data::bScaleChanged = false;


int  Data::gl_width = 10;
int  Data::gl_height = 10;
int  Data::npixels_x = 10;
int  Data::npixels_y = 10;
int  Data::cur_ts = 0;
int  Data::cur_win = 0;
int  Data::num_plots = 1;
int  Data::step_code = 0;
int  Data::zoom_state = 0;  // 0 = X&y .. 1=X ... -1=Y 
int  Data::zoom_factor = 1;

int  Data::jstep = 1;

float Data::xzoom_i = 1;
float Data::yzoom_i = 1;
float Data::xzoom_j = 1;
float Data::yzoom_j = 1;
float Data::xzoom_k = 1;
float Data::yzoom_k = 1;

float Data::xshift_i = 0;
float Data::yshift_i = 0;
float Data::xshift_j = 0;
float Data::yshift_j = 0;
float Data::xshift_k = 0;
float Data::yshift_k = 0;

float Data::table_last_time;

MyFont Data::myfont;

MfixData Data::mfix;

std::vector<PlotWindowInfo>   Data::vPlotWindowInfo;
std::vector<VariableRange>    Data::vVariableRange;
std::vector<MyTable>          Data::vTable;
std::vector<WindowVariable>   Data::vWindowVariable;
std::map<int,WindowVariable>  Data::map_windowVariables;
std::map<int,MfixVariable>    Data::map_mfixVariables;


namespace
{
  const int title_height = 75;
  const int time_height  = 75;
  
  float rra[] = { 1.00 , 0.95 , 0.90 , 0.90 , 0.90 , 0.90 , 
       0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 ,  
       0.75 , 0.65 , 0.55 , 0.50 , 0.40 , 0.40 , 0.40 , 
       0.35 , 0.35 , 0.30 , 0.30 , 0.30 , 0.30 , 0.20 , 
       0.20 , 0.20 , 0.10 , 0.00 , 0.00 , 1.00 }; 
 
   float gga[] = { 0.00 , 0.10 , 0.15 , 0.40 , 0.55 , 0.65 , 
       0.75 , 0.80 , 0.85 , 0.90 , 0.90 , 0.90 , 0.90 ,  
       0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 ,  
       0.90 , 0.90 , 0.90 , 0.85 , 0.80 , 0.75 , 0.65 , 
       0.55 , 0.40 , 0.15 , 0.10 , 0.00 , 1.0 };
 
   float bba[] = { 0.00 , 0.00 , 0.10 , 0.20 , 0.20 , 0.20 , 
       0.30 , 0.30 , 0.30 , 0.30 , 0.35 , 0.35 , 0.40 , 
       0.40 , 0.40 , 0.50 , 0.55 , 0.65 , 0.75 , 0.90 ,  
       0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 ,  
       0.90 , 0.90 , 0.90 , 0.95 , 1.00 , 1.0};
}


void Data::init()
{
   bSlice_only = false;
   bPaused    = true;
   bTable     = false;
   bScaleChanged = false;

   table_last_time = -2.0f;
   gl_width   = 10;
   gl_height  = 10;
   cur_ts     = 0;
   num_plots  = 1;
   cur_win    = 0;
   zoom_state = 0;
   xzoom_i     = 1;
   yzoom_i     = 1;
   xzoom_j     = 1;
   yzoom_j     = 1;
   xzoom_k     = 1;
   yzoom_k     = 1;
   zoom_factor = 1;
   
   jstep = 1;
   
   xshift_i = 0;
   yshift_i = 0;
   xshift_j = 0;
   yshift_j = 0;
   xshift_k = 0;
   yshift_k = 0;
   
   
   //  myfont.Read("pfSFuturaMono.c");
   //  myfont.Debug("font_debug.txt");
   myfont.Process();
}





int Data::CalculateColor(int win , int ijk , float vmin , float vmax)
{
  int icol = 31;
  
  if (vPlotWindowInfo[win].bContourFilled)
  {
    float rcol = 0.5f + 31.0f * (mfix.var[ijk] - vmin) / (vmax-vmin);
    icol = rcol;
    
    if (icol > 31) icol = 31;
    if (icol <  0) icol = 0;
  }
 
  return icol;
}

int Data::CalculateColor_v2(int ijk , float vmin , float vmax)
{
  int icol = 31;
  
  float rcol = 0.5f + 31.0f * (mfix.var[ijk] - vmin) / (vmax-vmin);
  icol = rcol;
  
  if (icol > 31) icol = 31;
  if (icol <  0) icol = 0;
 
  return icol;
}

void Data::Set_IJ_ortho(int win)
{  
   float z1 = -1;
   float z2 =  1;
   
   float x1 = 0.0f;
   float xtot = mfix.xlength;
   if (vPlotWindowInfo[win].bFlipX) xtot *= 2.0f;
   float xcen,x2,xhalf;
   
   float y1 = 0.0f;
   float y2 = mfix.ylength;
   
  // float ratio_screen = static_cast<float>(gl_width) / static_cast<float>(gl_height);
   float ratio_screen = static_cast<float>(npixels_x) / static_cast<float>(npixels_y);
   float ratio_world  = xtot / mfix.ylength;
   float ratio        = ratio_screen / ratio_world;
   
   if (ratio >= 1.0f)
   {
     y2 = mfix.ylength;
     float xhalf = 0.5f * ratio * xtot;
     if (vPlotWindowInfo[win].bFlipX)
       xcen = 0.0f;
     else
       xcen = 0.5 * xtot;
     
     x1 = xcen - xhalf;
     x2 = xcen + xhalf;
   }
   else
   {
     xhalf = 0.5f * xtot;
     if (vPlotWindowInfo[win].bFlipX)
       xcen = 0.0f;
     else
       xcen = 0.5 * xtot;
     
     x1 = xcen - xhalf;
     x2 = xcen + xhalf;
     
     y1 = 0.0f;
     y2 = mfix.ylength / ratio;
   }
   
   float ycen  = 0.5f * (y1+y2);
   float xdiff = x2 - x1;
   float ydiff = y2 - y1;
   
   x1 = xcen - 0.5f * Data::xzoom_k * xdiff;
   x2 = xcen + 0.5f * Data::xzoom_k * xdiff;
   y1 = ycen - 0.5f * Data::yzoom_k * ydiff;
   y2 = ycen + 0.5f * Data::yzoom_k * ydiff;
   
   y1 = y1 + yshift_k * yzoom_k * ydiff;
   y2 = y2 + yshift_k * yzoom_k * ydiff;
   
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   
   glOrtho(x1,x2,y1,y2,z1,z2);
      
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

void Data::Set_Viewport(int win)
{
  float tx = gl_width;
  float xnum = Data::num_plots; //vPlotWindowInfo.size();
  
  float xw = (0.95f - (xnum-1.0f)*0.05f) * tx / xnum;
  npixels_x = xw;
  npixels_y = gl_height - title_height - time_height;
    
  int ix1 = win*(npixels_x+gl_width/20)+gl_width/40;
  
  glViewport(ix1,0,npixels_x,npixels_y);
}

void Data::Set_Viewport_title(int win , int & ix , int & iy)
{
  float tx = gl_width;
  float xnum = Data::num_plots; //vPlotWindowInfo.size();
  
  float xw = (0.95f - (xnum-1.0f)*0.05f) * tx / xnum;
  int npixels_x = xw;
  int npixels_y = gl_height - time_height - title_height;
    
  int ix1 = win*(npixels_x+gl_width/20)+gl_width/40;
  
//  glViewport(ix1,npixels_y,npixels_x,title_height);
  glViewport(ix1,npixels_y,npixels_x,title_height);
  
  ix = npixels_x;
  iy = title_height;
}

void Data::Set_Viewport_time(int & ix , int & iy)
{
  int npixels_y = gl_height - title_height;
    
  glViewport(0,npixels_y,gl_width,title_height);
  
  ix = gl_width;
  iy = title_height;
}


void Data::Plot_IJ_slice(int win)
{  
   Set_Viewport(win);
   Set_IJ_ortho(win);
  
   size_t ijk = vPlotWindowInfo[win].win_k * mfix.imax2 * mfix.jmax2;
   
   map<int,WindowVariable>::iterator itf = map_windowVariables.find( vPlotWindowInfo[win].win_var   );
   map<int,MfixVariable>::iterator   itm =   map_mfixVariables.find( vPlotWindowInfo[win].win_var   );
	  
   if (itf == map_windowVariables.end() || itm==map_mfixVariables.end())
   {
      cout << "error # 1\n";
      return;
   }
   
   
 ///////////  float vmin = vVariableRange[ vPlotWindowInfo[win].win_var ].vmin;
 //////////  float vmax = vVariableRange[ vPlotWindowInfo[win].win_var ].vmax;
   
   double y1 = -mfix.DY[0];
   for (int j=0; j<mfix.jmax2; j=j+jstep)
   {
     double x1 = -mfix.DX[0];
     
     double yadd = 0;
     for (int jj=j; jj<j+jstep; ++jj)
     {
	if (jj < mfix.jmax2) yadd += mfix.DY[jj];
     }
          
     for (int i=0; i<mfix.imax2; ++i)
     {
       if (mfix.FLAG[ijk] == 1)
       {	
	 
	 if (ijk >= itm->second.colors.size())
	 {
	   cout << "error # 2\n";
	   cout << ijk << " : " << itm->second.colors.size() << "\n";
	   return;
	 }
	 
	  int icol = itm->second.colors[ijk];
        
          glBegin(GL_POLYGON);
       
          glColor3f(rra[icol],gga[icol],bba[icol]);
	  
       
          glVertex3f(x1            , y1      , 0.0f);
          glVertex3f(x1+mfix.DX[i] , y1      , 0.0f);
          glVertex3f(x1+mfix.DX[i] , y1+yadd , 0.0f);
          glVertex3f(x1            , y1+yadd , 0.0f);
       
          glEnd();
	  
	  
	  if (vPlotWindowInfo[win].bFlipX)
	  {
	    glBegin(GL_POLYGON);
       
       
	    glVertex3f(-x1            , y1            , 0.0f);
	    glVertex3f(-x1-mfix.DX[i] , y1            , 0.0f);
	    glVertex3f(-x1-mfix.DX[i] , y1+mfix.DY[j] , 0.0f);
	    glVertex3f(-x1            , y1+mfix.DY[j] , 0.0f);
       
	    glEnd();
	    
	    
	  }
	  glFlush();
       }
       
       x1 += mfix.DX[i];
       
       ++ijk;
     }
     y1 += yadd;
     ijk = ijk + (jstep-1)*mfix.imax2;
   }
}

void Data::DrawTitle(int win)
{
   int ix,iy;
  
   Set_Viewport_title(win,ix,iy);

   int NY = 80000;
   int NX = NY * ix / iy;
   
   if (NX < 500000)
   {
     NY *= 500000/NX;
     NX = 500000;
   }
   
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   
   glOrtho(0,NX,0,NY,-1,1);
         
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
  
   if (win == Data::cur_win)
     glColor3f(0.4f,0.4f,0.4f);
   else
     glColor3f(0.2f,0.2f,0.2f);
   
   glBegin(GL_POLYGON);
      glVertex3f(0,0,0);
      glVertex3f(NX,0,0);
      glVertex3f(NX,NY,0);
      glVertex3f(0,NY,0);
   glEnd();
   
   
   glColor3f( 1.0f , 1.0f , 1.0f );
   
  
    glLineWidth(2.0);
    string s =  mfix.variable_names[ vPlotWindowInfo[win].win_var ];
    if (s.size() < 10)
    {
      s.insert(0,(12-s.size())/2,' ');
    }
    myfont.DrawString(10,15000, s.c_str());
    glLineWidth(1.0);
    
    glFlush();
}

void Data::DrawTime()
{
   int ix,iy;
  
   Set_Viewport_time(ix,iy);
   
   const int NY = 80000;
   const int NX = NY * ix / iy;
  
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   
   glOrtho(0,NX,0,NY,-1,1);
      
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
  
   glColor3f(0.1f,0.1f,0.1f);
   
   glBegin(GL_POLYGON);
      glVertex3f(0,0,0);
      glVertex3f(NX,0,0);
      glVertex3f(NX,NY,0);
      glVertex3f(0,NY,0);
   glEnd();
   
 //  qgl->qglColor( QColor(255,255,255) );
 
   glColor3f(1.0,1.0,1.0);
   
   std::string s = Convert<float,std::string>(mfix.times[cur_ts]);
   myfont.DrawString(NX/3,15000,s.c_str());
   glFlush();

  
}


void Data::SetVariableRanges()
{
   vVariableRange.resize( mfix.Get_NVARS() );
   
   for (size_t i=0; i<vVariableRange.size(); ++i)
   {
      mfix.GetVariableAtTimestep(i,mfix.Get_NTIMES()-1);
      
      float vmin = 1.0e+32;
      float vmax = -vmin;
      
      for (int ijk=0; ijk<mfix.ijkmax2; ++ijk)
      {
	if (mfix.FLAG[ijk] == 1)
	{
	  if (mfix.var[ijk] < vmin) vmin = mfix.var[ijk];
	  if (mfix.var[ijk] > vmax) vmax = mfix.var[ijk];
	}
      }
      
      if (vmin >= vmax) vmax = vmin + 1.0f;
      
      vVariableRange[i].vmin = vmin;
      vVariableRange[i].vmax = vmax;
   }
  
}


 
PlotWindowInfo::PlotWindowInfo() :
     bIJ_plot(true) , 
     bIK_plot(false) ,
     bContourFilled(true) , 
     bContourLines(false) ,
     bPlotGasVectors(false) ,
     bPlotSolidVectors(false) ,
     bPlotGrid(false) ,
     bFlipX(false) ,
     bPlot(true) , 
     win_var(0) ,
     win_i(0) ,
     win_j(0) , 
     win_k(0) ,
     win_imin(0) ,
     win_jmin(0) ,
     win_kmin(0) ,
     win_imax(0) ,
     win_jmax(0) ,
     win_kmax(0) ,
     win_m(0) ,
     win_ng(0) ,
     win_ns(0) ,
     win_scalar(0) , 
     win_rr(0) ,
     contour_line_min(0.4f) ,
     contour_line_max(1.0f) ,
     n_contours(5)
     {
       if (Data::mfix.kmax2 > 1) win_k = 1; 
     }
     
     
void Data::Draw_IJ_vectors(int win , std::vector<float> & uvel , std::vector<float> & vvel)
{
      const float pi     = 3.14159f;
      const float two_pi = 2.0f*pi;
      
      const float y_factor = 1.0f;
            float x_factor = 1.0f;

      
   /////////////////////////////////   float xlen = mfix.DX[0] + mfix.xlength + mfix.DX[mfix.imax2-1];
      float ylen = mfix.DY[0] + mfix.ylength + mfix.DY[mfix.jmax2-1]; 
      
      float y1   = -mfix.DY[0];
      float y2   = (mfix.ylength + mfix.DY[mfix.jmax2-1])*y_factor; 

      float xpixels = npixels_x;
      float ypixels = npixels_y;
      
      float ratio = xpixels / ypixels;
      
      bool true_scale = true;
      
      
      if (true_scale) x_factor = y2 / ylen;
      
      
      float xtot = ylen*ratio;
      
      float x1 = -mfix.DX[0];
      
      float x2 = mfix.xlength + mfix.DX[mfix.imax2-1];
      
      float xcen = 0.0f;
      
      if (!vPlotWindowInfo[win].bFlipX) 
         xcen = 0.5 * (x1 + x2);
      
      
      x1 = xcen - 0.5f*xtot*x_factor;
      x2 = xcen + 0.5f*xtot*x_factor;
      
      float d_tran_x = 0.05f;
      float d_tran_y = 0.05f;
      
      float n_tran_x = 0.0f;
      float n_tran_y = 0.0f;
      
      
      float dtx = d_tran_x * (x2-x1) * n_tran_x; 
      float dty = d_tran_y * (y2-y1) * n_tran_y;
      float y2_last = y2 + dty;
      
      float xview_len = x2 - x1; 
      float yview_len = y2 - y1;
      
      float xview_cen = -999;
      float yview_cen = -999;
      
      
      if (xview_cen<-990.0f &&  yview_cen<-990.0f)
      {
         xview_cen = (x1+dtx) + 0.5*xview_len;
         yview_cen = (y1+dty) + 0.5*yview_len;
      }
      
      x1 = 0.1f*n_tran_x*xview_len + xview_cen - 0.5f*xview_len;
      y1 = 0.1f*n_tran_y*yview_len + yview_cen - 0.5f*yview_len;
      x2 = 0.1f*n_tran_x*xview_len + xview_cen + 0.5f*xview_len;
      y2 = 0.1f*n_tran_y*yview_len + yview_cen + 0.5f*yview_len;
      
      
            y2_last = y2;  
  //////////////////    float x1_org  = x1; 
 ////////////////////     float x2_org  = x2; 
 /////////////////////     float y1_org  = y1; 
 ////////////////////     float y2_org  = y2; 
      

      float dist1 = mfix.xlength / static_cast<float>(mfix.imax2);
      float dist2 = mfix.ylength / static_cast<float>(mfix.jmax2);
      
      bool scale_vectors = true;
      
      float dist;
      float vec_factor = 1.0f;
      float mag_scale = 0.01f;
      
      if (scale_vectors)
         dist  = mag_scale * vec_factor * std::min(dist1,dist2); 
      else
         dist  = std::min(dist1,dist2)  * vec_factor;// ! * mag_scale ;
      
 
 /*
c 
            if (color_code.eq.0 .or. colormode.eq.1 .or.
     &               colormode.eq.2) then 
               if (background.eq.0) then 
                  z1 = 1.0 
               else 
                  z1 = 0.0 
               end if 
               call set_color(z1,z1,z1) 
            else if (color_code.eq.1) then 
               if (background.eq.0) then 
                  z1 = 1.0 
                  z2 = 0.9 
                  z3 = 1.0 
               else 
                  z1 = 1.0 
                  z2 = 0.1 
                  z3 = 0.0 
               end if 
               call set_color(z1,z2,z3) 
            end if 
	    */
 
      glColor3f(1.0f,1.0f,1.0f);
 
      glLineWidth(2.0f);
 
 
 
      int count = vPlotWindowInfo[win].win_k*mfix.imax2*mfix.jmax2 - 1;
      y1 = -mfix.DY[0];
      
      float xc,yc,mag2,mag,theta,vector_scale,theta1,xv,yv;
      float th , dx2 , vel , uel , ydx , ydxa , phi;
      float xt,yt , phmth,phpth , xr1,yr1 , xr2,yr2 , x3,y3 , x4,y4;
      float ddist , ddis1;
      
      int nvec_x = 1;
      int nvec_y = 1;
            
      for (int j=0; j<mfix.jmax2; ++j)
      {
         x1 = -mfix.DX[0]; 
         yc = y1 + 0.5f*mfix.DY[j];
	 
         for (int i=0; i<mfix.imax2; ++i)
	 {
            ++count; 
            if (i%nvec_x != 0) goto L111;
            if (j%nvec_y != 0) goto L111;
            if (mfix.FLAG[count] ==  1)   goto L112;
            if (mfix.FLAG[count] == 10)   goto L112;
            if (mfix.FLAG[count] == 20)   goto L112; 
            goto L111; 
L112:       xc = x1 + 0.5f*mfix.DX[i];

            mag2 = uvel[count]*uvel[count] + vvel[count]*vvel[count];
	    mag  = sqrt(mag2);
	    
            if (scale_vectors)
               vector_scale = mag;  
            else 
               vector_scale = 1.0f;

            if (mag == 0.0f) goto L111;

            if (uvel[count] == 0.0f)
	    {
               if (vvel[count] > 0.0f) theta = 0.5f * pi; 
               if (vvel[count] < 0.0f) theta = 3.5f * pi; 
               goto L222;
	    }
	    
	    
            if (vvel[count] == 0.0f)
	    {
               if (uvel[count] > 0.0f) theta = 0.0f;
               if (uvel[count] < 0.0f) theta = pi; 
               goto L222;
	    }

            theta1 = atan(fabs(vvel[count]/uvel[count]));
	    
            if (uvel[count] > 0.0 && vvel[count]>0.0) theta = theta1;
            if (uvel[count] < 0.0 && vvel[count]>0.0) theta = pi - theta1;
            if (uvel[count] < 0.0 && vvel[count]<0.0) theta = pi + theta1;
            if (uvel[count] > 0.0 && vvel[count]<0.0) theta = two_pi - theta1;

L222:

            if (scale_vectors)
	    {
               xv = xc+dist*vec_factor*cos(theta)* vector_scale;
               yv = yc+dist*vec_factor*sin(theta)* vector_scale;
	    }
            else 
	    {
               xv = xc+dist*cos(theta)* vector_scale;
               yv = yc+dist*sin(theta)* vector_scale;
	    }
 
            glBegin(GL_LINES);
               glVertex3f(xc,yc,0.0);
               glVertex3f(xv,yv,0.0);
           glEnd(); 

// sym_x start 
            if (vPlotWindowInfo[win].bFlipX)
	    {
               glBegin(GL_LINES);
                  glVertex3f(-xc,yc,0.0); 
                  glVertex3f(-xv,yv,0.0);
               glEnd();
	    }
// sym_x end 


            th = 0.25f * pi;
            if (scale_vectors)
               dx2 = 0.15f*vec_factor*vector_scale*dist ;
            else 
               dx2 = 0.15f*vector_scale*dist ;
   
            vel = yv - yc;
            uel = xv - xc; 
            if (uel != 0.0)
	    {
               ydx = vel/uel ;
               ydxa = fabs(ydx) ;
               phi = atan(ydx) ;
	    }
            else
	    {
               if (vel < 0.0) phi = -0.5*pi ;
               if (vel > 0.0) phi = +0.5*pi ;
	    }
	    
            yt = yv ;
            xt = xv ;
            phmth = phi-th ;
            phpth = phi+th ;
            xr1 = dx2*cos(phmth) ;
            yr1 = dx2*sin(phmth) ;
            xr2 = dx2*cos(phpth) ;
            yr2 = dx2*sin(phpth) ;
            x3 = xt - xr1 ;
            y3 = yt - yr1 ;
            x4 = xt - xr2 ;
            y4 = yt - yr2 ;
	    
	    
            ddist = (xc-xt)*(xc-xt)+(yc-yt)*(yc-yt);
            ddis1 = (xc-x3)*(xc-x3)+(yc-y3)*(yc-y3);
            if (ddis1 >= ddist)
	    {
               x3 = xt + xr1 ;
               y3 = yt + yr1 ;
	    }
	    
            ddis1 = (xc-x4)*(xc-x4)+(yc-y4)*(yc-y4);
            if (ddis1 >= ddist)
	    {
               x4 = xt + xr2 ;
               y4 = yt + yr2 ;
	    }
            glBegin(GL_LINE_STRIP);
               glVertex3f(xt,yt,0.0);
               glVertex3f(x3,y3,0.0); 
               glVertex3f(x4,y4,0.0);
               glVertex3f(xt,yt,0.0);
            glEnd();

// sym_x start 

            if (vPlotWindowInfo[win].bFlipX)
	    {
               glBegin(GL_LINE_STRIP);
                 glVertex3f(-xt,yt,0.0) ;
                 glVertex3f(-x3,y3,0.0) ;
                 glVertex3f(-x4,y4,0.0) ;
                 glVertex3f(-xt,yt,0.0) ;
               glEnd();
	    }

// sym_x end 

 
L111:
 
 
            x1 += mfix.DX[i];
	 }
         y1 += mfix.DY[j];
         if (y1 > y2_last) goto end_loop; 
      }
      
      
end_loop:

      glLineWidth(1.0f);
      return;
}

int Data::FindVariableIndex(const char * var_name)
{
  using namespace std;
  
  vector<string>::iterator it = find(mfix.variable_names.begin(),mfix.variable_names.end(),var_name);
  
  if (it == mfix.variable_names.end())
  {
    return -1;
  }
  else
  {
    return distance( mfix.variable_names.begin() , it );
  }
}

struct VariableEquals
{
  VariableEquals(int vi) : variable_index(vi) {}

  bool operator () (const WindowVariable & wv) const
  {
    return wv.mfix_variable_index == variable_index;
  }
  
  int variable_index;
  
};

struct VariableUsed
{
  bool operator () (const WindowVariable & wv) const
  {
    return !wv.bUsed;
  }  
};


void Data::Update_vWindowVariable()
{
  using namespace std;  
  
  map<int,WindowVariable> newMap;
  
   map<int,MfixVariable>::iterator mit = map_mfixVariables.begin();
   while (mit != map_mfixVariables.end())
   {
     mit->second.bUsed = false;
     mit->second.bMustCalculateColor = false;
     ++mit;
   }
  
  set<int> sColorVariables;
  set<int> sNonColorVariables;
  
  for (int  i=0; i<num_plots; ++i)
  {
      sColorVariables.insert( vPlotWindowInfo[i].win_var );
      if (vPlotWindowInfo[i].bPlotGasVectors)
      {
	 int component_index = FindVariableIndex("U_g");
	 if  (component_index >=0) sNonColorVariables.insert(component_index);
	     component_index = FindVariableIndex("V_g");
	 if  (component_index >=0) sNonColorVariables.insert(component_index);
      }
  }
  
  set<int>::iterator sit = sColorVariables.begin();
  while (sit != sColorVariables.end())
  {
      mit = map_mfixVariables.find(*sit);
      if (mit != map_mfixVariables.end())
      {
	  mit->second.bUsed = true;
	  mit->second.bMustCalculateColor = true;
      }
      else
      {
	  // add variable to map_mfixVariables
	  MfixVariable mv;
	  mv.bUsed = true;
	  mit->second.bMustCalculateColor = true;
	  map_mfixVariables[*sit] = mv;
      }
      ++sit;
  }
  
  sit = sNonColorVariables.begin();
  while (sit != sNonColorVariables.end())
  {
      mit = map_mfixVariables.find(*sit);
      if (mit != map_mfixVariables.end())
      {
	  mit->second.bUsed = true;
      }
      else
      {
	  // add variable to map_mfixVariables
	  MfixVariable mv;
	  mv.bUsed = true;
	  mit->second.bMustCalculateColor = false;
	  map_mfixVariables[*sit] = mv;
      }
      ++sit;
  }

  
  for (int  i=0; i<num_plots; ++i)
  {
      int color_variable_index = vPlotWindowInfo[i].win_var;

      WindowVariable wv;

      map<int,WindowVariable>::iterator it = newMap.find(color_variable_index);

      if (it != newMap.end())
      {
	wv = it->second;
      }
      else
      {
	  map<int,WindowVariable>::iterator itf = map_windowVariables.find(color_variable_index);
	  
	  if (itf == map_windowVariables.end())
	  {
	    wv.ts                   = -1; // this causes the need to read everything
	    wv.bMustCalculateColors = true;
	    wv.bColorsCalculated    = false;
	    wv.bUsed                = true;
	    wv.bMustRead            = true;
	    wv.mfix_variable_index  = color_variable_index;
	  }
	  else
	  {
	    wv = itf->second;
	  }
      }
      
      newMap[color_variable_index] = wv;
  }
      
  map<int,WindowVariable>::iterator it = newMap.begin();
  while (it != newMap.end())
  {
    map<int,WindowVariable>::iterator itf = map_windowVariables.find(it->first);
    
    if (itf == map_windowVariables.end())
    {
     	  it->second.bMustRead            = true;
	  it->second.ts                   = cur_ts;
	  it->second.bMustCalculateColors = true;
	  it->second.bColorsCalculated    = false;
	  it->second.bUsed                = true;
//	  it->second.mfix_variable_index  = color_variable_index;
    }
    else
    {
      if (it->second.ts != cur_ts)
      {
	  it->second.ts                   = cur_ts;
	  it->second.bMustCalculateColors = true;
	  it->second.bColorsCalculated    = false;
	  it->second.bUsed                = true;
	  it->second.bMustRead            = true;
      }
      else
      {
	  it->second.ts                   = cur_ts;
	  it->second.bMustCalculateColors = false;
	  it->second.bColorsCalculated    = true;
	  it->second.bUsed                = true;
	  it->second.bMustRead            = false;
      }
    }
 
    
    
    ++it;
  }
  
  map_windowVariables = newMap;
  
}


void Data::ReadVariables()
{
      map<int,WindowVariable>::iterator it = map_windowVariables.begin();
      while (it !=  map_windowVariables.end())
      {
	if (it->second.bMustRead || Data::bScaleChanged)
	{
	    it->second.bMustRead = false;
	    it->second.bColorsCalculated = true;
	    
	    int var = it->first;
	    
	   // it->second.values.resize(mfix.ijkmax2);
	    it->second.ts = cur_ts;
	    
	    map<int,MfixVariable>::iterator mit = map_mfixVariables.find(var);
	    if (mit != map_mfixVariables.end())
	    {
	     //   cout << "read\n";
		mfix.GetVariableAtTimestep(var,cur_ts);
		//cout << "set\n";
		mit->second.values = mfix.var;  // should have option to read directly into variable
		
		mit->second.colors.resize(mfix.ijkmax2);
	    
		float vmin = vVariableRange[ var ].vmin;
		float vmax = vVariableRange[ var ].vmax;
	    
		//cout << "colors\n";
		for (int ijk=0; ijk<mfix.ijkmax2; ++ijk)
		{
		  mit->second.colors[ijk] = CalculateColor_v2(ijk ,  vmin ,  vmax);
		}
		//cout << "finished calculating colors for var = " << var << "\n";
	    }
	    else
	    {
	        cout << "error : position ABC\n";
	    }
	}

        Data::bScaleChanged = false;
	
	++it;
      }

}

