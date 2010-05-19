#include <QtGui>
#include <QtOpenGL>
#include <QMessageBox>
#include <QMenu>


#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>

#include "tetrahedron.h"
#include "data.h"

namespace
{
   Tetrahedron * tetra_window = 0;
   
   int image_number = 0;
   
//////////////   bool bRightButtonMenu = false;
   
   int draw_counter = 0;

 /////////////  bool bFirst = true;
   
   GLUquadricObj * quadric;
   
   GLfloat mat_specular[]   = { 0.5f , 0.5f , 0.5f , 0.5f };
   GLfloat mat_ambient[]    = { 100.0f/255.0f , 100.0f/255.0f , 100.0f/255.0f , 0.5f };
   GLfloat mat_shininess[]  = { 10.0f };
   GLfloat light_position[] = { 130.0f , 75.0f , -130.0f , 0.0f };
  ////////// GLfloat mat_emmision[]   = { 0.6f , 0.2f , 0.2f };
  ////////// GLfloat mat_emmision2[]  = { 0.2f , 0.2f , 0.6f };
   GLfloat mat_emmision_trace[]   = { 0.6f , 0.2f , 0.2f };

    const int mat_size = 4;

   GLfloat mat[mat_size][4] = { 
				{ 0.6f , 0.2f , 0.2f } , 
				{ 0.2f , 0.2f , 0.6f } ,
				{ 0.2f , 0.6f , 0.2f } ,
				{ 0.8f , 0.6f , 0.4f } 
			      };
}

Tetrahedron::Tetrahedron(QWidget *parent)
    : QGLWidget(parent)
{
    setFormat(QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer)); // | QGL::DepthBuffer));
   
    Data::bPaused = true;

    tetra_window = this;
    
}

void Tetrahedron::togglePaused()
{
   if (Data::bPaused)
      killTimer(timerID);
   else
      timerID = startTimer(200);
}



extern void RightClick_variable_change(int var);

void Tetrahedron::contextMenuEvent(QContextMenuEvent * /* event */)
{
  return;
}



void Tetrahedron::timerEvent(QTimerEvent * event)
{
   if (Data::bBusy) return;
   
   Data::ReadNext_XML_file();
  
   if (event != 0) updateGL();
}


void Tetrahedron::initializeGL()
{
    qglClearColor(Qt::black);
    glShadeModel(GL_SMOOTH);
    
    quadric = gluNewQuadric();
    gluQuadricNormals(quadric,GL_SMOOTH);
    
    glMaterialfv(GL_FRONT,GL_SPECULAR,mat_specular);
    glMaterialfv(GL_FRONT,GL_SHININESS,mat_shininess);
    glLightfv(GL_LIGHT0,GL_POSITION,light_position);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,mat_ambient);
}

void Tetrahedron::resizeGL(int width, int height)
{
    Data::bBusy = true;
    Data::gl_width  = width;
    Data::gl_height = height;
    glViewport(0, 0, width, height);
    Data::bBusy = false;
}

void Tetrahedron::paintGL()
{
    if (Data::bBusy) return;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // | GL_DEPTH_BUFFER_BIT);
    draw();
}

extern void CurrentWindowChanged(int cw);

void Tetrahedron::mousePressEvent(QMouseEvent *  /* event */)
{  
    return;
}


void Tetrahedron::keyPressEvent(QKeyEvent * /* event */)
{  
}

void DrawCircle(float x , float y , float d)
{
  float r = d/2.0f;
  
  glBegin(GL_POLYGON);
  
  for (int i=0; i<=360; i=i+20)
  {
    float radians = 3.14159f * (float)i / 180.0f;
    float xp = x + r * cos(radians);
    float yp = y + r * sin(radians);
    glVertex3f(xp,yp,0.0f);
  }
  glEnd();
}


extern void UpdateTimeSlider(int plot_ts);

void Tetrahedron::draw()
{    
   glViewport(0,0,Data::gl_width,Data::gl_height-75);
    
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   double yy1 = -0.075 * Data::ymax;
   double yy2 = Data::ymax - yy1;
   double ydiff = yy2 - yy1;
   
   double xdiff = ydiff * static_cast<float>(Data::gl_width)/static_cast<float>(Data::gl_height-75);
   double xcen  = 0.5 * (Data::xmin+Data::xmax);

   double xx1 = xcen - 0.5*xdiff;
   double xx2 = xcen + 0.5*xdiff;
   
   double dx = (xx2-xx1)*Data::xtrans*Data::scale;
   double dy = (yy2-yy1)*Data::ytrans*Data::scale;
   
    
   glOrtho(xx1+dx,xx2+dx,yy1+dy,yy2+dy,-1000,1000);
//   cout << xx1+dx << " : " << xx2+dx << " : " << yy1+dy
 //        << " : " << yy2+dy << "\n";
   
   
   glTranslatef(xcen , Data::ymax/2 , 0.0);
   
   glScalef(Data::scale,Data::scale,Data::scale);
   

   glRotatef(Data::xrot,1.0f,0.0f,0.0f);
   glRotatef(Data::yrot,0.0f,1.0f,0.0f);

   glTranslatef(-xcen , -Data::ymax/2 , 0.0);
     
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();


/*
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   
   double dx = (Data::xmax - Data::xmin) * 0.05;
   double dy = (Data::ymax - Data::ymin) * 0.05;
   
   double ylen = Data::ymax - Data::ymin + 2.0*dy;
   
   double xdiff = 0.5*ylen*static_cast<float>(Data::gl_width)/static_cast<float>(Data::gl_height-75);
   double xcenter = 0.5 * (Data::xmin + Data::xmax);
   double xx1 = xcenter - xdiff;
   double xx2 = xcenter + xdiff;
   
   xx1 += Data::xtrans;
   xx2 += Data::xtrans;
    
   glOrtho(xx1,xx2,Data::ymin-dy,Data::ymax+dy,-10000,10000);
   
   float xcenter_t = 0; //0.5 * ( Data::xmin + Data::xmax);
   float ycenter_t = -0.5 * ( Data::ymin + Data::ymax);
   
 //  glTranslatef(xcenter_t , ycenter_t , 0.0);
   
   glScalef(Data::scale,Data::scale,Data::scale);
   

   glRotatef(Data::xrot,1.0f,0.0f,0.0f);
   glRotatef(Data::yrot,0.0f,1.0f,0.0f);

//   glTranslatef(-xcenter_t , -ycenter_t , 0.0);
     
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
       */

    if (Data::b3D_bubbles)
    {
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	glColorMaterial(GL_FRONT,GL_EMISSION);
	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);
	
	glLightfv(GL_LIGHT0,GL_POSITION,light_position);
    }
    
    for (int i=0,n=Data::vP.size(); i<n; ++i)
    {
	int color = 0;
        if (Data::color_method == 1) color = Data::vColors_position[i];
        if (Data::color_method == 2) color = Data::vColors_diameter[i];

        if (Data::color_method == 3)
	{
            double vm = Data::vP[i].u * Data::vP[i].u + Data::vP[i].v * Data::vP[i].v + Data::vP[i].w * Data::vP[i].w;
            double xcol = 0.01 + 3.0 * vm / Data::vel_max;
            color = xcol;
 	}

        if (Data::color_method == 4)
	{
            if (Data::vP[i].v >= 0)
	      color = 0;
	    else
	      color = 1;
 	}

        if (color >= mat_size) color = 0;

	if (Data::b3D_bubbles)
	{
	    glPushMatrix();

	    glMaterialfv(GL_FRONT,GL_EMISSION,mat[color]);
	    
	    glTranslatef(Data::vP[i].x , Data::vP[i].y , Data::vP[i].z);
	    gluSphere(quadric,Data::vP[i].d/2.0f,Data::sphere_resolution,Data::sphere_resolution);
	    glPopMatrix();
	}
	else
	{
	    glColor3f(mat[color][0] , mat[color][1] , mat[color][2] );

	    DrawCircle( Data::vP[i].x , Data::vP[i].y , Data::vP[i].d );
	}
    }
    
   
    if (Data::time != Data::bAdd_time)
    {
      if (Data::time < Data::bAdd_time) Data::vPos.clear();

      Data::vPos.push_back(Data::vP[0]);
    }

    Data::bAdd_time = Data::time;

    Data::bAdd_to_vpos = false;
    
    const int n = Data::vPos.size();
    int istart =  n - 20;
    if (istart < 0) istart = 0;

    if (Data::bTrace)
    {
	  for (int i=istart; i<n; ++i)
	  {
	    float blue = static_cast<float>(i-istart) / static_cast<float>(n-1);
      
	    if (Data::b3D_bubbles)
	    {
		mat_emmision_trace[0] = 1.0f;
		mat_emmision_trace[1] = 1.0f;
		mat_emmision_trace[2] = blue;
		glPushMatrix();
      
		glMaterialfv(GL_FRONT,GL_EMISSION,mat_emmision_trace);
      
      
		glTranslatef(Data::vPos[i].x , Data::vPos[i].y , Data::vPos[i].z);
		gluSphere(quadric,Data::vPos[i].d/2.0f,10,10);
		glPopMatrix();
	    }
	    else
	    {
		blue = 0.0f;
		glColor3f(1.0f,1.0f,blue);
		DrawCircle( Data::vPos[i].x , Data::vPos[i].y , Data::vPos[i].d);
	    }
	  }
	      
	  glDisable(GL_LIGHTING);
	  glDisable(GL_LIGHT0);
	  glDisable(GL_DEPTH_TEST);
	  glDisable(GL_AUTO_NORMAL);
	  glDisable(GL_NORMALIZE);
		
	  glLineWidth(2.0f);
	  glBegin(GL_LINE_STRIP);
	  for (int i=istart; i<n; ++i)
	  {
	    float blue = static_cast<float>(i) / static_cast<float>(n-1);
	    glColor3f(1.0f,1.0f,blue);
	    glVertex3f(Data::vPos[i].x , Data::vPos[i].y , Data::vPos[i].z);
	  }
	  glEnd();
    }
    
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);

    DrawOutline();

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_AUTO_NORMAL);
    glDisable(GL_NORMALIZE);

    DrawTime();
    UpdateTimeSlider( Data::cur_ts );
    if (!Data::bPaused)  ++Data::cur_ts;

    if (Data::bImage)
    {
	std::stringstream ss;
        ss << "DEM_MFIX-plot_" << std::setfill('0') << std::setw(5) << image_number << ".bmp";
	QImage qim = grabFrameBuffer();
	qim.save(QString(ss.str().c_str()),"BMP");
	++image_number;
        Data::bImage = false;
    }
      
    
    ++draw_counter;
}



void MyTogglePause()
{
   if (tetra_window != 0)
   {
      tetra_window->togglePaused();
   }
}

void StepGL()
{
  if (tetra_window != 0) tetra_window->updateGL();
}

void SendKey(int /* key */ )
{
}
   
void Tetrahedron::DrawTime()
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
    
   glColor3f(1.0,1.0,1.0);
   
   std::string s = Convert<double,std::string>(Data::time);
   Data::myfont.DrawString(NX/3,15000,s.c_str());
   glFlush();
}

void Tetrahedron::Set_Viewport_time(int & ix , int & iy)
{
  glViewport(0,Data::gl_height-75,Data::gl_width,75);
  
  ix = Data::gl_width;
  iy = 75;
}

void Tetrahedron::DrawOutline()
{
    if (Data::mfix_kmax2==1 && Data::mfix_bCartesian)   // 2D cartesian case
    {
      glColor3f(0.4f,0.2f,0.6f);
      glBegin(GL_LINE_STRIP);
	  glVertex3f(0.0f,0.0f,0.0f);
	  glVertex3f(Data::mfix_xlength,0.0f,0.0f);
	  glVertex3f(Data::mfix_xlength,Data::mfix_ylength,0.0f);
	  glVertex3f(0.0f,Data::mfix_ylength,0.0f);
	  glVertex3f(0.0f,0.0f,0.0f);
      glEnd();
    }
}
