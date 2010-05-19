#include <QtGui>
#include <QtOpenGL>
#include <QMessageBox>
#include <QMenu>


#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include "tetrahedron.h"
#include "data.h"

namespace
{
   Tetrahedron * tetra_window = 0;
   
   int image_number = 0;
   
   bool bRightButtonMenu = false;
   
   int draw_counter = 0;

   bool bFirst = true;
}

Tetrahedron::Tetrahedron(QWidget *parent)
    : QGLWidget(parent)
{
    setFormat(QGLFormat(QGL::DoubleBuffer)); // | QGL::DepthBuffer));

    Data::bPaused = true;

    tetra_window = this;
}

void Tetrahedron::togglePaused()
{
   if (Data::bPaused)
      killTimer(timerID);
   else
      timerID = startTimer(5);
}



extern void RightClick_variable_change(int var);

void Tetrahedron::contextMenuEvent(QContextMenuEvent * /* event */)
{
    using namespace std;
    
    bRightButtonMenu = true;
    
    QMenu * qmenu = new QMenu();
    
    for (size_t i=0; i<Data::mfix.variable_names.size(); ++i)
    {
      qmenu->addAction(Data::mfix.variable_names[i].c_str());
    }
    
    QAction * paction = qmenu->addAction("gas vectors");
    paction->setCheckable(true);
    bool checked = Data::vPlotWindowInfo[Data::cur_win].bPlotGasVectors;
    paction->setChecked(checked);
    
    
    paction = qmenu->addAction("solid vectors");
    paction->setCheckable(true);
    checked = Data::vPlotWindowInfo[Data::cur_win].bPlotSolidVectors;
    paction->setChecked(checked);
    
    QAction * action = qmenu->exec(QCursor::pos());
    
    if (action != 0)
    {
      string s = action->text().toStdString();
      vector<string>::iterator it = std::find(Data::mfix.variable_names.begin() , Data::mfix.variable_names.end(),s);
      
      if (it != Data::mfix.variable_names.end())
      {
	RightClick_variable_change( std::distance(Data::mfix.variable_names.begin() , it ) );
      }
      else if (s == "gas vectors")
      {
	Data::vPlotWindowInfo[Data::cur_win].bPlotGasVectors = !Data::vPlotWindowInfo[Data::cur_win].bPlotGasVectors;
      }
      else if (s == "solid vectors")
      {
	Data::vPlotWindowInfo[Data::cur_win].bPlotSolidVectors = !Data::vPlotWindowInfo[Data::cur_win].bPlotSolidVectors;
      }
    }
    
    delete qmenu;
    
    bRightButtonMenu = false;
}


void Tetrahedron::timerEvent(QTimerEvent * event)
{
  if (Data::bBusy) return;
  
   if (event != 0) updateGL();
}


void Tetrahedron::initializeGL()
{
    qglClearColor(Qt::black);
    glShadeModel(GL_FLAT);
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
    glClear(GL_COLOR_BUFFER_BIT); // | GL_DEPTH_BUFFER_BIT);
    draw();
}

extern void CurrentWindowChanged(int cw);

void Tetrahedron::mousePressEvent(QMouseEvent *event)
{  
  Data::bBusy = true;
     
  if (event == 0)
  {
     Data::bBusy = false;
     return;
  }
  
  if (event->button() == Qt::LeftButton)
  {
    if (Data::num_plots < 2) 
    {
      Data::bBusy = false;
      return;
    }
  
    setFocus();
  
    float x = 0.1f + Data::num_plots * static_cast<float>(event->x()) / static_cast<float>(Data::gl_width);
    int ix = x;
  
    if (ix == Data::cur_win) 
    {
      Data::bBusy = false;
      return;
    }
  
    CurrentWindowChanged(ix);  
    Data::bBusy = false;    
    updateGL();
  }
  Data::bBusy = false;
}
/*
void Tetrahedron::mouseMoveEvent(QMouseEvent *event)
{
 
    GLfloat dx = GLfloat(event->x() - lastPos.x()) / width();
    GLfloat dy = GLfloat(event->y() - lastPos.y()) / height();

    if (event->buttons() & Qt::LeftButton) {
        rotationX += 180 * dy;
        rotationY += 180 * dx;
        updateGL();
    } else if (event->buttons() & Qt::RightButton) {
        rotationX += 180 * dy;
        rotationZ += 180 * dx;
        updateGL();
    }
    lastPos = event->pos();
 
}*/

// void Tetrahedron::mouseDoubleClickEvent(QMouseEvent *event)
// {
//  setFocus();
// }

void Tetrahedron::keyPressEvent(QKeyEvent * event)
{  
   switch (event->key())
   {
     case 's' :
     case 'S' :
       
       Data::step_code = 1;
       updateGL();
       break;
       
     case 'x' :
     case 'X' :
       
      if (Data::zoom_state != -1) 
       {
	 if (Data::zoom_factor == 1)
	    Data::xzoom_k /= 0.9f;
	 else
	    Data::xzoom_k /= 0.5f;
       }
       if (Data::zoom_state !=  1) 
       {
	 if (Data::zoom_factor == 1)
	    Data::yzoom_k /= 0.9f;
	 else
	    Data::yzoom_k /= 0.5f;
       }
       updateGL();
       break;
       
     case 'z' :
     case 'Z' :
       
       if (Data::zoom_state != -1) 
       {
	 if (Data::zoom_factor == 1)
	    Data::xzoom_k *= 0.9f;
	 else
	    Data::xzoom_k *= 0.5f;
       }
       if (Data::zoom_state !=  1) 
       {
	 if (Data::zoom_factor == 1)
	    Data::yzoom_k *= 0.9f;
	 else
	    Data::yzoom_k *= 0.5f;
       }
       updateGL();
       break;
       
     case 'r' :
     case 'R' :
       
       Data::xzoom_k  = 1;
       Data::yzoom_k  = 1;
       Data::yshift_k = 0;

       updateGL();
       break;
       
     case 'i' :
     case 'I' :
     {
       
       std::string fname = "animate_" + Convert<int,std::string>(image_number) + ".bmp";
       QImage qim = grabFrameBuffer();
       qim.save(QString(fname.c_str()),"BMP");
       ++image_number;
      
       
       break;
     }
       
     default :
       
       QWidget::keyPressEvent(event);
   }
}


extern void UpdateTimeSlider(int plot_ts);


void Tetrahedron::draw()
{    
    if (bFirst)
    {
       bFirst = false;
       return;
    }

    if (bRightButtonMenu) return;
    if (Data::bBusy) return;
    Data::bBusy = true;
    
//    cout << "start draw : " << draw_counter << "\n";
    
    Data::Update_vWindowVariable();
    Data::ReadVariables();
    
    for (int i=0; i<Data::num_plots; ++i)
    {
      if (Data::vPlotWindowInfo[i].bPlot)
      {
       // Data::mfix.GetVariableAtTimestep(Data::vPlotWindowInfo[i].win_var,Data::cur_ts);
	
	if (Data::vPlotWindowInfo[i].bIJ_plot)
	{  
	   Data::Plot_IJ_slice(i);
	   
	   if (Data::vPlotWindowInfo[i].bPlotGasVectors)
	   {
	      int index_vel = Data::FindVariableIndex("U_g");
	      
	      if (index_vel >=0)
	      {
		Data::mfix.GetVariableAtTimestep(index_vel,Data::cur_ts);	   
		std::vector<float> uvel = Data::mfix.var;
	   
		Data::mfix.GetVariableAtTimestep(index_vel+1,Data::cur_ts);	   
		std::vector<float> vvel = Data::mfix.var;
	   
		Data::Draw_IJ_vectors(i ,uvel , vvel);
	      }
	   }
	   
	   if (Data::vPlotWindowInfo[i].bPlotSolidVectors)
	   {
	      int index_vel = Data::FindVariableIndex("U_s_1");
	      
	      if (index_vel >=0)
	      {
		Data::mfix.GetVariableAtTimestep(index_vel,Data::cur_ts);	   
		std::vector<float> uvel = Data::mfix.var;
	   
		Data::mfix.GetVariableAtTimestep(index_vel+1,Data::cur_ts);	   
		std::vector<float> vvel = Data::mfix.var;
	   
		Data::Draw_IJ_vectors(i ,uvel , vvel);
	      }
	   }
	}
	
	Data::DrawTitle(i);
      }
      Data::DrawTime();
    }
      
    UpdateTimeSlider( Data::cur_ts );

    if (!Data::bPaused)  ++Data::cur_ts;
 //   Data::cur_ts += Data::step_code;
    
    if (Data::cur_ts >= Data::mfix.Get_NTIMES()) Data::cur_ts = 0;
    if (Data::cur_ts <  0                      ) Data::cur_ts = Data::mfix.Get_NTIMES() - 1;
    
    if (Data::bPaused)
      Data::step_code = 0;
    else
      Data::step_code = 1;
    
    ++draw_counter;

    Data::bBusy = false;
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

void SendKey(int key)
{
  if (tetra_window != 0)
  {
     unsigned short count = 1;
     
     QKeyEvent *event = new QKeyEvent ( QEvent::KeyPress, key, Qt::NoModifier ,QString(""),false,count);
     QCoreApplication::postEvent (tetra_window, event);
  }
}
   

