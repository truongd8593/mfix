#include <QtGui>
#include <QProgressBar>


#include <iostream>
#include <vector>

#include "left_dialog.h"

#include "data.h"

extern void MyTogglePause();
extern void StepGL();

leftDialog::leftDialog(QWidget * parent) : QDialog(parent)
{
   setupUi(this);
   connect(inButton             , SIGNAL(clicked())                    , this  , SLOT(IN_clicked()));
   connect(outButton            , SIGNAL(clicked())                    , this  , SLOT(OUT_clicked()));
   connect(resetButton          , SIGNAL(clicked())                    , this  , SLOT(RESET_clicked()));
   connect(upButton             , SIGNAL(clicked())                    , this  , SLOT(UP_clicked()));
   connect(spin_plus            , SIGNAL(clicked())                    , this  , SLOT(SPIN_PLUS_clicked()));
   connect(spin_minus           , SIGNAL(clicked())                    , this  , SLOT(SPIN_MINUS_clicked()));
   connect(tilt_plus            , SIGNAL(clicked())                    , this  , SLOT(TILT_PLUS_clicked()));
   connect(tilt_minus           , SIGNAL(clicked())                    , this  , SLOT(TILT_MINUS_clicked()));
   connect(downButton           , SIGNAL(clicked())                    , this  , SLOT(DOWN_clicked()));
   connect(checkBox_3d          , SIGNAL(clicked())                    , this  , SLOT(C3D_clicked()));
   connect(checkBox_trace       , SIGNAL(clicked())                    , this  , SLOT(TRACE_clicked()));
   connect(pushButton_image     , SIGNAL(clicked())                    , this  , SLOT(IMAGE_clicked()));
   connect(comboBox_color       , SIGNAL(currentIndexChanged(int))     , this  , SLOT(COLOR_METHOD_changed(int)));
   connect(comboBox_sphere_res  , SIGNAL(currentIndexChanged(int))     , this  , SLOT(SPHERE_RES_changed(int)));
   
   comboBox_sphere_res->setCurrentIndex(2);
 }

leftDialog::~leftDialog()
{

}


void leftDialog::RESET_clicked()
{
  Data::scale  = 1.0f;
  Data::xtrans = 0.0f;
  Data::xrot   = 0.0f;
  Data::yrot   = 0.0f;
  Data::xtrans = 0.0f;
  Data::ytrans = 0.0f;
  StepGL();
}

void leftDialog::C3D_clicked()
{
    Data::b3D_bubbles = !Data::b3D_bubbles;
    StepGL();
}

void leftDialog::TRACE_clicked()
{
    Data::bTrace = !Data::bTrace;
    StepGL();
}

void leftDialog::IN_clicked()
{
  Data::scale *= 1.2f;
  StepGL();
}


void leftDialog::OUT_clicked()
{
  Data::scale *= 0.8f;
  StepGL();
}



void leftDialog::UP_clicked()
{
 //  Data::ytrans += (Data::ymax - Data::ymin) * Data::scale * 0.1f;
   Data::ytrans -= 0.05;
   StepGL();
}

void leftDialog::DOWN_clicked()
{
 //   Data::ytrans -= (Data::ymax - Data::ymin) * Data::scale * 0.1f;
    Data::ytrans += 0.05;
    StepGL();
}


void leftDialog::COLOR_METHOD_changed(int method)
{
   Data::bBusy = true;
   
   Data::color_method = method;

   Data::bBusy = false;
   StepGL();
}

void leftDialog::SPHERE_RES_changed(int res)
{
   Data::bBusy = true;
   
   if (res == 0) Data::sphere_resolution =  3;
   if (res == 1) Data::sphere_resolution =  5;
   if (res == 2) Data::sphere_resolution =  7;
   if (res == 3) Data::sphere_resolution = 10;
   if (res == 4) Data::sphere_resolution = 15;

   Data::bBusy = false;
   StepGL();
}

void leftDialog::SPIN_PLUS_clicked()
{
  Data::yrot -= 5;
  StepGL();
}

void leftDialog::SPIN_MINUS_clicked()
{
  Data::yrot += 5;
  StepGL();
}

void leftDialog::TILT_PLUS_clicked()
{
  Data::xrot -= 5;
  StepGL();
}

void leftDialog::TILT_MINUS_clicked()
{
  Data::xrot += 5;
  StepGL();
}

void leftDialog::IMAGE_clicked()
{
  Data::bImage = true;
  StepGL();
}
