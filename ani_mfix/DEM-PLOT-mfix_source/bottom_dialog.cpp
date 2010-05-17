#include <QtGui>

#include <iostream>
#include <string>

#include "bottom_dialog.h"

#include "data.h"

extern void MyTogglePause();
extern void StepGL();

namespace
{
  bottomDialog * dlg = 0;
}

bottomDialog::bottomDialog(QWidget * parent) : QDialog(parent)
{
   setupUi(this);
   
   
   connect(pauseButton        , SIGNAL(clicked())                    , this  , SLOT(PPP()));
   connect(sliderTime         , SIGNAL(sliderReleased())             , this  , SLOT(TIME_changed()));
   connect(sliderTime         , SIGNAL(sliderMoved(int))             , this  , SLOT(TIME_moved(int)));
   connect(sliderTime         , SIGNAL(valueChanged(int))            , this  , SLOT(TIME_value(int)));
   connect(sliderTime         , SIGNAL(sliderPressed())              , this  , SLOT(TIME_pressed()));
    

  
   pauseButton->setText("go");
    
   sliderTime->setRange(0,Data::vTimes.size()-1);
   
   timeLabel->setText( Convert<double,std::string>( Data::vTimes[0] ).c_str() );
     
   timeSliderPressed = false;
   
   dlg = this;
}

void bottomDialog::PPP()
{
  Data::bPaused = !Data::bPaused;
  
  if (Data::bPaused)
    pauseButton->setText("go");
  else
    pauseButton->setText("pause");
  
  MyTogglePause();
}

void bottomDialog::Variable_changed(int /* index */)
{

}

void bottomDialog::JSTEP_changed(int /* index */)
{

}
void bottomDialog::VMIN_changed(const QString & /* s */)
{

}

void bottomDialog::VMAX_changed(const QString & /* s */)
{

}

void bottomDialog::TIME_changed()
{
  timeSliderPressed = false;
  Data::cur_ts = sliderTime->value();
  Data::ReadNext_XML_file();
    Data::bAdd_to_vpos = true;
  StepGL();
}

void bottomDialog::TIME_moved(int ts)
{
  timeLabel->setText( Convert<double,std::string>( Data::vTimes[ts] ).c_str() );
}

void bottomDialog::TIME_value(int ts)
{
  if (!timeSliderPressed) 
  {
    Data::bAdd_to_vpos = true;
    Data::cur_ts = ts;
    Data::ReadNext_XML_file();
    StepGL();
    timeLabel->setText( Convert<float,std::string>( Data::vTimes[ts] ).c_str() );
  }
}

void bottomDialog::TIME_pressed()
{
  timeSliderPressed = true;
}

void UpdateTimeSlider(int plot_ts)
{
  if (dlg != 0) dlg->sliderTime->setValue(plot_ts);
}


void bottomDialog::NPLOTS_changed(int /* nplots */)
{

}

void RightClick_variable_change(int /* var */)
{

}


void bottomDialog::KSLICE_changed(int /* k */)
{  
 
}



