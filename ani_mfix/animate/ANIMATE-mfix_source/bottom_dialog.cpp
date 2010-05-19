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
   connect(comboBox_variables , SIGNAL(currentIndexChanged(int))     , this  , SLOT(Variable_changed(int)));
   connect(comboBox_jstep     , SIGNAL(currentIndexChanged(int))     , this  , SLOT(JSTEP_changed(int)));
   connect(comboBox_kslice    , SIGNAL(currentIndexChanged(int))     , this  , SLOT(KSLICE_changed(int)));
   connect(vmin               , SIGNAL(textChanged(const QString &)) , this  , SLOT(VMIN_changed(const QString &)));
   connect(vmax               , SIGNAL(textChanged(const QString &)) , this  , SLOT(VMAX_changed(const QString &)));
   connect(sliderTime         , SIGNAL(sliderReleased())             , this  , SLOT(TIME_changed()));
   connect(sliderTime         , SIGNAL(sliderMoved(int))             , this  , SLOT(TIME_moved(int)));
   connect(sliderTime         , SIGNAL(valueChanged(int))            , this  , SLOT(TIME_value(int)));
   connect(sliderTime         , SIGNAL(sliderPressed())              , this  , SLOT(TIME_pressed()));
   connect(spinNPLOTS         , SIGNAL(valueChanged(int))            , this  , SLOT(NPLOTS_changed(int)));
    
   if (Data::mfix.kmax2 > 1)
   {     
     for (int k=2; k<=Data::mfix.kmax2-1; ++k)
     {
	comboBox_kslice->addItem( Convert<int,std::string>(k).c_str() );
     }
   }
   else
   {
     label_K->hide();
     comboBox_kslice->hide();
   }
   
   if (Data::mfix.jmax2 < 500)
   {
     label_jstep->hide();
     comboBox_jstep->hide();
   }
   
   
   
   pauseButton->setText("go");
   
   for (int i=0; i<Data::mfix.Get_NVARS(); ++i) 
   {
     comboBox_variables->addItem( Data::mfix.variable_names[i].c_str() );
   }
   

   
   sliderTime->setRange(0,Data::mfix.Get_NTIMES()-1);
   
   timeLabel->setText( Convert<float,std::string>( Data::mfix.times[0] ).c_str() );
   
   timeSliderPressed = false;

   spinNPLOTS->setRange(1,9);
   
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

void bottomDialog::Variable_changed(int index)
{
   Data::bBusy = true;
   
   Data::vPlotWindowInfo[Data::cur_win].win_var = index;
   
   vmin->setText( Convert<float,string>( Data::vVariableRange[index].vmin ).c_str() );
   vmax->setText( Convert<float,string>( Data::vVariableRange[index].vmax ).c_str() );
   
   Data::bBusy = false;
   StepGL();
}

void bottomDialog::JSTEP_changed(int index)
{
   Data::bBusy = true;
   
   if (index == 0) Data::jstep =  1;
   if (index == 1) Data::jstep =  2;
   if (index == 2) Data::jstep =  4;
   if (index == 3) Data::jstep =  8;
   if (index == 4) Data::jstep = 12;
   if (index == 5) Data::jstep = 16;
   
   
   Data::bBusy = false;
   StepGL();
}
void bottomDialog::VMIN_changed(const QString & s)
{
   int var =  Data::vPlotWindowInfo[Data::cur_win].win_var;
   
   Data::vVariableRange[var].vmin = s.toFloat();
   
   StepGL();
}

void bottomDialog::VMAX_changed(const QString & s)
{
   int var =  Data::vPlotWindowInfo[Data::cur_win].win_var;
   
   Data::vVariableRange[var].vmax = s.toFloat();
   
   StepGL();
}

void bottomDialog::TIME_changed()
{
  timeSliderPressed = false;
  Data::cur_ts = sliderTime->value();
  StepGL();
}

void bottomDialog::TIME_moved(int ts)
{
  timeLabel->setText( Convert<float,std::string>( Data::mfix.times[ts] ).c_str() );
}

void bottomDialog::TIME_value(int ts)
{
  if (!timeSliderPressed) 
  {
    Data::cur_ts = ts;
    StepGL();
    timeLabel->setText( Convert<float,std::string>( Data::mfix.times[ts] ).c_str() );
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


void CurrentWindowChanged(int cw);

void bottomDialog::NPLOTS_changed(int nplots)
{
  Data::bBusy = true;
  
//   cout << "NPLOTS_changed : start\n";
   if (nplots > (int)Data::vPlotWindowInfo.size())
   {
      Data::vPlotWindowInfo.resize(nplots);
   }

   if (Data::cur_win >= nplots) 
   {
     CurrentWindowChanged(0);
   }
   else
   {
     CurrentWindowChanged(nplots-1);
   }
//   cout << "NPLOTS_changed : bbb\n";

   Data::num_plots = nplots;
   Data::bBusy = false;
   StepGL();
//   cout << "NPLOTS_changed : end\n";
}

void CurrentWindowChanged(int cw)
{
    Data::cur_win = cw;
    
    if (dlg != 0)
    {
       dlg->comboBox_variables->setCurrentIndex(Data::vPlotWindowInfo[cw].win_var);
       dlg->comboBox_kslice->setCurrentIndex( Data::vPlotWindowInfo[Data::cur_win].win_k - 1 );
    }
}

void RightClick_variable_change(int var)
{
  Data::vPlotWindowInfo[Data::cur_win].win_var = var;
  CurrentWindowChanged( Data::cur_win );
}


void bottomDialog::KSLICE_changed(int k)
{  
  Data::vPlotWindowInfo[Data::cur_win].win_k = k + 1;
  StepGL();
}



