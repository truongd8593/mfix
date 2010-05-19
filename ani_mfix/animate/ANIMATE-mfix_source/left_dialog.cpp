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
   connect(zoom_xy            , SIGNAL(clicked())                    , this  , SLOT(XY_clicked()));
   connect(zoom_x             , SIGNAL(clicked())                    , this  , SLOT(X_clicked()));
   connect(zoom_y             , SIGNAL(clicked())                    , this  , SLOT(Y_clicked()));
   connect(inButton           , SIGNAL(clicked())                    , this  , SLOT(IN_clicked()));
   connect(biginButton        , SIGNAL(clicked())                    , this  , SLOT(BIGIN_clicked()));
   connect(bigoutButton       , SIGNAL(clicked())                    , this  , SLOT(BIGOUT_clicked()));
   connect(bigupButton        , SIGNAL(clicked())                    , this  , SLOT(BIGUP_clicked()));
   connect(bigdownButton      , SIGNAL(clicked())                    , this  , SLOT(BIGDOWN_clicked()));
   connect(outButton          , SIGNAL(clicked())                    , this  , SLOT(OUT_clicked()));
   connect(resetButton        , SIGNAL(clicked())                    , this  , SLOT(RESET_clicked()));
   connect(upButton           , SIGNAL(clicked())                    , this  , SLOT(UP_clicked()));
   connect(downButton         , SIGNAL(clicked())                    , this  , SLOT(DOWN_clicked()));
   connect(checkBox_xsym      , SIGNAL(clicked())                    , this  , SLOT(XSYM_clicked()));
   connect(button_createTable , SIGNAL(clicked())                    , this  , SLOT(CREATE_TABLE_clicked()));
 
   zoom_xy->setChecked(true);
   
   spin_itable->setRange(1,Data::mfix.imax2);
   spin_jtable->setRange(1,Data::mfix.jmax2);
   spin_ktable->setRange(1,Data::mfix.kmax2);
   
   spin_itable->setValue(1);
   spin_jtable->setValue(1);
   spin_ktable->setValue(1);
   
   progressBar_table->reset();
   progressBar_table->setMinimum(0);
   progressBar_table->setMaximum(Data::mfix.times.size());
   progressBar_table->setValue(0);
   
   qtable = 0;
}

leftDialog::~leftDialog()
{
   for (size_t i=0; i<vQI.size(); ++i)
   {
     delete vQI[i];
   }
   delete qtable;
}


void leftDialog::XY_clicked()
{
  Data::zoom_state = 0;
}

void leftDialog::X_clicked()
{
  Data::zoom_state = 1;
}

void leftDialog::Y_clicked()
{
  Data::zoom_state = -1;
}

//////////////////////////////////////////////
                                            //
extern void SendKey(int key);

void leftDialog::RESET_clicked()
{
  SendKey(Qt::Key_R);
}

void leftDialog::IN_clicked()
{
  Data::zoom_factor = 1;
  SendKey(Qt::Key_Z);
}

void leftDialog::BIGIN_clicked()
{
  Data::zoom_factor = 2;
  SendKey(Qt::Key_Z);
}

void leftDialog::OUT_clicked()
{
  Data::zoom_factor = 1;
  SendKey(Qt::Key_X);
}

void leftDialog::BIGOUT_clicked()
{
  Data::zoom_factor = 2;
  SendKey(Qt::Key_X);
}
                                            //
//////////////////////////////////////////////


void leftDialog::UP_clicked()
{
  Data::yshift_k -= 0.05f;
  StepGL();
}

void leftDialog::DOWN_clicked()
{
  Data::yshift_k += 0.05f;
  StepGL();
}

void leftDialog::BIGUP_clicked()
{
  Data::yshift_k -= 0.25f;
  StepGL();
}

void leftDialog::BIGDOWN_clicked()
{
  Data::yshift_k += 0.25f;
  StepGL();
}

void leftDialog::XSYM_clicked()
{
  bool bChecked = checkBox_xsym->isChecked();
  
  for (size_t i=0; i<Data::vPlotWindowInfo.size(); ++i)
  {
     Data::vPlotWindowInfo[i].bFlipX = bChecked;
  }
  
  StepGL();
}

void leftDialog::CREATE_TABLE_clicked()
{ 
   // delete previously created QTableWidgetItem * elements
   
   for (size_t i=0; i<vQI.size(); ++i)
   {
     delete vQI[i];
   }
   vQI.clear();
   
   //
  
//   delete qtable;  
   if (qtable ==0) qtable = new QTableWidget;
   
   qtable->setColumnCount(Data::num_plots+1);
   qtable->setRowCount(Data::mfix.times.size()+1);

   progressBar_table->reset();
   progressBar_table->setMinimum(0);
   progressBar_table->setMaximum(Data::mfix.times.size());
   progressBar_table->setValue(0);
   
   
   MyTable table;
   
   table.values.resize( Data::num_plots );
   
   Data::vTable.clear();
   
   int itable = spin_itable->value() - 1;
   int jtable = spin_jtable->value() - 1;
   int ktable = spin_ktable->value() - 1;
   
   int ijk = ktable * Data::mfix.imax2 * Data::mfix.jmax2 + jtable*Data::mfix.imax2 + itable;
   
   int previous_variable = -1;
   
   for (int ts=0; ts<Data::mfix.times.size(); ++ts)
   {
     table.time = Data::mfix.times[ts];
     
     for (int win=0; win<Data::num_plots; ++win)
     {
       int current_variable = Data::vPlotWindowInfo[win].win_var;
       
       if (current_variable != previous_variable)
       {
          Data::mfix.GetVariableAtTimestep(current_variable,ts);
          previous_variable = -1; //current_variable; (time change )
       }

       table.values[win] = Data::mfix.var[ijk];
     }
     Data::vTable.push_back(table);
     progressBar_table->setValue(ts+1);     
   }
   
   QString qs = table_fname->text();
   
   std::string fname = qs.toStdString();
   
   ofstream out(fname.c_str());
   
   int row = 0;
   out << "time , ";
   QTableWidgetItem * qi = new QTableWidgetItem("time");
   qtable->setItem(0,0,qi);
   vQI.push_back(qi);
   for (int win=0; win<Data::num_plots; ++win)
   {
       QTableWidgetItem * qi = new QTableWidgetItem(Data::mfix.variable_names[  Data::vPlotWindowInfo[win].win_var ].c_str());
       qtable->setItem(row,win+1,qi);
       vQI.push_back(qi);

       out << Data::mfix.variable_names[  Data::vPlotWindowInfo[win].win_var ];
       if (win != Data::num_plots-1) out << " , ";
   }
   out << "\n";
   
   for (size_t i=0; i<Data::vTable.size(); ++i)
   {
     ++row;
     out << Data::vTable[i].time << " , ";

     QTableWidgetItem * qi = new QTableWidgetItem(Convert<float,string>(Data::vTable[i].time).c_str());
     qtable->setItem(row,0,qi);
     vQI.push_back(qi);
     
     for (size_t j=0; j<Data::vTable[i].values.size(); ++j)
     {
       out << Data::vTable[i].values[j];
       QTableWidgetItem * qi = new QTableWidgetItem(Convert<float,string>(Data::vTable[i].values[j]).c_str());
       qtable->setItem(row,j+1,qi);
       vQI.push_back(qi);
       if (j != Data::vTable[i].values.size()-1)
	 out << " , ";
     }
     out << "\n";
   }
     
   Data::vTable.clear();
   qtable->show();
}


