#ifndef BOTTOM_DIALOG_H
#define BOTTOM_DIALOG_H

#include "ui_bottom_dialog.h"

class bottomDialog : public QDialog , public Ui::bottomDialog
{
    Q_OBJECT

public:

  bottomDialog(QWidget * parent = 0);
  
private:
    
  bool timeSliderPressed;
  
public slots:
  
  void Variable_changed(int);
    

private slots:

  void PPP();
  void VMIN_changed(const QString & s);
  void VMAX_changed(const QString & s);
  void TIME_changed();
  void TIME_moved(int ts);
  void TIME_value(int ts);
  void TIME_pressed();
  void NPLOTS_changed(int nplots);
  void JSTEP_changed(int);
  void KSLICE_changed(int);
  
};

#endif
