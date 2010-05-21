#ifndef LEFT_DIALOG_H
#define LEFT_DIALOG_H

#include "ui_left_dialog.h"

class QTableWidgetItem;
class QTableWidget;

class leftDialog : public QDialog , public Ui::leftDialog
{
    Q_OBJECT

public:

   leftDialog(QWidget * parent = 0);
   ~leftDialog();
   
private:
    
   std::vector<QTableWidgetItem*> vQI;
   QTableWidget * qtable;
    

private slots:
  
  
  void IN_clicked();
  void OUT_clicked();
  void RESET_clicked();
  
  void UP_clicked();
  void SPIN_PLUS_clicked();
  void SPIN_MINUS_clicked();
  void TILT_PLUS_clicked();
  void TILT_MINUS_clicked();
  void DOWN_clicked();
  void C3D_clicked();

  void IMAGE_clicked();

  void TRACE_clicked();
  void COLOR_METHOD_changed(int index);

  void SPHERE_RES_changed(int index);
  
};

#endif
