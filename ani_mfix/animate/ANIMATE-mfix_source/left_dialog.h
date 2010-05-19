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
  
  void XY_clicked();
  void X_clicked();
  void Y_clicked();
  
  void XSYM_clicked();
  
  void IN_clicked();
  void BIGIN_clicked();
  void BIGOUT_clicked();
  void BIGUP_clicked();
  void BIGDOWN_clicked();
  void OUT_clicked();
  void RESET_clicked();
  
  void UP_clicked();
  void DOWN_clicked();
  
  void CREATE_TABLE_clicked();
};

#endif
