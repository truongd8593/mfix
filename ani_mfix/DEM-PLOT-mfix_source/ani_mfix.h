#ifndef ANI_MFIX_H
#define ANI_MFIX_H

#include <QMainWindow>

class QSplitter;
class QDialog;
class QScrollArea;
class Tetrahedron;
class leftDialog;

class AniMfix : public QMainWindow
{
    Q_OBJECT

public:
   AniMfix();
   ~AniMfix();

protected:
//    void showEvent(QShowEvent *event);

    

private:
    void readSettings();
    void writeSettings();

    QSplitter *mainSplitter;
    QSplitter *rightSplitter;
    leftDialog * dlg_l;
    QDialog * dlg_b;
    QDialog * dlg_gl;
    Tetrahedron * tetra;

    QScrollArea * scrollArea_l;
    QScrollArea * scrollArea_b;
};

#endif
