#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include <QGLWidget>

class Tetrahedron : public QGLWidget
{
    Q_OBJECT

public:
    Tetrahedron(QWidget *parent = 0);

    void togglePaused();
    void DrawTime();
    void DrawOutline();
    void Set_Viewport_time(int & ix, int & iy);

protected:
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();
    void mousePressEvent(QMouseEvent *event);
    void timerEvent(QTimerEvent * event);
    void keyPressEvent(QKeyEvent * event);
    void contextMenuEvent(QContextMenuEvent * event);

private:
    void draw();

    int  timerID;
};

#endif
