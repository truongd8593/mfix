#include <QtGui>

#include <iostream>

#include "ani_mfix.h"

#include "left_dialog.h"
#include "bottom_dialog.h"
#include "tetrahedron.h"

AniMfix::AniMfix()
{
   dlg_l  = new leftDialog;
   dlg_b  = new bottomDialog;
   tetra  = new Tetrahedron;

   scrollArea_l = new QScrollArea;
   scrollArea_l->setWidget(dlg_l);

   scrollArea_b = new QScrollArea;
   scrollArea_b->setWidget(dlg_b);

   rightSplitter = new QSplitter(Qt::Vertical);
   rightSplitter->addWidget(tetra);
   rightSplitter->addWidget(scrollArea_b);
   rightSplitter->setStretchFactor(1,1);

   mainSplitter = new QSplitter(Qt::Horizontal);
   mainSplitter->addWidget(scrollArea_l);
   mainSplitter->addWidget(rightSplitter);
   mainSplitter->setStretchFactor(1,1);

   setCentralWidget(mainSplitter);

   setWindowTitle(tr("animate mfix"));

   readSettings();

   tetra->setFocus();
}


AniMfix::~AniMfix()
{
   writeSettings();

   delete dlg_l;
   delete dlg_b;
   delete tetra;
   delete  scrollArea_b;
   delete rightSplitter;
   delete mainSplitter;
}


void AniMfix::readSettings()
{
    QSettings settings("MFIX", "animate mfix");

    settings.beginGroup("mainWindow");
    restoreGeometry(settings.value("geometry").toByteArray());
    mainSplitter->restoreState(
            settings.value("mainSplitter").toByteArray());
    rightSplitter->restoreState(
            settings.value("rightSplitter").toByteArray());
    settings.endGroup();
}

void AniMfix::writeSettings()
{
    QSettings settings("MFIX", "animate mfix");

    settings.beginGroup("mainWindow");
    settings.setValue("geometry", saveGeometry());
    settings.setValue("mainSplitter", mainSplitter->saveState());
    settings.setValue("rightSplitter", rightSplitter->saveState());
    settings.endGroup();
}

