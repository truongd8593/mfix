/********************************************************************************
** Form generated from reading UI file 'bottom_dialog.ui'
**
** Created: Thu Apr 29 11:01:58 2010
**      by: Qt User Interface Compiler version 4.6.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_BOTTOM_DIALOG_H
#define UI_BOTTOM_DIALOG_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSlider>

QT_BEGIN_NAMESPACE

class Ui_bottomDialog
{
public:
    QPushButton *pauseButton;
    QSlider *sliderTime;
    QLabel *timeLabel;
    QLabel *label;

    void setupUi(QDialog *bottomDialog)
    {
        if (bottomDialog->objectName().isEmpty())
            bottomDialog->setObjectName(QString::fromUtf8("bottomDialog"));
        bottomDialog->resize(602, 133);
        pauseButton = new QPushButton(bottomDialog);
        pauseButton->setObjectName(QString::fromUtf8("pauseButton"));
        pauseButton->setGeometry(QRect(20, 10, 85, 27));
        sliderTime = new QSlider(bottomDialog);
        sliderTime->setObjectName(QString::fromUtf8("sliderTime"));
        sliderTime->setGeometry(QRect(20, 80, 551, 21));
        sliderTime->setOrientation(Qt::Horizontal);
        timeLabel = new QLabel(bottomDialog);
        timeLabel->setObjectName(QString::fromUtf8("timeLabel"));
        timeLabel->setGeometry(QRect(260, 46, 81, 20));
        timeLabel->setFrameShape(QFrame::Box);
        timeLabel->setFrameShadow(QFrame::Sunken);
        timeLabel->setAlignment(Qt::AlignCenter);
        label = new QLabel(bottomDialog);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(270, 20, 55, 16));
        label->setAlignment(Qt::AlignCenter);

        retranslateUi(bottomDialog);

        QMetaObject::connectSlotsByName(bottomDialog);
    } // setupUi

    void retranslateUi(QDialog *bottomDialog)
    {
        bottomDialog->setWindowTitle(QApplication::translate("bottomDialog", "bottom dialog", 0, QApplication::UnicodeUTF8));
        pauseButton->setText(QApplication::translate("bottomDialog", "pause/go", 0, QApplication::UnicodeUTF8));
        timeLabel->setText(QString());
        label->setText(QApplication::translate("bottomDialog", "time", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class bottomDialog: public Ui_bottomDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_BOTTOM_DIALOG_H
