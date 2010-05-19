/********************************************************************************
** Form generated from reading UI file 'bottom_dialog.ui'
**
** Created: Mon Apr 26 09:00:51 2010
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
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QSlider>
#include <QtGui/QSpinBox>

QT_BEGIN_NAMESPACE

class Ui_bottomDialog
{
public:
    QPushButton *pauseButton;
    QComboBox *comboBox_variables;
    QLineEdit *vmin;
    QLineEdit *vmax;
    QSlider *sliderTime;
    QLabel *timeLabel;
    QLabel *label;
    QSpinBox *spinNPLOTS;
    QLabel *label_2;
    QLabel *label_K;
    QComboBox *comboBox_jstep;
    QLabel *label_jstep;
    QComboBox *comboBox_kslice;

    void setupUi(QDialog *bottomDialog)
    {
        if (bottomDialog->objectName().isEmpty())
            bottomDialog->setObjectName(QString::fromUtf8("bottomDialog"));
        bottomDialog->resize(719, 199);
        pauseButton = new QPushButton(bottomDialog);
        pauseButton->setObjectName(QString::fromUtf8("pauseButton"));
        pauseButton->setGeometry(QRect(20, 10, 85, 27));
        comboBox_variables = new QComboBox(bottomDialog);
        comboBox_variables->setObjectName(QString::fromUtf8("comboBox_variables"));
        comboBox_variables->setGeometry(QRect(140, 10, 161, 28));
        vmin = new QLineEdit(bottomDialog);
        vmin->setObjectName(QString::fromUtf8("vmin"));
        vmin->setGeometry(QRect(190, 50, 113, 24));
        vmax = new QLineEdit(bottomDialog);
        vmax->setObjectName(QString::fromUtf8("vmax"));
        vmax->setGeometry(QRect(190, 90, 113, 24));
        sliderTime = new QSlider(bottomDialog);
        sliderTime->setObjectName(QString::fromUtf8("sliderTime"));
        sliderTime->setGeometry(QRect(130, 150, 551, 21));
        sliderTime->setOrientation(Qt::Horizontal);
        timeLabel = new QLabel(bottomDialog);
        timeLabel->setObjectName(QString::fromUtf8("timeLabel"));
        timeLabel->setGeometry(QRect(20, 149, 81, 20));
        timeLabel->setFrameShape(QFrame::Box);
        timeLabel->setFrameShadow(QFrame::Sunken);
        timeLabel->setAlignment(Qt::AlignCenter);
        label = new QLabel(bottomDialog);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(30, 123, 55, 16));
        label->setAlignment(Qt::AlignCenter);
        spinNPLOTS = new QSpinBox(bottomDialog);
        spinNPLOTS->setObjectName(QString::fromUtf8("spinNPLOTS"));
        spinNPLOTS->setGeometry(QRect(40, 70, 41, 31));
        QSizePolicy sizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(spinNPLOTS->sizePolicy().hasHeightForWidth());
        spinNPLOTS->setSizePolicy(sizePolicy);
        label_2 = new QLabel(bottomDialog);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(40, 50, 61, 16));
        label_K = new QLabel(bottomDialog);
        label_K->setObjectName(QString::fromUtf8("label_K"));
        label_K->setGeometry(QRect(350, 20, 41, 16));
        comboBox_jstep = new QComboBox(bottomDialog);
        comboBox_jstep->setObjectName(QString::fromUtf8("comboBox_jstep"));
        comboBox_jstep->setGeometry(QRect(520, 110, 80, 28));
        label_jstep = new QLabel(bottomDialog);
        label_jstep->setObjectName(QString::fromUtf8("label_jstep"));
        label_jstep->setGeometry(QRect(540, 80, 41, 16));
        comboBox_kslice = new QComboBox(bottomDialog);
        comboBox_kslice->setObjectName(QString::fromUtf8("comboBox_kslice"));
        comboBox_kslice->setGeometry(QRect(340, 50, 80, 28));

        retranslateUi(bottomDialog);

        QMetaObject::connectSlotsByName(bottomDialog);
    } // setupUi

    void retranslateUi(QDialog *bottomDialog)
    {
        bottomDialog->setWindowTitle(QApplication::translate("bottomDialog", "bottom dialog", 0, QApplication::UnicodeUTF8));
        pauseButton->setText(QApplication::translate("bottomDialog", "pause/go", 0, QApplication::UnicodeUTF8));
        timeLabel->setText(QString());
        label->setText(QApplication::translate("bottomDialog", "time", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("bottomDialog", "# of Plots", 0, QApplication::UnicodeUTF8));
        label_K->setText(QApplication::translate("bottomDialog", "K-Slice", 0, QApplication::UnicodeUTF8));
        comboBox_jstep->clear();
        comboBox_jstep->insertItems(0, QStringList()
         << QApplication::translate("bottomDialog", "1", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("bottomDialog", "2", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("bottomDialog", "4", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("bottomDialog", "8", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("bottomDialog", "12", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("bottomDialog", "16", 0, QApplication::UnicodeUTF8)
        );
        label_jstep->setText(QApplication::translate("bottomDialog", "jstep", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class bottomDialog: public Ui_bottomDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_BOTTOM_DIALOG_H
