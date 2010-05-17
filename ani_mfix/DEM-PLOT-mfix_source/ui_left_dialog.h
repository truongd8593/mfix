/********************************************************************************
** Form generated from reading UI file 'left_dialog.ui'
**
** Created: Mon May 17 10:02:59 2010
**      by: Qt User Interface Compiler version 4.6.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_LEFT_DIALOG_H
#define UI_LEFT_DIALOG_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>

QT_BEGIN_NAMESPACE

class Ui_leftDialog
{
public:
    QGroupBox *groupBox;
    QPushButton *inButton;
    QPushButton *outButton;
    QPushButton *resetButton;
    QPushButton *upButton;
    QPushButton *downButton;
    QCheckBox *checkBox_3d;
    QCheckBox *checkBox_trace;
    QComboBox *comboBox_color;
    QLabel *label;
    QPushButton *spin_plus;
    QPushButton *spin_minus;
    QPushButton *tilt_plus;
    QPushButton *tilt_minus;
    QLabel *label_2;
    QLabel *label_3;
    QPushButton *pushButton_image;
    QComboBox *comboBox_sphere_res;

    void setupUi(QDialog *leftDialog)
    {
        if (leftDialog->objectName().isEmpty())
            leftDialog->setObjectName(QString::fromUtf8("leftDialog"));
        leftDialog->resize(157, 704);
        groupBox = new QGroupBox(leftDialog);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setGeometry(QRect(20, 20, 111, 161));
        groupBox->setStyleSheet(QString::fromUtf8("color: rgb(255, 0, 0);"));
        groupBox->setAlignment(Qt::AlignCenter);
        inButton = new QPushButton(groupBox);
        inButton->setObjectName(QString::fromUtf8("inButton"));
        inButton->setGeometry(QRect(10, 30, 85, 27));
        outButton = new QPushButton(groupBox);
        outButton->setObjectName(QString::fromUtf8("outButton"));
        outButton->setGeometry(QRect(10, 70, 85, 27));
        resetButton = new QPushButton(groupBox);
        resetButton->setObjectName(QString::fromUtf8("resetButton"));
        resetButton->setGeometry(QRect(10, 110, 85, 27));
        upButton = new QPushButton(leftDialog);
        upButton->setObjectName(QString::fromUtf8("upButton"));
        upButton->setGeometry(QRect(30, 180, 85, 27));
        downButton = new QPushButton(leftDialog);
        downButton->setObjectName(QString::fromUtf8("downButton"));
        downButton->setGeometry(QRect(30, 220, 85, 27));
        checkBox_3d = new QCheckBox(leftDialog);
        checkBox_3d->setObjectName(QString::fromUtf8("checkBox_3d"));
        checkBox_3d->setGeometry(QRect(30, 390, 91, 21));
        checkBox_3d->setChecked(true);
        checkBox_trace = new QCheckBox(leftDialog);
        checkBox_trace->setObjectName(QString::fromUtf8("checkBox_trace"));
        checkBox_trace->setGeometry(QRect(20, 480, 101, 41));
        comboBox_color = new QComboBox(leftDialog);
        comboBox_color->setObjectName(QString::fromUtf8("comboBox_color"));
        comboBox_color->setGeometry(QRect(10, 580, 121, 28));
        label = new QLabel(leftDialog);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(40, 550, 55, 16));
        spin_plus = new QPushButton(leftDialog);
        spin_plus->setObjectName(QString::fromUtf8("spin_plus"));
        spin_plus->setGeometry(QRect(60, 270, 31, 27));
        spin_minus = new QPushButton(leftDialog);
        spin_minus->setObjectName(QString::fromUtf8("spin_minus"));
        spin_minus->setGeometry(QRect(100, 270, 31, 27));
        tilt_plus = new QPushButton(leftDialog);
        tilt_plus->setObjectName(QString::fromUtf8("tilt_plus"));
        tilt_plus->setGeometry(QRect(60, 320, 31, 27));
        tilt_minus = new QPushButton(leftDialog);
        tilt_minus->setObjectName(QString::fromUtf8("tilt_minus"));
        tilt_minus->setGeometry(QRect(100, 320, 31, 27));
        label_2 = new QLabel(leftDialog);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(20, 270, 31, 16));
        label_3 = new QLabel(leftDialog);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(20, 320, 21, 16));
        pushButton_image = new QPushButton(leftDialog);
        pushButton_image->setObjectName(QString::fromUtf8("pushButton_image"));
        pushButton_image->setGeometry(QRect(20, 630, 111, 27));
        comboBox_sphere_res = new QComboBox(leftDialog);
        comboBox_sphere_res->setObjectName(QString::fromUtf8("comboBox_sphere_res"));
        comboBox_sphere_res->setGeometry(QRect(30, 420, 91, 28));

        retranslateUi(leftDialog);

        QMetaObject::connectSlotsByName(leftDialog);
    } // setupUi

    void retranslateUi(QDialog *leftDialog)
    {
        leftDialog->setWindowTitle(QApplication::translate("leftDialog", "left dialog", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("leftDialog", "   zoom option", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        inButton->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        inButton->setText(QApplication::translate("leftDialog", "zoom in", 0, QApplication::UnicodeUTF8));
        outButton->setText(QApplication::translate("leftDialog", "zoom out", 0, QApplication::UnicodeUTF8));
        resetButton->setText(QApplication::translate("leftDialog", "reset", 0, QApplication::UnicodeUTF8));
        upButton->setText(QApplication::translate("leftDialog", "shift up", 0, QApplication::UnicodeUTF8));
        downButton->setText(QApplication::translate("leftDialog", "shift down", 0, QApplication::UnicodeUTF8));
        checkBox_3d->setText(QApplication::translate("leftDialog", "3D particles", 0, QApplication::UnicodeUTF8));
        checkBox_trace->setText(QApplication::translate("leftDialog", "    trace\n"
"first particle", 0, QApplication::UnicodeUTF8));
        comboBox_color->clear();
        comboBox_color->insertItems(0, QStringList()
         << QApplication::translate("leftDialog", "constant", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("leftDialog", "initial position : L/R", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("leftDialog", "diameter", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("leftDialog", "velocity magnitude", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("leftDialog", "V velocity : +/-", 0, QApplication::UnicodeUTF8)
        );
        label->setText(QApplication::translate("leftDialog", "color by", 0, QApplication::UnicodeUTF8));
        spin_plus->setText(QApplication::translate("leftDialog", "+", 0, QApplication::UnicodeUTF8));
        spin_minus->setText(QApplication::translate("leftDialog", "-", 0, QApplication::UnicodeUTF8));
        tilt_plus->setText(QApplication::translate("leftDialog", "+", 0, QApplication::UnicodeUTF8));
        tilt_minus->setText(QApplication::translate("leftDialog", "-", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("leftDialog", "spin", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("leftDialog", "tilt", 0, QApplication::UnicodeUTF8));
        pushButton_image->setText(QApplication::translate("leftDialog", "create image file", 0, QApplication::UnicodeUTF8));
        comboBox_sphere_res->clear();
        comboBox_sphere_res->insertItems(0, QStringList()
         << QApplication::translate("leftDialog", "low", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("leftDialog", "fair", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("leftDialog", "average", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("leftDialog", "good", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("leftDialog", "high", 0, QApplication::UnicodeUTF8)
        );
    } // retranslateUi

};

namespace Ui {
    class leftDialog: public Ui_leftDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_LEFT_DIALOG_H
