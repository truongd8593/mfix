/********************************************************************************
** Form generated from reading UI file 'left_dialog.ui'
**
** Created: Fri Apr 23 10:25:32 2010
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
#include <QtGui/QDialog>
#include <QtGui/QFrame>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_leftDialog
{
public:
    QGroupBox *groupBox;
    QWidget *layoutWidget;
    QVBoxLayout *verticalLayout;
    QRadioButton *zoom_xy;
    QRadioButton *zoom_x;
    QRadioButton *zoom_y;
    QPushButton *inButton;
    QPushButton *outButton;
    QPushButton *resetButton;
    QPushButton *upButton;
    QPushButton *downButton;
    QCheckBox *checkBox_xsym;
    QPushButton *button_createTable;
    QSpinBox *spin_itable;
    QSpinBox *spin_jtable;
    QSpinBox *spin_ktable;
    QFrame *line;
    QLabel *label;
    QLabel *label_2;
    QLabel *label_3;
    QProgressBar *progressBar_table;
    QLineEdit *table_fname;
    QPushButton *biginButton;
    QPushButton *bigupButton;
    QPushButton *bigdownButton;
    QPushButton *bigoutButton;

    void setupUi(QDialog *leftDialog)
    {
        if (leftDialog->objectName().isEmpty())
            leftDialog->setObjectName(QString::fromUtf8("leftDialog"));
        leftDialog->resize(194, 805);
        groupBox = new QGroupBox(leftDialog);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setGeometry(QRect(20, 20, 111, 281));
        groupBox->setStyleSheet(QString::fromUtf8("color: rgb(255, 0, 0);"));
        groupBox->setAlignment(Qt::AlignCenter);
        layoutWidget = new QWidget(groupBox);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(20, 40, 75, 77));
        verticalLayout = new QVBoxLayout(layoutWidget);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        zoom_xy = new QRadioButton(layoutWidget);
        zoom_xy->setObjectName(QString::fromUtf8("zoom_xy"));

        verticalLayout->addWidget(zoom_xy);

        zoom_x = new QRadioButton(layoutWidget);
        zoom_x->setObjectName(QString::fromUtf8("zoom_x"));

        verticalLayout->addWidget(zoom_x);

        zoom_y = new QRadioButton(layoutWidget);
        zoom_y->setObjectName(QString::fromUtf8("zoom_y"));

        verticalLayout->addWidget(zoom_y);

        inButton = new QPushButton(groupBox);
        inButton->setObjectName(QString::fromUtf8("inButton"));
        inButton->setGeometry(QRect(10, 140, 85, 27));
        outButton = new QPushButton(groupBox);
        outButton->setObjectName(QString::fromUtf8("outButton"));
        outButton->setGeometry(QRect(10, 180, 85, 27));
        resetButton = new QPushButton(groupBox);
        resetButton->setObjectName(QString::fromUtf8("resetButton"));
        resetButton->setGeometry(QRect(10, 220, 85, 27));
        upButton = new QPushButton(leftDialog);
        upButton->setObjectName(QString::fromUtf8("upButton"));
        upButton->setGeometry(QRect(30, 320, 85, 27));
        downButton = new QPushButton(leftDialog);
        downButton->setObjectName(QString::fromUtf8("downButton"));
        downButton->setGeometry(QRect(30, 370, 85, 27));
        checkBox_xsym = new QCheckBox(leftDialog);
        checkBox_xsym->setObjectName(QString::fromUtf8("checkBox_xsym"));
        checkBox_xsym->setGeometry(QRect(30, 420, 91, 61));
        button_createTable = new QPushButton(leftDialog);
        button_createTable->setObjectName(QString::fromUtf8("button_createTable"));
        button_createTable->setGeometry(QRect(40, 510, 85, 27));
        spin_itable = new QSpinBox(leftDialog);
        spin_itable->setObjectName(QString::fromUtf8("spin_itable"));
        spin_itable->setGeometry(QRect(70, 570, 46, 24));
        spin_jtable = new QSpinBox(leftDialog);
        spin_jtable->setObjectName(QString::fromUtf8("spin_jtable"));
        spin_jtable->setGeometry(QRect(70, 620, 46, 24));
        spin_ktable = new QSpinBox(leftDialog);
        spin_ktable->setObjectName(QString::fromUtf8("spin_ktable"));
        spin_ktable->setGeometry(QRect(70, 670, 46, 24));
        line = new QFrame(leftDialog);
        line->setObjectName(QString::fromUtf8("line"));
        line->setGeometry(QRect(10, 490, 118, 3));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);
        label = new QLabel(leftDialog);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(10, 570, 41, 20));
        label_2 = new QLabel(leftDialog);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(10, 620, 41, 20));
        label_3 = new QLabel(leftDialog);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(10, 670, 51, 20));
        progressBar_table = new QProgressBar(leftDialog);
        progressBar_table->setObjectName(QString::fromUtf8("progressBar_table"));
        progressBar_table->setGeometry(QRect(20, 720, 118, 23));
        progressBar_table->setValue(24);
        table_fname = new QLineEdit(leftDialog);
        table_fname->setObjectName(QString::fromUtf8("table_fname"));
        table_fname->setGeometry(QRect(20, 760, 113, 24));
        biginButton = new QPushButton(leftDialog);
        biginButton->setObjectName(QString::fromUtf8("biginButton"));
        biginButton->setGeometry(QRect(130, 160, 41, 27));
        bigupButton = new QPushButton(leftDialog);
        bigupButton->setObjectName(QString::fromUtf8("bigupButton"));
        bigupButton->setGeometry(QRect(130, 320, 41, 27));
        bigdownButton = new QPushButton(leftDialog);
        bigdownButton->setObjectName(QString::fromUtf8("bigdownButton"));
        bigdownButton->setGeometry(QRect(130, 370, 41, 27));
        bigoutButton = new QPushButton(leftDialog);
        bigoutButton->setObjectName(QString::fromUtf8("bigoutButton"));
        bigoutButton->setGeometry(QRect(130, 200, 41, 27));

        retranslateUi(leftDialog);

        QMetaObject::connectSlotsByName(leftDialog);
    } // setupUi

    void retranslateUi(QDialog *leftDialog)
    {
        leftDialog->setWindowTitle(QApplication::translate("leftDialog", "left dialog", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("leftDialog", "   zoom option", 0, QApplication::UnicodeUTF8));
        zoom_xy->setText(QApplication::translate("leftDialog", "X and Y", 0, QApplication::UnicodeUTF8));
        zoom_x->setText(QApplication::translate("leftDialog", "X only", 0, QApplication::UnicodeUTF8));
        zoom_y->setText(QApplication::translate("leftDialog", "Y only", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        inButton->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        inButton->setText(QApplication::translate("leftDialog", "zoom in", 0, QApplication::UnicodeUTF8));
        outButton->setText(QApplication::translate("leftDialog", "zoom out", 0, QApplication::UnicodeUTF8));
        resetButton->setText(QApplication::translate("leftDialog", "reset", 0, QApplication::UnicodeUTF8));
        upButton->setText(QApplication::translate("leftDialog", "shift up", 0, QApplication::UnicodeUTF8));
        downButton->setText(QApplication::translate("leftDialog", "shift down", 0, QApplication::UnicodeUTF8));
        checkBox_xsym->setText(QApplication::translate("leftDialog", "K slices - \n"
"X symetry", 0, QApplication::UnicodeUTF8));
        button_createTable->setText(QApplication::translate("leftDialog", "create table", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("leftDialog", "I-index", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("leftDialog", "J-index", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("leftDialog", "K-index", 0, QApplication::UnicodeUTF8));
        table_fname->setText(QApplication::translate("leftDialog", "table.csv", 0, QApplication::UnicodeUTF8));
        biginButton->setText(QApplication::translate("leftDialog", ">>", 0, QApplication::UnicodeUTF8));
        bigupButton->setText(QApplication::translate("leftDialog", ">>", 0, QApplication::UnicodeUTF8));
        bigdownButton->setText(QApplication::translate("leftDialog", ">>", 0, QApplication::UnicodeUTF8));
        bigoutButton->setText(QApplication::translate("leftDialog", ">>", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class leftDialog: public Ui_leftDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_LEFT_DIALOG_H
