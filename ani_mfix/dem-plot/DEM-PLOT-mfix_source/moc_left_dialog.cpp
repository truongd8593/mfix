/****************************************************************************
** Meta object code from reading C++ file 'left_dialog.h'
**
** Created: Mon May 17 10:04:12 2010
**      by: The Qt Meta Object Compiler version 62 (Qt 4.6.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "left_dialog.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'left_dialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.6.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_leftDialog[] = {

 // content:
       4,       // revision
       0,       // classname
       0,    0, // classinfo
      14,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x08,
      25,   11,   11,   11, 0x08,
      39,   11,   11,   11, 0x08,
      55,   11,   11,   11, 0x08,
      68,   11,   11,   11, 0x08,
      88,   11,   11,   11, 0x08,
     109,   11,   11,   11, 0x08,
     129,   11,   11,   11, 0x08,
     150,   11,   11,   11, 0x08,
     165,   11,   11,   11, 0x08,
     179,   11,   11,   11, 0x08,
     195,   11,   11,   11, 0x08,
     217,  211,   11,   11, 0x08,
     243,  211,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_leftDialog[] = {
    "leftDialog\0\0IN_clicked()\0OUT_clicked()\0"
    "RESET_clicked()\0UP_clicked()\0"
    "SPIN_PLUS_clicked()\0SPIN_MINUS_clicked()\0"
    "TILT_PLUS_clicked()\0TILT_MINUS_clicked()\0"
    "DOWN_clicked()\0C3D_clicked()\0"
    "IMAGE_clicked()\0TRACE_clicked()\0index\0"
    "COLOR_METHOD_changed(int)\0"
    "SPHERE_RES_changed(int)\0"
};

const QMetaObject leftDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_leftDialog,
      qt_meta_data_leftDialog, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &leftDialog::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *leftDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *leftDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_leftDialog))
        return static_cast<void*>(const_cast< leftDialog*>(this));
    if (!strcmp(_clname, "Ui::leftDialog"))
        return static_cast< Ui::leftDialog*>(const_cast< leftDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int leftDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: IN_clicked(); break;
        case 1: OUT_clicked(); break;
        case 2: RESET_clicked(); break;
        case 3: UP_clicked(); break;
        case 4: SPIN_PLUS_clicked(); break;
        case 5: SPIN_MINUS_clicked(); break;
        case 6: TILT_PLUS_clicked(); break;
        case 7: TILT_MINUS_clicked(); break;
        case 8: DOWN_clicked(); break;
        case 9: C3D_clicked(); break;
        case 10: IMAGE_clicked(); break;
        case 11: TRACE_clicked(); break;
        case 12: COLOR_METHOD_changed((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 13: SPHERE_RES_changed((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 14;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
