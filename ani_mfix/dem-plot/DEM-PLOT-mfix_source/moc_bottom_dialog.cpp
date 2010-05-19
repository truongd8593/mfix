/****************************************************************************
** Meta object code from reading C++ file 'bottom_dialog.h'
**
** Created: Mon May 17 08:29:26 2010
**      by: The Qt Meta Object Compiler version 62 (Qt 4.6.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "bottom_dialog.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'bottom_dialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.6.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_bottomDialog[] = {

 // content:
       4,       // revision
       0,       // classname
       0,    0, // classinfo
      11,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      14,   13,   13,   13, 0x0a,
      36,   13,   13,   13, 0x08,
      44,   42,   13,   13, 0x08,
      66,   42,   13,   13, 0x08,
      88,   13,   13,   13, 0x08,
     106,  103,   13,   13, 0x08,
     122,  103,   13,   13, 0x08,
     138,   13,   13,   13, 0x08,
     160,  153,   13,   13, 0x08,
     180,   13,   13,   13, 0x08,
     199,   13,   13,   13, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_bottomDialog[] = {
    "bottomDialog\0\0Variable_changed(int)\0"
    "PPP()\0s\0VMIN_changed(QString)\0"
    "VMAX_changed(QString)\0TIME_changed()\0"
    "ts\0TIME_moved(int)\0TIME_value(int)\0"
    "TIME_pressed()\0nplots\0NPLOTS_changed(int)\0"
    "JSTEP_changed(int)\0KSLICE_changed(int)\0"
};

const QMetaObject bottomDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_bottomDialog,
      qt_meta_data_bottomDialog, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &bottomDialog::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *bottomDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *bottomDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_bottomDialog))
        return static_cast<void*>(const_cast< bottomDialog*>(this));
    if (!strcmp(_clname, "Ui::bottomDialog"))
        return static_cast< Ui::bottomDialog*>(const_cast< bottomDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int bottomDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: Variable_changed((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: PPP(); break;
        case 2: VMIN_changed((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 3: VMAX_changed((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 4: TIME_changed(); break;
        case 5: TIME_moved((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: TIME_value((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: TIME_pressed(); break;
        case 8: NPLOTS_changed((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: JSTEP_changed((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: KSLICE_changed((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 11;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
