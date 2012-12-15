/****************************************************************************
** Meta object code from reading C++ file 'mymainwindow.h'
**
** Created: Sat Dec 15 04:02:06 2012
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "mymainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mymainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MyMainWindow[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      35,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      13,   28,   28,   28, 0x0a,
      29,   28,   28,   28, 0x0a,
      42,   28,   28,   28, 0x0a,
      52,   28,   28,   28, 0x0a,
      65,   82,   28,   28, 0x0a,
      88,   28,   28,   28, 0x0a,
     100,   28,   28,   28, 0x0a,
     117,   28,   28,   28, 0x0a,
     139,   28,   28,   28, 0x0a,
     159,   28,   28,   28, 0x0a,
     179,   28,   28,   28, 0x0a,
     205,   28,   28,   28, 0x0a,
     233,   28,   28,   28, 0x0a,
     259,   28,   28,   28, 0x0a,
     270,   28,   28,   28, 0x0a,
     283,   28,   28,   28, 0x0a,
     298,   28,   28,   28, 0x0a,
     310,   28,   28,   28, 0x0a,
     325,   28,   28,   28, 0x0a,
     337,   28,   28,   28, 0x0a,
     354,   28,   28,   28, 0x0a,
     376,   28,   28,   28, 0x0a,
     398,   28,   28,   28, 0x0a,
     420,   28,   28,   28, 0x0a,
     440,   28,   28,   28, 0x0a,
     461,   28,   28,   28, 0x0a,
     470,   28,   28,   28, 0x0a,
     480,   28,   28,   28, 0x0a,
     490,   28,   28,   28, 0x0a,
     506,  529,   28,   28, 0x0a,
     533,   28,   28,   28, 0x0a,
     549,   28,   28,   28, 0x0a,
     568,   28,   28,   28, 0x0a,
     580,   28,   28,   28, 0x0a,
     593,  610,   28,   28, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_MyMainWindow[] = {
    "MyMainWindow\0saveCurvesAs()\0\0saveCurves()\0"
    "saveXML()\0saveImages()\0drawChanged(int)\0"
    "state\0pickColor()\0pickPenTexture()\0"
    "pickCharcoalTexture()\0pickPastelTexture()\0"
    "pickPencilTexture()\0pickPencilStrokeTexture()\0"
    "pickCharcoalStrokeTexture()\0"
    "pickPastelStrokeTexture()\0loadCalc()\0"
    "writeTrack()\0improveTrack()\0redoTrack()\0"
    "trackForward()\0copyTrack()\0nlevelsChanged()\0"
    "smooth0DerivChanged()\0smooth2DerivChanged()\0"
    "smooth1DerivChanged()\0edgeWeightChanged()\0"
    "shape2DerivChanged()\0runCan()\0rotoNow()\0"
    "drawNow()\0deleteAllTime()\0"
    "aToolChanged(QAction*)\0act\0selectChanged()\0"
    "currColorChanged()\0setCursor()\0"
    "saveMattes()\0frameChange(int)\0i\0"
};

void MyMainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MyMainWindow *_t = static_cast<MyMainWindow *>(_o);
        switch (_id) {
        case 0: _t->saveCurvesAs(); break;
        case 1: _t->saveCurves(); break;
        case 2: _t->saveXML(); break;
        case 3: _t->saveImages(); break;
        case 4: _t->drawChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: _t->pickColor(); break;
        case 6: _t->pickPenTexture(); break;
        case 7: _t->pickCharcoalTexture(); break;
        case 8: _t->pickPastelTexture(); break;
        case 9: _t->pickPencilTexture(); break;
        case 10: _t->pickPencilStrokeTexture(); break;
        case 11: _t->pickCharcoalStrokeTexture(); break;
        case 12: _t->pickPastelStrokeTexture(); break;
        case 13: _t->loadCalc(); break;
        case 14: _t->writeTrack(); break;
        case 15: _t->improveTrack(); break;
        case 16: _t->redoTrack(); break;
        case 17: _t->trackForward(); break;
        case 18: _t->copyTrack(); break;
        case 19: _t->nlevelsChanged(); break;
        case 20: _t->smooth0DerivChanged(); break;
        case 21: _t->smooth2DerivChanged(); break;
        case 22: _t->smooth1DerivChanged(); break;
        case 23: _t->edgeWeightChanged(); break;
        case 24: _t->shape2DerivChanged(); break;
        case 25: _t->runCan(); break;
        case 26: _t->rotoNow(); break;
        case 27: _t->drawNow(); break;
        case 28: _t->deleteAllTime(); break;
        case 29: _t->aToolChanged((*reinterpret_cast< QAction*(*)>(_a[1]))); break;
        case 30: _t->selectChanged(); break;
        case 31: _t->currColorChanged(); break;
        case 32: _t->setCursor(); break;
        case 33: _t->saveMattes(); break;
        case 34: _t->frameChange((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData MyMainWindow::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject MyMainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MyMainWindow,
      qt_meta_data_MyMainWindow, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MyMainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MyMainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MyMainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyMainWindow))
        return static_cast<void*>(const_cast< MyMainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MyMainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 35)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 35;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
