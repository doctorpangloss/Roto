/****************************************************************************
** Meta object code from reading C++ file 'centralWidget.h'
**
** Created: Sat Dec 15 04:02:07 2012
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "centralWidget.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'centralWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_CentralWidget[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      31,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       4,       // signalCount

 // signals: signature, parameters, type, tag, flags
      14,   39,   41,   41, 0x05,
      42,   41,   41,   41, 0x05,
      52,   41,   41,   41, 0x05,
      62,   41,   41,   41, 0x05,

 // slots: signature, parameters, type, tag, flags
      78,   39,   41,   41, 0x0a,
      95,  112,   41,   41, 0x0a,
     114,  112,   41,   41, 0x0a,
     136,  112,   41,   41, 0x0a,
     157,  112,   41,   41, 0x0a,
     180,   41,   41,   41, 0x0a,
     191,  112,   41,   41, 0x0a,
     206,  112,   41,   41, 0x0a,
     221,  112,   41,   41, 0x0a,
     236,  112,   41,   41, 0x0a,
     251,  112,   41,   41, 0x0a,
     266,  112,   41,   41, 0x0a,
     282,  112,   41,   41, 0x0a,
     297,  112,   41,   41, 0x0a,
     312,  112,   41,   41, 0x0a,
     328,  112,   41,   41, 0x0a,
     344,  112,   41,   41, 0x0a,
     360,   41,   41,   41, 0x0a,
     374,   41,   41,   41, 0x0a,
     391,  414,   41,   41, 0x0a,
     418,  414,   41,   41, 0x0a,
     440,   41,   41,   41, 0x0a,
     452,  414,   41,   41, 0x0a,
     469,   41,   41,   41, 0x0a,
     476,   41,   41,   41, 0x0a,
     484,   41,   41,   41, 0x0a,
     501,   41,  513,   41, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_CentralWidget[] = {
    "CentralWidget\0internalFrameChange(int)\0"
    "i\0\0rotoNow()\0drawNow()\0restoreCursor()\0"
    "frameChange(int)\0showChanged(int)\0s\0"
    "imageShowChanged(int)\0propModeChanged(int)\0"
    "renderModeChanged(int)\0calcLerp()\0"
    "cwToggle(bool)\0plToggle(bool)\0"
    "usToggle(bool)\0unToggle(bool)\0"
    "psToggle(bool)\0upsToggle(bool)\0"
    "dwToggle(bool)\0prToggle(bool)\0"
    "stpToggle(bool)\0fixToggle(bool)\0"
    "fhdToggle(bool)\0corrPropAll()\0"
    "rotoToDrawPath()\0imageAlphaChanged(int)\0"
    "val\0drawAlphaChanged(int)\0resetZoom()\0"
    "sensChanged(int)\0copy()\0paste()\0"
    "wireValChanged()\0zoomOutIn()\0int\0"
};

void CentralWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        CentralWidget *_t = static_cast<CentralWidget *>(_o);
        switch (_id) {
        case 0: _t->internalFrameChange((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->rotoNow(); break;
        case 2: _t->drawNow(); break;
        case 3: _t->restoreCursor(); break;
        case 4: _t->frameChange((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: _t->showChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: _t->imageShowChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: _t->propModeChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 8: _t->renderModeChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: _t->calcLerp(); break;
        case 10: _t->cwToggle((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 11: _t->plToggle((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 12: _t->usToggle((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 13: _t->unToggle((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 14: _t->psToggle((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 15: _t->upsToggle((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 16: _t->dwToggle((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 17: _t->prToggle((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 18: _t->stpToggle((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 19: _t->fixToggle((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 20: _t->fhdToggle((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 21: _t->corrPropAll(); break;
        case 22: _t->rotoToDrawPath(); break;
        case 23: _t->imageAlphaChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 24: _t->drawAlphaChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 25: _t->resetZoom(); break;
        case 26: _t->sensChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 27: _t->copy(); break;
        case 28: _t->paste(); break;
        case 29: _t->wireValChanged(); break;
        case 30: { int _r = _t->zoomOutIn();
            if (_a[0]) *reinterpret_cast< int*>(_a[0]) = _r; }  break;
        default: ;
        }
    }
}

const QMetaObjectExtraData CentralWidget::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject CentralWidget::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_CentralWidget,
      qt_meta_data_CentralWidget, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &CentralWidget::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *CentralWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *CentralWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_CentralWidget))
        return static_cast<void*>(const_cast< CentralWidget*>(this));
    return QGLWidget::qt_metacast(_clname);
}

int CentralWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 31)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 31;
    }
    return _id;
}

// SIGNAL 0
void CentralWidget::internalFrameChange(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void CentralWidget::rotoNow()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}

// SIGNAL 2
void CentralWidget::drawNow()
{
    QMetaObject::activate(this, &staticMetaObject, 2, 0);
}

// SIGNAL 3
void CentralWidget::restoreCursor()
{
    QMetaObject::activate(this, &staticMetaObject, 3, 0);
}
QT_END_MOC_NAMESPACE
