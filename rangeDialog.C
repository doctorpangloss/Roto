/* 

Copyright (C) 2004, Aseem Agarwala, roto@agarwala.org

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

*/



#include "rangeDialog.h"

RangeDialog::RangeDialog(QWidget *parent, int lowVal, int highVal) : QDialog (parent,0,TRUE) {
  //QIntValidator v(lowVal, highVal, this);
  QString num;
  lowEdit = new QLineEdit(this);
  highEdit = new QLineEdit(this);
  lowEdit->setMaxLength(3);
  //lowEdit->setValidator(&v);
  num.setNum(lowVal);
  lowEdit->setText(num);
  highEdit->setMaxLength(3);
  //highEdit->setValidator(&v);
  num.setNum(highVal);
  highEdit->setText(num);
  
  lay = new QGridLayout(this,2,2);
  layout()->add(lowEdit);
  layout()->add(highEdit);

  _ok = new QPushButton("Ok",this);
  _ok->setDefault(true);
  connect(_ok,SIGNAL(clicked()), this, SLOT(accept()));
  _cancel = new QPushButton("Cancel",this);
  connect(_cancel,SIGNAL(clicked()), this, SLOT(reject()));
  layout()->add(_ok);
  layout()->add(_cancel);
}

Vec2i RangeDialog::getRange() const {
  bool ok;
  Vec2i res;
  QString out = lowEdit->text();
  res.set_x(out.toInt(&ok));
  assert(ok);
  out = highEdit->text();
  res.set_y(out.toInt(&ok));
  assert(ok);
  //assert(res.x() >= lowVal && res.y() <= highVal);
  return res;
}
