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



#ifndef RANGEDIALOG_H
#define RANGEDIALOG_H

#include <qdialog.h>
#include <qlineedit.h>
#include <qlayout.h>
#include <qvalidator.h>
#include <qpushbutton.h>
#include "jl_vectors.h"

class RangeDialog : public QDialog {

  Q_OBJECT

 public:

  RangeDialog(QWidget *parent, int lowVal, int highVal);
    
  Vec2i getRange() const;

  ~RangeDialog() {
    delete lowEdit; delete highEdit; delete lay;
  }

 private:
  QLineEdit *lowEdit, *highEdit;
  QGridLayout* lay;
  QPushButton *_ok, *_cancel;
};

#endif
