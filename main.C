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



//USEAGE: main fileRoot(#) startFrame(#) endFrame(#)

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <qapplication.h>
#include <qgl.h>
#include <qfiledialog.h>

#include "mymainwindow.h"
#include "imgSequence.h"

int main( int argc, char *argv[] ) {

  printf("hello world\n");
  QApplication a( argc, argv );
  printf("hello world 2\n");
  char **myArgv = a.argv();
  int myArgc = a.argc();
  ImgSequence* is;
  bool fromFile = false;
  FILE* fp;
  QFile qf;
  
  if ( !QGLFormat::hasOpenGL() ) {
    qWarning( "This system has no OpenGL support. Exiting." );
    return -1;
  }
  printf("hello world 3\n");

  if (myArgc==4)
    is = new ImgSequence(myArgv[1],atoi(myArgv[2]),atoi(myArgv[3]));
  else if (myArgc==2) {
    fromFile = true;
    is = new ImgSequence();
    fp = fopen(myArgv[1],"r");
    if (!fp) {
      printf("Can't find specified param file\n");
      exit(0);
    }
    QFileInfo fi(myArgv[1]);
    qf.setName(myArgv[1]);
    qf.open(IO_ReadOnly);
    QString s = fi.dirPath();
    s.append("/frame");
    is->loadInitialDataqt(&qf,s.ascii()); //HERE
    //is->loadInitialData(fp,s.ascii());
  }
  else if (myArgc==1) {
    QString s = QFileDialog::getExistingDirectory(NULL,
						  0, 0,  "Choose clip directory",
						  TRUE );
    if (!s.isEmpty()) {
      QString root = s;
      root.append("/frame");
      s.append(QString("frame.param"));
      fromFile = true;
      is = new ImgSequence();
      fp = fopen(s.ascii(),"r");
      qf.setName(myArgv[1]);
      qf.open(IO_ReadOnly);
      is->loadInitialDataqt(&qf,root.ascii()); // HERE
      //is->loadInitialData(fp, root.ascii());
    }
    else
      exit(0);
  }
  else {
    printf("useage: main fileRoot startFrame endFrame\n");
    printf("or: main file.param\n");
    printf("or: main\n");
    exit(0);
  }

  MyMainWindow* mainWin = new MyMainWindow(is);
  a.setMainWidget(mainWin);
  mainWin->setFocus();
  mainWin->show();
  //mainWin->extWin->undock();

  if (fromFile) { // HERE
    is->loadRestDataqt(&qf);
    //is->loadRestData(fp);
    //fclose(fp);
    qf.close();
  }
  
  
  a.connect( &a, SIGNAL( lastWindowClosed() ), &a, SLOT( quit() ) );
  int result = a.exec();
  delete mainWin;
  delete is;
  return result;
}
  
