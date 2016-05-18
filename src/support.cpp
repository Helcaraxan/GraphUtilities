// GraphUtilities library - Directed graph manipulation and querying
// Copyright (C) 2016 - Duco van Amstel
//
// For more license information please see the 'LICENSE' file at the top-level
// of the source code tree.

#include <iostream>

#include "unistd.h"
#include "sys/ioctl.h"

#include "graph-utilities/implementation/support.hpp"

using namespace std;


// Variables for progress bars

int batchFlag = 0;
int progressBarFinish = 1;

static int pastRatio = -1;
static int terminalWidth = 0;
static int progressBarTitleWidth = 7;
static string progressBarTitle = "";


// Progress bar functions

void configureProgressBar(string * title, int finish) {
  if (finish != 0) {
    pastRatio = -1;
    progressBarFinish = finish;
  }

  if (title) {
    progressBarTitle = *title;
    progressBarTitleWidth = title->length() + 7;
  }
}

void resultProgressBar(int progress) {
  int intRatio;
  int refreshModulo;
  double floatRatio;

  if (batchFlag == 1)
    return;

  if (terminalWidth == 0) {
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

    terminalWidth = w.ws_col;
  }

  refreshModulo = progressBarFinish / (terminalWidth - progressBarTitleWidth);
  if (refreshModulo == 0)
    refreshModulo = 1;

  if (progress == progressBarFinish) {
    cout << progressBarTitle << "100% ["
      << string(terminalWidth - progressBarTitleWidth, '=') << "]\r" << flush;
    progressBarTitle.clear();
    return;
  }

  if (progress % refreshModulo != 0)
    return;

  floatRatio = ((double) progress) / ((double) progressBarFinish);
  intRatio = floatRatio * (terminalWidth - progressBarTitleWidth);

  if (intRatio == pastRatio)
    return;
  else
    pastRatio = intRatio;

  cout << progressBarTitle;
  cout.width(3);
  cout << right << (int) (floatRatio * 100) << "% [";

  string completeString(intRatio, '=');
  cout << string(intRatio, '=')
    << string(terminalWidth - progressBarTitleWidth - intRatio, ' ') << "]\r"
    << flush;
}
