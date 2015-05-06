// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Version
//  1.01	Claes Rolén		2014-11-30	First version
//  1.02	Claes Rolén		2015-01-11	Added 'angle','specgram','fd_filter','gplot'
//                                      Changed file structure
//  1.03	Claes Rolén		2015-04-26	Added 'parser' class, 'err_handler','wrn_handler'
//                                      'freqz','phasez'


#ifndef ARMA_INCLUDES
#include <armadillo>
#endif

#ifndef SIGPACK_H
#define SIGPACK_H
#include "base/base.h"
#include "window/window.h"
#include "filter/filter.h"
#include "resampling/resampling.h"
#include "spectrum/spectrum.h"
#include "timing/timing.h"
#include "gplot/gplot.h"
#include "parser/parser.h"
#endif