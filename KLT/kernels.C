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



#include "kernels.h"

Kernels::Kernels(const float sigma) {
  const float factor = 0.01;   /* for truncating tail */
  int i;
  
  assert(MAX_KERNEL_WIDTH % 2 == 1);
  assert(sigma >= 0.0);
  _sigma = sigma;

  /* Compute kernels, and automatically determine widths */
  {
    const int hw = MAX_KERNEL_WIDTH / 2;
    float max_gauss = 1.0f, max_gaussderiv = sigma*exp(-0.5f);
	
    /* Compute gauss and deriv */
    for (i = -hw ; i <= hw ; i++)  {
      _gauss.data[i+hw]      = exp(-i*i / (2*sigma*sigma));
      _gauss_deriv.data[i+hw] = -i * _gauss.data[i+hw];
    }

    /* Compute widths */
    _gauss.width = MAX_KERNEL_WIDTH;
    for (i = -hw ; fabs(_gauss.data[i+hw] / max_gauss) < factor ; 
         i++, _gauss.width -= 2);
    _gauss_deriv.width = MAX_KERNEL_WIDTH;
    for (i = -hw ; fabs(_gauss_deriv.data[i+hw] / max_gaussderiv) < factor ; 
         i++, _gauss_deriv.width -= 2);
    if (_gauss.width == MAX_KERNEL_WIDTH || 
        _gauss_deriv.width == MAX_KERNEL_WIDTH)
      KLTError("(_computeKernels) MAX_KERNEL_WIDTH %d is too small for "
               "a sigma of %f", MAX_KERNEL_WIDTH, sigma);
  }

  /* Shift if width less than MAX_KERNEL_WIDTH */
  for (i = 0 ; i < _gauss.width ; i++)
    _gauss.data[i] = _gauss.data[i+(MAX_KERNEL_WIDTH-_gauss.width)/2];
  for (i = 0 ; i < _gauss_deriv.width ; i++)
    _gauss_deriv.data[i] = _gauss_deriv.data[i+(MAX_KERNEL_WIDTH-_gauss_deriv.width)/2];
  /* Normalize gauss and deriv */
  {
    const int hw = _gauss_deriv.width / 2;
    float den;
			
    den = 0.0;
    for (i = 0 ; i < _gauss.width ; i++)  den += _gauss.data[i];
    for (i = 0 ; i < _gauss.width ; i++)  _gauss.data[i] /= den;
    den = 0.0;
    for (i = -hw ; i <= hw ; i++)  den -= i*_gauss_deriv.data[i+hw];
    for (i = -hw ; i <= hw ; i++)  _gauss_deriv.data[i+hw] /= den;
  }

  
}
