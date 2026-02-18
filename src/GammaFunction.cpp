//
//  GammaFunction.cpp
//  biclust
//
//  Created by Thierry Chekouo Tekougang on 2021-02-20.
//  Copyright © 2021 Thierry Chekouo Tekougang. All rights reserved.
//

#include "GammaFunction.hpp"
#include <cmath>
extern "C" float LGammaF(float);
  float LGammaF(float x)
  {
      return (float)std::lgamma(x);
      // ...
  }
