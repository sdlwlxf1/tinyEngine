//
//  utils.h
//  tinyEngine
//
//  Created by lixuefeng on 2017/3/9.
//  Copyright © 2017年 sdlwlxf. All rights reserved.
//

#ifndef utils_h
#define utils_h

#include "tiny3D.h"

const char* getFilePath(const char* name, const char* type);

int load_obj(float bmin[3], float bmax[3], const char* name, const char* type);

int load_png_image( const char *name, const char *type, unsigned int **bits, unsigned int *width, unsigned int *height);

#endif /* utils_h */
