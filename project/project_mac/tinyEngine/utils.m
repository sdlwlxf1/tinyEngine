//
//  utils.m
//  tinyEngine
//
//  Created by tuyoo on 2017/3/9.
//  Copyright © 2017年 sdlwlxf. All rights reserved.
//

#import <Foundation/Foundation.h>

const char* getFilePath(const char* name, const char* type) {
    NSString* filePath = [[NSBundle mainBundle] pathForResource:[NSString stringWithUTF8String:name]
                                                         ofType:[NSString stringWithUTF8String:type]];
    return [filePath UTF8String];
}
