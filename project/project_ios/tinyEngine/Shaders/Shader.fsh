//
//  Shader.fsh
//  tinyEngine
//
//  Created by 李雪峰 on 2017/3/5.
//  Copyright © 2017年 李雪峰. All rights reserved.
//

varying lowp vec4 colorVarying;

void main()
{
    gl_FragColor = colorVarying;
}
