## 大规模重构后，shader效果编写更加灵活，完整加载obj模型，渲染管线更接近标准。只差法线贴图没有加了！
## 经过两天debug，终于发现了问题，原来是计算重心坐标没有考虑透视插值，这个问题解决后，phong光照的问题也解决了。效果比以前更加自然：）
# tinyEngine帮助文档
#### 工程目录：

- **external**：依赖库文件
- **project**：多平台项目工程以及跟平台相关代码
- **resources**：核心代码和资源

#### 依赖库：

- [**SDL2**](https://www.libsdl.org/download-2.0.php)

跨平台的多媒体开发函数库，非常好用，cocos2dx也使用到了它哦

版本是SDL2-2.0.5，我已经交叉编译成了静态库，可以直接使用

- [**libpng**](http://www.libpng.org/pub/png/libpng.html)

跨平台的png格式读取处理库，如果静态库有问题，可以下载源代码[下载地址](https://sourceforge.net/projects/libpng/files/libpng16/1.6.28/)，自行编译

#### 多平台支持：

目前支持mac、ios、windows三个平台，android就暂时不弄了，过程太繁琐。

（windows的release好像有问题）

### 下面是重点！！

想要通过该工程学习光栅化渲染的同学，请务必备好以下两本书籍：

- **《3D数学基础：图形与游戏开发》**（建议认真看完）
- **《3D游戏编程大师技巧》上下册**

还有其他一些其他资料可以辅助理解

《[计算机图形学].(美国)Peter.Shirley》

[learnopengl-cn](http://learnopengl-cn.readthedocs.io/zh/latest/)一个学习opengl的网站，讲解光影部分挺好

#### 代码简介

##### tiny3D.h 

核心的数学算法和光栅化渲染实现，掌握了它可以说你就彻底弄明白了光栅化渲染

数学部分，全部是重点，特别是矩阵操作等，透视矩阵和旋转矩阵，是重中之重。

光栅化部分，为了实现phong光照加入了重心坐标计算，可能还有更好的实现方法，有能力的同学可以重构之。

##### main.c

非核心代码，包括sdl和纹理处理，有能力的同学可以自己实现并移植到任何平台

#### 关于学习方法

因为我也是通过参考大神韦易笑的代码来学习图形学的。这里我可以提供一些我的学习经验

首先要端正自己的学习态度，不要认为软渲染有多难，看整个实现也不过1000行左右代码，当然也不要被这微小的代码量所迷惑，量少就意味每一行代码每一个函数都潜藏着巨大的知识量。如果有不明白的，一定要查阅上面提到的两本书，比如[透视](http://blog.csdn.net/xiaowei_cqu/article/details/26471527)和旋转的算法公式，最好都要弄明白推导过程。然后自己实现一遍代码，通过实践加深对原理的理解。

**数学部分在于理解，光栅化部分在于实现**

请认真实现自己的光栅化部分，每个人的实现方法都不同，也许你的实现会更高效

个人认为软渲染是学习图形学的一个关键点。

等你也实现了自己的软渲染引擎，再看unity3D和shader，就会有一种醍醐灌顶的感觉，甚至可以自行脑补其实现方法呢



下面是截图：

![](https://github.com/sdlwlxf1/sdlwlxf1.github.io/blob/master/images/tinyEngine/screenshot%202017-03-12%20下午11.15.23.png)

![](https://github.com/sdlwlxf1/sdlwlxf1.github.io/blob/master/images/tinyEngine/screenshot%202017-03-12%20下午11.16.16.png)

![](https://github.com/sdlwlxf1/sdlwlxf1.github.io/blob/master/images/tinyEngine/screenshot%202017-03-12%20下午11.16.43.png)
