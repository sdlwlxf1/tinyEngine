/*This source code copyrighted by Lazy Foo' Productions (2004-2015)
and may not be redistributed without written permission.*/

//Using SDL, SDL_image, standard IO, math, and strings
#include <stdbool.h>
#include "SDL.h"
#include "png.h"
#include <stdio.h>
#include "tiny3D.h"
#include "utils.h"


#define PI 3.141592653
#define angle_to_radian(X) ((X)/180*PI)
#define radian_to_angle(X) ((X)/PI*180)

//Screen dimension constants
const int SCREEN_WIDTH = 600;
const int SCREEN_HEIGHT = 600;

const int REAL_WIDTH = 600;
const int REAL_HEIGHT = 600;


//The window we'll be rendering to
SDL_Window* gWindow = NULL;

//The window renderer
SDL_Renderer* gRenderer = NULL;

bool sdl_init(int width, int height, const char *title)
{
	bool success = true;
	if( SDL_Init( SDL_INIT_VIDEO ) < 0 )
	{
		printf( "SDL could not initialize! SDL Error: %s\n", SDL_GetError() );
		success = false;
	}
	else
	{
		//Set texture filtering to linear
		if( !SDL_SetHint( SDL_HINT_RENDER_SCALE_QUALITY, "1" ) )
		{
			printf( "Warning: Linear texture filtering not enabled!" );
		}
		//Create window
		gWindow = SDL_CreateWindow( title, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_SHOWN );
		if( gWindow == NULL )
		{
			printf( "Window could not be created! SDL Error: %s\n", SDL_GetError() );
			success = false;
		}
		else
		{
			//Create renderer for window
			gRenderer = SDL_CreateRenderer( gWindow, -1, SDL_RENDERER_ACCELERATED );
			if( gRenderer == NULL )
			{
				printf( "Renderer could not be created! SDL Error: %s\n", SDL_GetError() );
				success = false;
			}
			else
			{
				//Initialize renderer color
				SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
			}
		}
	}

	return success;
}

void sdl_close()
{
	//Destroy window	
	SDL_DestroyRenderer( gRenderer );
	SDL_DestroyWindow( gWindow );
	gWindow = NULL;
	gRenderer = NULL;

	//Quit SDL subsystems
	//IMG_Quit();
	SDL_Quit();
}


vertex_t ground_mesh[6] = {
    // Positions                  // Texture Coords  //color           //rhw // Normals
    {{-0.5f,  0.0f, -0.5f, 1.0f}, {0.0f,  5.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.0f,  0.5f, 1.0f},  {0.0f,  0.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}},
    {{0.5f,  0.0f,  0.5f, 1.0f},  {5.0f,  0.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 ,  { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.0f,  0.5f, 1.0f},  {5.0f,  0.0f}, { 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.0f, -0.5f, 1.0f},  {5.0f,  5.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 ,  { 0.0f,1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.0f, -0.5f, 1.0f},  {0.0f,  5.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}}
};

vertex_t box_mesh[36] = {
    // Positions                  // Texture Coords  //color           //rhw // Normals
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  0.0f},{ 1.0f, 0.2f, 0.2f, 1.0f }, 1 , { 0.0f,  0.0f,-1.0f,0.0f}},
    {{-0.5f,  0.5f, -0.5f, 1.0f},{ 0.0f,  1.0f},{ 1.0f, 0.2f, 0.2f, 1.0f }, 1 , { 0.0f,  0.0f,-1.0f,0.0f}},
    {{0.5f,  0.5f, -0.5f, 1.0f}, {1.0f,  1.0f}, { 1.0f, 0.2f, 0.2f, 1.0f }, 1 , {0.0f,  0.0f,-1.0f ,0.0f}},
    {{0.5f,  0.5f, -0.5f, 1.0f}, { 1.0f,  1.0f}, { 1.0f, 0.2f, 0.2f, 1.0f }, 1 , {0.0f,  0.0f,-1.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f, 1.0f}, {1.0f,  0.0f}, { 1.0f, 0.2f, 0.2f, 1.0f }, 1 , {0.0f,  0.0f,-1.0f ,0.0f}},
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  0.0f},{ 1.0f, 0.2f, 0.2f, 1.0f }, 1 , { 0.0f,  0.0f,-1.0f,0.0f}},
    
    {{-0.5f, -0.5f,  0.5f, 1.0f},{ 0.0f,  0.0f},{ 0.2f, 1.0f, 0.2f, 1.0f }, 1 , { 0.0f,  0.0f, 1.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f, 1.0f},{ 1.0f,  0.0f}, { 0.2f, 1.0f, 0.2f, 1.0f }, 1 , {0.0f,  0.0f,  1.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},{ 1.0f,  1.0f}, { 0.2f, 1.0f, 0.2f, 1.0f }, 1 , {0.0f,  0.0f,  1.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},{ 1.0f,  1.0f}, { 0.2f, 1.0f, 0.2f, 1.0f }, 1 , {0.0f,  0.0f,  1.0f,0.0f}},
    {{-0.5f,  0.5f,  0.5f, 1.0f},{ 0.0f,  1.0f},{ 0.2f, 1.0f, 0.2f, 1.0f }, 1 , { 0.0f,  0.0f,  1.0f,0.0f}},
    {{-0.5f, -0.5f,  0.5f, 1.0f},{ 0.0f,  0.0f},{ 0.2f, 1.0f, 0.2f, 1.0f }, 1 , { 0.0f,  0.0f,  1.0f,0.0f}},
    
    {{-0.5f,  0.5f,  0.5f, 1.0f}, { 1.0f,  0.0f},{ 0.2f, 0.2f, 1.0f, 1.0f}, 1 , {-1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f, -0.5f, 1.0f},{ 1.0f,  1.0f},{ 0.2f, 0.2f, 1.0f, 1.0f }, 1 , { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  1.0f},{ 0.2f, 0.2f, 1.0f, 1.0f }, 1 , { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  1.0f},{ 0.2f, 0.2f, 1.0f, 1.0f }, 1 , { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f,  0.5f, 1.0f},{ 0.0f,  0.0f},{ 0.2f, 0.2f, 1.0f, 1.0f}, 1 , { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f,  0.5f,1.0f},{ 1.0f,  0.0f},{ 0.2f, 0.2f, 1.0f, 1.0f }, 1 , { -1.0f,  0.0f,  0.0f,0.0f}},
    
    {{0.5f,  0.5f,  0.5f,1.0f}, { 1.0f,  0.0f}, { 1.0f, 0.2f, 1.0f, 1.0f }, 1 , {1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f,1.0f},{ 0.0f,  0.0f}, { 1.0f, 0.2f, 1.0f, 1.0f }, 1 , { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f,1.0f},{ 0.0f,  1.0f}, { 1.0f, 0.2f, 1.0f, 1.0f }, 1 , { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f,1.0f},{ 0.0f,  1.0f}, { 1.0f, 0.2f, 1.0f, 1.0f }, 1 , { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f, -0.5f,1.0f},{ 1.0f,  1.0f}, { 1.0f, 0.2f, 1.0f, 1.0f }, 1 , { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f,1.0f},{ 1.0f,  0.0f}, { 1.0f, 0.2f, 1.0f, 1.0f }, 1 , { 1.0f,  0.0f,  0.0f,0.0f}},
    
    {{-0.5f, -0.5f, -0.5f,1.0f},{  0.0f,  1.0f},{ 1.0f, 1.0f, 0.2f, 1.0f }, 1 , {  0.0f, -1.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f,1.0f}, { 1.0f,  1.0f},{ 1.0f, 1.0f, 0.2f, 1.0f }, 1 ,  { 0.0f, -1.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f,1.0f}, { 1.0f,  0.0f}, { 1.0f, 1.0f, 0.2f, 1.0f }, 1 , { 0.0f, -1.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f,1.0f}, { 1.0f,  0.0f},{ 1.0f, 1.0f, 0.2f, 1.0f }, 1 ,  { 0.0f, -1.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f,  0.5f,1.0f},{ 0.0f,  0.0f},{ 1.0f, 1.0f, 0.2f, 1.0f }, 1 , {  0.0f, -1.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f, -0.5f,1.0f},{ 0.0f,  1.0f},{ 1.0f, 1.0f, 0.2f, 1.0f }, 1 , {  0.0f, -1.0f,  0.0f,0.0f}},

    {{-0.5f,  0.5f, -0.5f, 1.0f}, {0.0f,  1.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f,  0.5f, 1.0f},  {0.0f,  0.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},  {1.0f,  0.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 ,  { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},  {1.0f,  0.0f}, { 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f, -0.5f, 1.0f},  {1.0f,  1.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 ,  { 0.0f,1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f, -0.5f, 1.0f},  {0.0f,  1.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}}
};

typedef struct {
    vertex_t *mesh;
    int mesh_num;
    int material_id;
    int texture_id;
    bool shadow;
    
    bool dirty;
    point_t pos;
    vector_t scale;
    
    vector_t axis;
    float theta;
    
    matrix_t matrix;
} object_t;
#define MAX_NUM_OBJECT 100
object_t objects[MAX_NUM_OBJECT];
int object_count = 0;

typedef struct {
    IUINT32 **datas;           // 纹理数据
    IUINT32 datas_len;
    bool use_mipmap;        // 是否开启mipmap
    IUINT32 width;              // 纹理宽度
    IUINT32 height;             // 纹理高度
} texture_t;
#define MAX_NUM_TEXTURE 100
texture_t textures[MAX_NUM_TEXTURE];
int texture_count = 0;

int generate_mipmaps(texture_t *texture, float gamma) {
    IUINT32 **mipmaps = NULL;
    int num_mip_levels = logbase2ofx(texture->width) + 1;
    texture->datas_len = num_mip_levels;
    mipmaps = (IUINT32**)malloc(num_mip_levels * sizeof(IUINT32*));
    mipmaps[0] = texture->datas[0];
    int mip_width = texture->width;
    int mip_height = texture->height;
    for(int mip_level = 1; mip_level < num_mip_levels; mip_level++) {
        mip_width = mip_width >> 1;
        mip_height = mip_height >> 1;
        mipmaps[mip_level] = (IUINT32*)malloc(mip_width * mip_height * sizeof(IUINT32));
        IUINT32 *src_buffer = mipmaps[mip_level-1];
        IUINT32 *dest_buffer = mipmaps[mip_level];
        for(int x = 0; x < mip_width; x++)
        {
            for(int y = 0; y < mip_height; y++)
            {
                float r0, g0, b0, a0,
                r1, g1, b1, a1,
                r2, g2, b2, a2,
                r3, g3, b3, a3;
                int r_avg, g_avg, b_avg, a_avg;
                
                IUINT32 c = src_buffer[(x*2+0) + (y*2+0)*mip_width*2];
                b0 = c & 0xff;
                g0 = (c >> 8) & 0xff;
                r0 = (c >> 16) & 0xff;
                a0 = (c >> 24) & 0xff;
                
                c = src_buffer[(x*2+1) + (y*2+0)*mip_width*2];
                b1 = c & 0xff;
                g1 = (c >> 8) & 0xff;
                r1 = (c >> 16) & 0xff;
                a1 = (c >> 24) & 0xff;
                
                c = src_buffer[(x*2+0) + (y*2+1)*mip_width*2];
                b2 = c & 0xff;
                g2 = (c >> 8) & 0xff;
                r2 = (c >> 16) & 0xff;
                a2 = (c >> 24) & 0xff;
                
                c = src_buffer[(x*2+1) + (y*2+1)*mip_width*2];
                b3 = c & 0xff;
                g3 = (c >> 8) & 0xff;
                r3 = (c >> 16) & 0xff;
                a3 = (c >> 24) & 0xff;
                
                r_avg = (IUINT32)(0.5f + gamma*(r0+r1+r2+r3)/4);
                g_avg = (IUINT32)(0.5f + gamma*(g0+g1+g2+g3)/4);
                b_avg = (IUINT32)(0.5f + gamma*(b0+b1+b2+b3)/4);
                a_avg = (IUINT32)(0.5f + gamma*(b0+b1+b2+b3)/4);
                
                int R = CMID(r_avg, 0, 255);
                int G = CMID(g_avg, 0, 255);
                int B = CMID(b_avg, 0, 255);
                int A = CMID(a_avg, 0, 255);
                
                dest_buffer[x+y*mip_width] = (A << 24) | (R << 16) | (G << 8) | B;
            }
        }
    }
    free(texture->datas);
    texture->datas = mipmaps;
    return num_mip_levels;
}

#define PNG_BYTES_TO_CHECK 4
int load_png_image( const char *name, const char *type, texture_t *texture )
{
    FILE *fp;
    png_structp png_ptr;
    png_infop info_ptr;
    png_bytep* row_pointers;
    char buf[PNG_BYTES_TO_CHECK];
    int w, h, x, y, temp, color_type;
    
	char path[100];
    fp = fopen( getFilePath(name, type, path), "rb" );
    if( fp == NULL ) {
        return 1; /* 返回值 */
    }
    
    png_ptr = png_create_read_struct( PNG_LIBPNG_VER_STRING, 0, 0, 0 );
    info_ptr = png_create_info_struct( png_ptr );
    
    setjmp( png_jmpbuf(png_ptr) );
    /* 读取PNG_BYTES_TO_CHECK个字节的数据 */
    temp = (int)fread( buf, 1, PNG_BYTES_TO_CHECK, fp );
    /* 若读到的数据并没有PNG_BYTES_TO_CHECK个字节 */
    if( temp < PNG_BYTES_TO_CHECK ) {
        fclose(fp);
        png_destroy_read_struct( &png_ptr, &info_ptr, 0);
        return 2;/* 返回值 */
    }
    /* 检测数据是否为PNG的签名 */
    temp = png_sig_cmp( (png_bytep)buf, (png_size_t)0, PNG_BYTES_TO_CHECK );
    /* 如果不是PNG的签名，则说明该文件不是PNG文件 */
    if( temp != 0 ) {
        fclose(fp);
        png_destroy_read_struct( &png_ptr, &info_ptr, 0);
        return 3;/* 返回值 */
    }
    
    /* 复位文件指针 */
    rewind( fp );
    /* 开始读文件 */
    png_init_io( png_ptr, fp );
    /* 读取PNG图片信息 */
    png_read_png( png_ptr, info_ptr, PNG_TRANSFORM_EXPAND, 0 );
    /* 获取图像的色彩类型 */
    color_type = png_get_color_type( png_ptr, info_ptr );
    /* 获取图像的宽高 */
    w = png_get_image_width( png_ptr, info_ptr );
    h = png_get_image_height( png_ptr, info_ptr );

    IUINT32 *bits = (IUINT32*)malloc(sizeof(IUINT32) * w * h);
    
    /* 获取图像的所有行像素数据，row_pointers里边就是rgba数据 */
    row_pointers = png_get_rows( png_ptr, info_ptr );
    /* 根据不同的色彩类型进行相应处理 */
    switch( color_type ) {
        case PNG_COLOR_TYPE_RGB_ALPHA:
            for( y=0; y<h; ++y ) {
                for( x=0; x < w; ++x ) {
                    bits[y*w+x] = 0;
                    /* 以下是RGBA数据，需要自己补充代码，保存RGBA数据 */
                    bits[y*w+x] |= row_pointers[y][4*x+0] << 16; // red
                    bits[y*w+x] |= row_pointers[y][4*x+1] << 8; // green
                    bits[y*w+x] |= row_pointers[y][4*x+2]; // blue
                    bits[y*w+x] |= row_pointers[y][4*x+3] << 24; // alpha
                }
            }
            break;
            
        case PNG_COLOR_TYPE_RGB:
            for( y=0; y<h; ++y ) {
                for( x=0; x<w; ++x ) {
                    bits[y*w+x] = 0xff000000;
                    bits[y*w+x] |= row_pointers[y][3*x+0] << 16; // red
                    bits[y*w+x] |= row_pointers[y][3*x+1] << 8; // green
                    bits[y*w+x] |= row_pointers[y][3*x+2]; // blue
                }
            }
            break;
            /* 其它色彩类型的图像就不读了 */
        default:
            fclose(fp);
            png_destroy_read_struct( &png_ptr, &info_ptr, 0);
            return 4/* 返回值 */;
    }
    png_destroy_read_struct( &png_ptr, &info_ptr, 0);
    
    texture->datas = (IUINT32**)malloc(1 * sizeof(IUINT32*));
    texture->datas[0] = bits;
    texture->height = h;
    texture->width = w;
    
    return 0;
}

void init_texture() {
    
    // 自建棋盘纹理
    int width = 256, height = 256;
    texture_t *texture = &textures[0];
    IUINT32 *bits = (IUINT32*)malloc(sizeof(IUINT32) * width * height);
    int i, j;
    for (j = 0; j < height; j++) {
        for (i = 0; i < width; i++) {
            int x = i / 32, y = j / 32;
            bits[j*width+i] = ((x + y) & 1)? 0xffffffff : 0xff3fbcef;
        }
    }
    texture->datas_len = 1;
    texture->datas = (IUINT32**)malloc(1 * sizeof(IUINT32*));
    texture->datas[0] = bits;
    texture->width = width;
    texture->height = height;
    texture->use_mipmap = true;
    generate_mipmaps(texture, 1.01);
    texture_count++;
    
    // libpng读取外部纹理

    texture = &textures[1];
    if(load_png_image("mabu", "png", texture) == 0) {
        texture->use_mipmap = true;
        generate_mipmaps(texture, 1.01);
        texture_count++;
    } else {
        printf("load image error, exit!");
        exit(0);
    }
    
    texture = &textures[2];
    if(load_png_image("dimian", "png", texture) == 0) {
        texture->use_mipmap = true;
        generate_mipmaps(texture, 1.01);
        texture_count++;
    } else {
        printf("load image error, exit!");
        exit(0);
    }
}

void free_texture() {
    for(int i = 0; i < texture_count; i++) {
        texture_t *texture = &textures[i];
        for(int j = 0; j < texture->datas_len; j++) {
            IUINT32 *data = texture->datas[j];
            free(data);
        }
        free(texture->datas);
        texture->datas = NULL;
    }
}

void draw_object(device_t *device, object_t *objects, int obj_cnt) {
    for(int i = 0; i < obj_cnt; i++)
    {
        object_t *object = &objects[i];
        if(object->dirty == true) {
            matrix_set_rotate_translate_scale(&object->matrix, &object->axis, object->theta, &object->pos, &object->scale);
            object->dirty = false;
        }
        // 切换纹理
        texture_t *texture = &textures[object->texture_id];
        device_set_texture(device, texture->datas, texture->width, texture->height, false);
        device->material = materials[object->material_id];
        device->transform.world = object->matrix;
        transform_update(&device->transform);
        vertex_t *mesh = object->mesh;
        for(int i = 0; i < object->mesh_num; i+=3)
            clip_polys(device, &mesh[i], &mesh[i+1], &mesh[i+2], false);
    }
}
/*
void draw_shadow(device_t *device, object_t *objects, int obj_cnt) {
    for(int i = 0; i < pointlight_cnt; i++) {
        vector_t pl = pointLights[i].pos;
        for(int j = 0; j < obj_cnt; j++)
        {
            object_t *object = &objects[j];
            if(object->shadow == false)
                continue;
            vertex_t *mesh = object->mesh;
            vertex_t shadow_mesh[object->mesh_num];
            for(int k = 0; k < object->mesh_num; k++) {
                shadow_mesh[k] = mesh[k];
                vector_t vi;
                matrix_apply(&vi, &shadow_mesh[k].pos, &object->matrix);
                
                float t0 = -pl.y / (vi.y - pl.y);
                shadow_mesh[k].pos.x = pl.x + t0 * (vi.x - pl.x);
                shadow_mesh[k].pos.y = 0.02f;
                shadow_mesh[k].pos.z = pl.z + t0 * (vi.z - pl.z);
                shadow_mesh[k].pos.w = 1.0f;
                
                shadow_mesh[k].normal = (vector_t){0.0f, 1.0f, 0.0f, 0.0f};
            }
            for(int k = 0; k < object->mesh_num; k+=3)
                clip_polys(device, &shadow_mesh[k], &shadow_mesh[k+1], &shadow_mesh[k+2], true);
        }
    }
}
*/
void camera_at_zero(device_t *device, const point_t *eye, const vector_t *at, const vector_t *up) {
    matrix_set_lookat(&device->transform.view, eye, at, up);
    transform_update(&device->transform);
}

int screen_keys[512];	// 当前键盘按下状态
float deltaTime = 0.0f;
Uint32 lastFrame = 0;

int main(int argc, char * argv[])
{
    //Start up SDL and create window
    if( !sdl_init(SCREEN_WIDTH, SCREEN_HEIGHT, "lixuefeng") )
    {
        printf( "Failed to initialize!\n" );
    }
    else
    {
        //Main loop flag
        bool quit = false;
        
        device_t device;
        int states[] = { RENDER_STATE_TEXTURE, RENDER_STATE_COLOR, RENDER_STATE_WIREFRAME };
        int indicator = 0;
        int kbhit = 0;
        
        float c_yaw = 0.0f;
        float c_pitch = 0.0f;
        vector_t c_pos = {0.0f, 1.0f, -2.0f, 1.0f};
        vector_t c_front = {0.0f, 0.0f, 1.0f, 0.0f};
        vector_t c_up = {0.0f, 1.0f, 0.0f, 0.0f};
        vector_t c_right = {1.0f, 0.0f, 0.0f, 0.0f};
        vector_t c_worldup = {0.0f, 1.0f, 0.0f, 0.0f};
        float c_movementspeed = 2.0f;
        float c_mouse_sensitivity = 0.7f;
        float c_zoom = 45.0f;
        float c_lastX = SCREEN_WIDTH >> 1, c_lastY = SCREEN_HEIGHT >> 1;
        bool firstMouse = true;
        bool c_dirty = true;
        
        memset(screen_keys, 0, sizeof(int) * 512);
        device_init(&device, REAL_WIDTH, REAL_HEIGHT, 3.1415926 * 0.5f, 0.1f, 500.0f, NULL);
        
        init_texture();
        device.render_state = RENDER_STATE_TEXTURE;
        
        materials[0] = (material_t){0.2f, 0.2f, 0.2f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 32.0f};
        material_cnt++;
        materials[1] = (material_t){0.2f, 0.2f, 0.2f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 32.0f};
        material_cnt++;
        materials[2] = (material_t){0.2f, 0.2f, 0.2f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 32.0f};
        material_cnt++;
        materials[3] = (material_t){0.2f, 0.2f, 0.2f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 32.0f};
        material_cnt++;
        
        dirLight = (dirlight_t){0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.05f, 0.05f, 0.05f, 1.0f, 0.4f, 0.4f, 0.4f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f};
//            dirLight = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        
        int i = 0;
        for(i = 0; i < 4; i++)
        {
            pointLights[i] = (pointlight_t){0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
            pointlight_cnt++;
        }
        
        pointLights[0] = (pointlight_t){{0.0f, 6.0f, -1.0f, 1.0f}, {0.0f, 0.0f, 0.0f, 1.0f}, 1.0f, 0.09f, 0.032f, {0.05f, 0.05f, 0.05f, 1.0f}, {0.4f, 0.4f, 0.4f, 1.0f}, {0.5f, 0.5f, 0.5f, 1.0f}};
        pointLights[1] = (pointlight_t){-1.0f, 6.0f, -1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.09f, 0.032f, 0.05f, 0.05f, 0.05f, 1.0f, 0.4f, 0.4f, 0.4f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f};
        pointLights[2] = (pointlight_t){7.0f, -1.0f, -6.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.09f, 0.032f, 0.05f, 0.05f, 0.05f, 1.0f, 0.4f, 0.4f, 0.4f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f};
        pointLights[3] = (pointlight_t){0.0f, 0.0f, -1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.09f, 0.032f, 0.05f, 0.05f, 0.05f, 1.0f, 0.4f, 0.4f, 0.4f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f};

        // init object
        // ground
        object_t *ground = &objects[0];
        ground->pos = (point_t){0, 0, 0, 1};
        ground->scale = (vector_t){10, 1, 10, 0};
        ground->axis = (vector_t){0, 0, 0, 1};
        ground->theta = 0.0f;
        ground->mesh = ground_mesh;
        ground->mesh_num = 6;
        ground->material_id = 0;
        ground->texture_id = 2;
        ground->shadow = false;
        ground->dirty = true;
        object_count++;
        
        // box
        object_t *box = &objects[1];
        box->pos = (point_t){0, 1, 0, 1};
        box->scale = (vector_t){1, 1, 1, 0};
        box->axis = (vector_t){0, 5, 2, 1};
        box->theta = 0.0f;
        box->mesh = box_mesh;
        box->mesh_num = 36;
        box->material_id = 0;
        box->texture_id = 1;
        box->shadow = true;
        box->dirty = true;
        object_count++;
        
        // box
        object_t *box1 = &objects[2];
        box1->pos = (point_t){0, 2, -1, 1};
        box1->scale = (vector_t){0.5, 0.5, 0.5, 0};
        box1->axis = (vector_t){1, 0, 1, 1};
        box1->theta = 0.0f;
        box1->mesh = box_mesh;
        box1->mesh_num = 36;
        box1->material_id = 0;
        box1->texture_id = 0;
        box1->shadow = false;
        box1->dirty = true;
        object_count++;

        //Event handler
        SDL_Event e;
        
        //While application is running
        while( !quit )
        {
            // Set frame time
            Uint32 currentFrame = SDL_GetTicks();
            deltaTime = (currentFrame - lastFrame) * 1.0f /1000;
            lastFrame = currentFrame;
            
            //Handle events on queue
            while( SDL_PollEvent( &e ) != 0 )
            {
                //User requests quit
                if( e.type == SDL_QUIT )
                {
                    quit = true;
                }
                //User presses a key
                else if( e.type == SDL_KEYDOWN )
                {
                    screen_keys[e.key.keysym.scancode] = 1;
                }
                else if(e.type == SDL_KEYUP)
                {
                    screen_keys[e.key.keysym.scancode] = 0;
                }
                else if(e.type == SDL_MOUSEMOTION)
                {
                    if(firstMouse) {
                        c_lastX = e.motion.x;
                        c_lastY = e.motion.y;
                        firstMouse = false;
                    }
                    float xoffset = e.motion.x - c_lastX;
                    float yoffset = e.motion.y - c_lastY;
                    c_lastX = e.motion.x;
                    c_lastY = e.motion.y;
                    
                    xoffset *= c_mouse_sensitivity;
                    yoffset *= c_mouse_sensitivity;
                    
                    c_yaw += xoffset;
                    c_pitch += yoffset;
                    if(c_pitch > 89.0f)
                        c_pitch = 89.0f;
                    if(c_pitch < -89.0f)
                        c_pitch = -89.0f;
                    
                    //c_front = (vector_t){e.motion.x*2.0/SCREEN_WIDTH-1, 1-e.motion.y*2.0/SCREEN_HEIGHT, 1, 1};
                    
                    c_dirty = true;
                }
            }
            
            if (screen_keys[SDL_SCANCODE_W]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = c_front;
                vector_scale(&temp, velocity);
                vector_add(&c_pos, &c_pos, &temp);
                c_dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_S]) {
                
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = c_front;
                vector_scale(&temp, velocity);
                vector_sub(&c_pos, &c_pos, &temp);
                c_dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_A]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp;
                vector_crossproduct(&temp, &c_front, &c_up);
                vector_normalize(&temp);
                vector_scale(&temp, velocity);
                vector_add(&c_pos, &c_pos, &temp);
                c_dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_D]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp;
                vector_crossproduct(&temp, &c_front, &c_up);
                vector_normalize(&temp);
                vector_scale(&temp, velocity);
                vector_sub(&c_pos, &c_pos, &temp);
                c_dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_Q]) {
                box->theta -= 0.04f;
                box->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_E]) {
                box->theta += 0.04f;
                box->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_UP]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = c_up;
                vector_scale(&temp, velocity);
                vector_add(&box->pos, &box->pos, &temp);
                box->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_LEFT]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = c_right;
                vector_scale(&temp, velocity);
                vector_sub(&box->pos, &box->pos, &temp);
                box->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_DOWN]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = c_up;
                vector_scale(&temp, velocity);
                vector_sub(&box->pos, &box->pos, &temp);
                box->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_RIGHT]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = c_right;
                vector_scale(&temp, velocity);
                vector_add(&box->pos, &box->pos, &temp);
                box->dirty = true;
            }
            
            if (screen_keys[SDL_SCANCODE_SPACE]) {
                if (kbhit == 0) {
                    kbhit = 1;
                    if (++indicator >= 3) indicator = 0;
                    device.render_state = states[indicator];
                }
            }   else {
                kbhit = 0;
            }
            
            // box auto rotate
            box->theta -= 0.04f;
            box->dirty = true;
            
            // box auto rotate
            box1->theta += 0.04f;
            box1->dirty = true;
            
            //Clear screen
            SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
            SDL_RenderClear( gRenderer );
            
            device_clear(&device, 1);
            
            if(c_dirty == true) {
                // 利用欧拉角原理来实现摄像机旋转
                c_front.x = sin(angle_to_radian(c_yaw)) * cos(angle_to_radian(c_pitch));
                c_front.y = -sin(angle_to_radian(c_pitch));
                c_front.z = cos(angle_to_radian(c_yaw)) * cos(angle_to_radian(c_pitch));
                
                vector_normalize(&c_front);
                vector_crossproduct(&c_right, &c_worldup, &c_front);
                vector_normalize(&c_right);
                vector_crossproduct(&c_up, &c_front, &c_right);
                vector_normalize(&c_up);
                vector_t at;
                vector_add(&at, &c_pos, &c_front);
                camera_at_zero(&device, &c_pos, &at, &c_up);
                
                matrix_apply(&dirLight.vdir, &dirLight.dir, &device.transform.view);
                for(i = 0; i < pointlight_cnt; i++)
                {
                    matrix_apply(&pointLights[i].vpos, &pointLights[i].pos, &device.transform.view);
                }
                
                c_dirty = false;
            }
            
            draw_object(&device, objects, object_count);
            
            // 渲染阴影
            // draw_shadow(&device, objects, object_count);

            for(int i = 0; i < SCREEN_WIDTH; i++)
            {
                for(int j = 0; j < SCREEN_HEIGHT; j++)
                {
                    IUINT32 color = device.framebuffer[j][i];
                    SDL_SetRenderDrawColor( gRenderer, (0xff<<16&color)>>16, (0xff<<8&color)>>8, 0xff&color, (0xff<<24&color)>>24);
                    SDL_RenderDrawPoint( gRenderer, i, j);
                }
            }

            //Update screen
            SDL_RenderPresent( gRenderer );
        }
        free_texture();
	}
	sdl_close();
	return 0;
}

