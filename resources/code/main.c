/*This source code copyrighted by Lazy Foo' Productions (2004-2015)
and may not be redistributed without written permission.*/

//Using SDL, SDL_image, standard IO, math, and strings
#include <stdbool.h>
#include <stdio.h>
#include "SDL.h"
#include "utils.h"


//Screen dimension constants
const int SCREEN_WIDTH = 400;
const int SCREEN_HEIGHT = 400;

const int REAL_WIDTH = 400;
const int REAL_HEIGHT = 400;


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
    {{-0.5f,  0.0f, -0.5f, 1.0f}, {0.0f,  8.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.0f,  0.5f, 1.0f},  {0.0f,  0.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}},
    {{0.5f,  0.0f,  0.5f, 1.0f},  {8.0f,  0.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 ,  { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.0f,  0.5f, 1.0f},  {8.0f,  0.0f}, { 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.0f, -0.5f, 1.0f},  {8.0f,  8.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 ,  { 0.0f,1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.0f, -0.5f, 1.0f},  {0.0f,  8.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}}
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
        
        for(int i = 0; i < object->mesh_num; i+=3) {
//            vertex_t p1 = mesh[i], p2 = mesh[i+1], p3 = mesh[i+2];
//            matrix_apply(&p1.pos, &p1.pos, &device->transform.transform_wv);
//            matrix_apply(&p2.pos, &p2.pos, &device->transform.transform_wv);
//            matrix_apply(&p3.pos, &p3.pos, &device->transform.transform_wv);
//            vertex_t h1 = mesh[i], h2 = mesh[i+1], h3 = mesh[i+2];
//            matrix_apply(&h1.pos, &h1.pos, &device->transform.world);
//            matrix_apply(&h2.pos, &h2.pos, &device->transform.world);
//            matrix_apply(&h3.pos, &h3.pos, &device->transform.world);
//            device_draw_primitive(device, &p1, &p2, &p3, &h1, &h2, &h3);
            clip_polys(device, &mesh[i], &mesh[i+1], &mesh[i+2], false);
        }
        
        
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
        device_init(&device);
        int states[] = { RENDER_STATE_TEXTURE, RENDER_STATE_COLOR, RENDER_STATE_WIREFRAME };
        int indicator = 0;
        int kbhit = 0;
        
        // 缓存
        IUINT32 framebuffer[REAL_HEIGHT][REAL_WIDTH];
        float zbuffer[REAL_HEIGHT][REAL_WIDTH];
        float shadowbuffer[REAL_HEIGHT][REAL_WIDTH];
        pshadowbuffer = (float*)shadowbuffer;

        float c_yaw = 0.0f;
        float c_pitch = 0.0f;
        float c_movementspeed = 2.0f;
        float c_mouse_sensitivity = 0.7f;
        float c_lastX = SCREEN_WIDTH >> 1, c_lastY = SCREEN_HEIGHT >> 1;
        bool firstMouse = true;
        
        memset(screen_keys, 0, sizeof(int) * 512);
        
        init_texture();
        
        materials[0] = (material_t){0.2f, 0.2f, 0.2f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f, 0.2f, 0.2f, 0.2f, 1.0f, 32.0f};
        material_cnt++;
        materials[1] = (material_t){0.2f, 0.2f, 0.2f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f, 0.05f, 0.05f, 0.05f, 1.0f, 32.0f};
        material_cnt++;
        materials[2] = (material_t){0.2f, 0.2f, 0.2f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f, 0.2f, 0.2f, 0.2f, 1.0f, 32.0f};
        material_cnt++;
        materials[3] = (material_t){0.2f, 0.2f, 0.2f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f, 0.2f, 0.2f, 0.2f, 1.0f, 32.0f};
        material_cnt++;
        
        dirLight = (dirlight_t){{0.0f, -1.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 1.0f}, {1.0f, 1.0f, 1.0f, 1.0f}, {0.0f, 0.0f, 0.0f, 1.0f}, true};
        if(dirLight.shadow == true)
        {
            // 影子摄像机
            camera_t *camera = &cameras[camera_count];
            camera->pos = (vector_t){0.0f, 5.0f, -5.0f, 1.0f};
            camera->front = dirLight.dir;
            camera->worldup = (vector_t){0.0f, 1.0f, 0.0f, 0.0f};
            camera->fovy = 3.1415926 * 0.5f;
            camera->zn = 0.1f;
            camera->zf = 20.0f;
            camera->width = REAL_WIDTH;
            camera->height = REAL_HEIGHT;
            camera->aspect = (float)REAL_WIDTH / (float)REAL_HEIGHT;
            camera->projection = orthographic;
            
            camera->left = -5.0f;
            camera->right = 5.0f;
            camera->bottom = -5.0f;
            camera->top = 5.0f;
            camera->dirty = true;
            
            camera_init_projection(camera);
            
            camera_count++;
        }
//        dirLight = (dirlight_t){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        
        for(int i = 0; i < 4; i++)
        {
            pointLights[i] = (pointlight_t){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, false};
            pointlight_cnt++;
        }
        
        pointLights[0] = (pointlight_t){{0.0f, 6.0f, -1.0f, 1.0f}, {0.0f, 0.0f, 0.0f, 1.0f}, 1.0f, 0.09f, 0.032f, {0.05f, 0.05f, 0.05f, 1.0f}, {0.4f, 0.4f, 0.4f, 1.0f}, {0.2f, 0.2f, 0.2f, 1.0f}, false};
//        pointLights[1] = (pointlight_t){-1.0f, 6.0f, -1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.09f, 0.032f, 0.05f, 0.05f, 0.05f, 1.0f, 0.4f, 0.4f, 0.4f, 1.0f, 0.2f, 0.2f, 0.2f, 1.0f};
//        pointLights[2] = (pointlight_t){7.0f, -1.0f, -6.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.09f, 0.032f, 0.05f, 0.05f, 0.05f, 1.0f, 0.4f, 0.4f, 0.4f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f};
//        pointLights[3] = (pointlight_t){0.0f, 0.0f, -1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.09f, 0.032f, 0.05f, 0.05f, 0.05f, 1.0f, 0.4f, 0.4f, 0.4f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f};
        
        // 主摄像机
        camera_t *camera = &cameras[camera_count];
        camera->main = true;
        camera->pos = (vector_t){0.0f, 2.0f, -3.0f, 1.0f};
        camera->front = (vector_t){0.0f, 0.0f, 1.0f, 0.0f};
        camera->worldup = (vector_t){0.0f, 1.0f, 0.0f, 0.0f};
        camera->fovy = 3.1415926 * 0.5f;
        camera->zn = 0.1f;
        camera->zf = 500.0f;
        camera->width = REAL_WIDTH;
        camera->height = REAL_HEIGHT;
        camera->aspect = (float)REAL_WIDTH / (float)REAL_HEIGHT;
        camera->projection = perspective;
        camera->left = -1.0f;
        camera->right = 1.0f;
        camera->bottom = -1.0f;
        camera->top = 1.0f;
        camera_init_projection(camera);
        camera->dirty = true;
        camera_count++;
        
        device_set_framebuffer(&device, (IUINT32*)framebuffer);
        device_set_zbuffer(&device, (float*)zbuffer);
        device_set_shadowbuffer(&device, (float*)shadowbuffer);
        
        device_set_camera(&device, camera);
        transform_update(&device.transform);

        // init object
        // ground
        object_t *ground = &objects[0];
        ground->pos = (point_t){0, 0, 0, 1};
        ground->scale = (vector_t){20, 1, 20, 0};
        ground->axis = (vector_t){0, 0, 0, 1};
        ground->theta = 0.0f;
        ground->mesh = ground_mesh;
        ground->mesh_num = 6;
        ground->material_id = 1;
        ground->texture_id = 1;
        ground->shadow = false;
        ground->dirty = true;
        object_count++;
        
        // box
        object_t *box = &objects[1];
        box->pos = (point_t){0, 0, 0, 1};
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
            deltaTime = (currentFrame - lastFrame) * 1.0f / 1000;
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
                    
                    camera->dirty = true;
                }
            }
            
            if (screen_keys[SDL_SCANCODE_W]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = camera->front;
                vector_scale(&temp, velocity);
                vector_add(&camera->pos, &camera->pos, &temp);
                camera->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_S]) {
                
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = camera->front;
                vector_scale(&temp, velocity);
                vector_sub(&camera->pos, &camera->pos, &temp);
                camera->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_A]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp;
                vector_crossproduct(&temp, &camera->front, &camera->worldup);
                vector_normalize(&temp);
                vector_scale(&temp, velocity);
                vector_add(&camera->pos, &camera->pos, &temp);
                camera->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_D]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp;
                vector_crossproduct(&temp, &camera->front, &camera->worldup);
                vector_normalize(&temp);
                vector_scale(&temp, velocity);
                vector_sub(&camera->pos, &camera->pos, &temp);
                camera->dirty = true;
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
                vector_t temp = {0.0f, 1.0f, 0.0f, 0.0f};
                vector_scale(&temp, velocity);
                vector_add(&box->pos, &box->pos, &temp);
                box->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_LEFT]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = {-1.0f, 0.0f, 0.0f, 0.0f};
                vector_scale(&temp, velocity);
                vector_add(&box->pos, &box->pos, &temp);
                box->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_DOWN]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = {0.0f, -1.0f, 0.0f, 0.0f};
                vector_scale(&temp, velocity);
                vector_add(&box->pos, &box->pos, &temp);
                box->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_RIGHT]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = {1.0f, 0.0f, 0.0f, 0.0f};
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
            } else {
                kbhit = 0;
            }
            
            // box auto rotate
            box->theta -= 0.04f;
            box->dirty = true;
            
            // box auto rotate
            box1->theta += 0.04f;
            box1->dirty = true;
            
            // Clear screen
            SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
            SDL_RenderClear( gRenderer );
            
            // shadowbuffer在这里设置是为了清空buffer
            device_set_shadowbuffer(&device, (float*)shadowbuffer);
            
            device_clear(&device);
            
            if(camera->dirty)
                camera_init_by_euler(camera, c_yaw, c_pitch);
            
            for(int i = 0; i < camera_count; i++)
            {
                camera_t *camera = &cameras[i];
                if(camera->main == true) {
                    device.cull = 1;
                    device_set_framebuffer(&device, (IUINT32*)framebuffer);
                    device_set_zbuffer(&device, (float*)zbuffer);
                    device_set_shadowbuffer(&device, NULL);
                } else {
                    device.cull = 2;
                    device_set_framebuffer(&device, NULL);
                    device_set_zbuffer(&device, NULL);
                    device_set_shadowbuffer(&device, (float*)shadowbuffer);
                }
                device_set_camera(&device, camera);
                if(camera->dirty == true) {
                    camera_update(camera);
                    device.transform.view = camera->view_matrix;
                    transform_update(&device.transform);
                    
                    if(camera->main == true) {
                        matrix_apply(&dirLight.vdir, &dirLight.dir, &device.transform.view);
                        vector_normalize(&dirLight.vdir);
                        for(int j = 0; j < pointlight_cnt; j++)
                        {
                            matrix_apply(&pointLights[i].vpos, &pointLights[i].pos, &device.transform.view);
                            vector_normalize(&pointLights[i].vpos);
                        }
                    }
    
                    camera->dirty = false;
                } else {
                    transform_update(&device.transform);
                }
                draw_object(&device, objects, object_count);
            }
 
            // 渲染阴影
            // draw_shadow(&device, objects, object_count);

            
            for(int y = 0; y < SCREEN_HEIGHT; y++)
            {
                for(int x = 0; x < SCREEN_WIDTH; x++)
                {
                    IUINT32 color = framebuffer[y][x];// * 255
                    SDL_SetRenderDrawColor(gRenderer, (0xff<<16&color)>>16, (0xff<<8&color)>>8, 0xff&color, (0xff<<24&color)>>24);
                    SDL_RenderDrawPoint(gRenderer, x, y);
                    //shadowbuffer[y][x] = 1.0f;
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

