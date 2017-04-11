/*This source code copyrighted by Lazy Foo' Productions (2004-2015)
and may not be redistributed without written permission.*/

//Using SDL, SDL_image, standard IO, math, and strings
#include "SDL.h"
#include "utils.h"
#include "tiny3D.h"
#include <stdio.h>

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
    {{-0.5f,  0.0f, -0.5f, 1.0f}, {0.0f,  8.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, { 0.0f, 1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.0f,  0.5f, 1.0f},  {0.0f,  0.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, { 0.0f, 1.0f,  0.0f,0.0f}},
    {{0.5f,  0.0f,  0.5f, 1.0f},  {8.0f,  0.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.0f,  0.5f, 1.0f},  {8.0f,  0.0f}, { 0.2f, 1.0f, 1.0f, 1.0f }, { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.0f, -0.5f, 1.0f},  {8.0f,  8.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, { 0.0f,1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.0f, -0.5f, 1.0f},  {0.0f,  8.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, { 0.0f, 1.0f,  0.0f,0.0f}}
};

vertex_t box_mesh[36] = {
    // Positions                  // Texture Coords  //color           //rhw // Normals
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  0.0f},{ 1.0f, 0.2f, 0.2f, 1.0f }, { 0.0f,  0.0f,-1.0f,0.0f}},
    {{-0.5f,  0.5f, -0.5f, 1.0f},{ 0.0f,  1.0f},{ 1.0f, 0.2f, 0.2f, 1.0f }, { 0.0f,  0.0f,-1.0f,0.0f}},
    {{0.5f,  0.5f, -0.5f, 1.0f}, {1.0f,  1.0f}, { 1.0f, 0.2f, 0.2f, 1.0f }, {0.0f,  0.0f,-1.0f ,0.0f}},
    {{0.5f,  0.5f, -0.5f, 1.0f}, { 1.0f,  1.0f}, { 1.0f, 0.2f, 0.2f, 1.0f }, {0.0f,  0.0f,-1.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f, 1.0f}, {1.0f,  0.0f}, { 1.0f, 0.2f, 0.2f, 1.0f }, {0.0f,  0.0f,-1.0f ,0.0f}},
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  0.0f},{ 1.0f, 0.2f, 0.2f, 1.0f }, { 0.0f,  0.0f,-1.0f,0.0f}},
    
    {{-0.5f, -0.5f,  0.5f, 1.0f},{ 0.0f,  0.0f},{ 0.2f, 1.0f, 0.2f, 1.0f }, { 0.0f,  0.0f, 1.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f, 1.0f},{ 1.0f,  0.0f}, { 0.2f, 1.0f, 0.2f, 1.0f }, {0.0f,  0.0f,  1.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},{ 1.0f,  1.0f}, { 0.2f, 1.0f, 0.2f, 1.0f },  {0.0f,  0.0f,  1.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},{ 1.0f,  1.0f}, { 0.2f, 1.0f, 0.2f, 1.0f },  {0.0f,  0.0f,  1.0f,0.0f}},
    {{-0.5f,  0.5f,  0.5f, 1.0f},{ 0.0f,  1.0f},{ 0.2f, 1.0f, 0.2f, 1.0f }, { 0.0f,  0.0f,  1.0f,0.0f}},
    {{-0.5f, -0.5f,  0.5f, 1.0f},{ 0.0f,  0.0f},{ 0.2f, 1.0f, 0.2f, 1.0f },  { 0.0f,  0.0f,  1.0f,0.0f}},
    
    {{-0.5f,  0.5f,  0.5f, 1.0f}, { 1.0f,  0.0f},{ 0.2f, 0.2f, 1.0f, 1.0f},  {-1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f, -0.5f, 1.0f},{ 1.0f,  1.0f},{ 0.2f, 0.2f, 1.0f, 1.0f }, { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  1.0f},{ 0.2f, 0.2f, 1.0f, 1.0f },  { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  1.0f},{ 0.2f, 0.2f, 1.0f, 1.0f }, { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f,  0.5f, 1.0f},{ 0.0f,  0.0f},{ 0.2f, 0.2f, 1.0f, 1.0f},  { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f,  0.5f,1.0f},{ 1.0f,  0.0f},{ 0.2f, 0.2f, 1.0f, 1.0f }, { -1.0f,  0.0f,  0.0f,0.0f}},
    
    {{0.5f,  0.5f,  0.5f,1.0f}, { 1.0f,  0.0f}, { 1.0f, 0.2f, 1.0f, 1.0f },  {1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f,1.0f},{ 0.0f,  0.0f}, { 1.0f, 0.2f, 1.0f, 1.0f }, { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f,1.0f},{ 0.0f,  1.0f}, { 1.0f, 0.2f, 1.0f, 1.0f }, { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f,1.0f},{ 0.0f,  1.0f}, { 1.0f, 0.2f, 1.0f, 1.0f }, { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f, -0.5f,1.0f},{ 1.0f,  1.0f}, { 1.0f, 0.2f, 1.0f, 1.0f }, { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f,1.0f},{ 1.0f,  0.0f}, { 1.0f, 0.2f, 1.0f, 1.0f }, { 1.0f,  0.0f,  0.0f,0.0f}},
    
    {{-0.5f, -0.5f, -0.5f,1.0f},{  0.0f,  1.0f},{ 1.0f, 1.0f, 0.2f, 1.0f },  {  0.0f, -1.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f,1.0f}, { 1.0f,  1.0f},{ 1.0f, 1.0f, 0.2f, 1.0f },  { 0.0f, -1.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f,1.0f}, { 1.0f,  0.0f}, { 1.0f, 1.0f, 0.2f, 1.0f }, { 0.0f, -1.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f,1.0f}, { 1.0f,  0.0f},{ 1.0f, 1.0f, 0.2f, 1.0f },  { 0.0f, -1.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f,  0.5f,1.0f},{ 0.0f,  0.0f},{ 1.0f, 1.0f, 0.2f, 1.0f }, {  0.0f, -1.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f, -0.5f,1.0f},{ 0.0f,  1.0f},{ 1.0f, 1.0f, 0.2f, 1.0f },  {  0.0f, -1.0f,  0.0f,0.0f}},

    {{-0.5f,  0.5f, -0.5f, 1.0f}, {0.0f,  1.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, { 0.0f, 1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f,  0.5f, 1.0f},  {0.0f,  0.0f},{ 0.2f, 1.0f, 1.0f, 1.0f },{ 0.0f, 1.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},  {1.0f,  0.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},  {1.0f,  0.0f}, { 0.2f, 1.0f, 1.0f, 1.0f }, { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f, -0.5f, 1.0f},  {1.0f,  1.0f},{ 0.2f, 1.0f, 1.0f, 1.0f }, { 0.0f,1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f, -0.5f, 1.0f},  {0.0f,  1.0f},{ 0.2f, 1.0f, 1.0f, 1.0f },{ 0.0f, 1.0f,  0.0f,0.0f}}
};

void init_texture() {
    // 自建棋盘纹理
    int width = 256, height = 256;
    texture_t *texture = &textures[texture_count++];
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
    
    // libpng读取外部纹理
    make_texture_by_png("mabu", true);
    make_texture_by_png("dimian", true);
}

void free_textures() {
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

void free_materials() {
    for(int i = 0; i < material_cnt; i++) {
        free_material(&materials[i]);
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
        
        device->transform.model = object->matrix;
        transform_update(&device->transform);
        vertex_t *mesh = object->mesh;
        
        for(int i = 0; i < object->mesh_num; i+=3) {
            // 切换材质组
            if(object->material_ids == NULL)
                device->material = materials[0];
            else
                device->material = materials[object->material_ids[i/3]];
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
    bool sdl_ret = sdl_init(SCREEN_WIDTH, SCREEN_HEIGHT, "lixuefeng");
    
    //Start up SDL and create window
    if(sdl_ret == false)
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
        
        // 加载初始化纹理和网格和材质数据
        init_texture();
        materials[material_cnt++] = (material_t){NULL, {0.2f, 0.2f, 0.2f}, {0.5f, 0.5f, 0.5f}, {0.2f, 0.2f, 0.2f}, {0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, 0.5f}, 32.0f, 1.0f, 1.0f, 1, 1, NULL, -1, NULL, 2, NULL, -1, NULL, -1, NULL, -1, NULL, -1, NULL, -1};
//        materials[material_cnt++] = (material_t){"mabu", {0.2f, 0.2f, 0.2f}, {0.5f, 0.5f, 0.5f}, {0.2f, 0.2f, 0.2f}, {0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, 0.5f}, 32.0f, 1.0f, 1.0f, 1, 1, NULL, -1, "", 1, NULL, -1, NULL, -1, NULL, -1, NULL, -1, NULL, -1};
        
        vertex_t *mesh_nan;
        unsigned long mesh_num_nan;
        int *material_ids_nan;
        unsigned long material_ids_num_nan;
        make_mesh_and_material_by_obj(&mesh_nan, &mesh_num_nan, &material_ids_nan, &material_ids_num_nan, "nanosuit");

        // 缓存
        IUINT32 *framebuffer = (IUINT32*)malloc(REAL_WIDTH * REAL_HEIGHT * sizeof(IUINT32));
        float *zbuffer = (float*)malloc(REAL_HEIGHT * REAL_WIDTH * sizeof(float));
        float *shadowbuffer = (float*)malloc(REAL_HEIGHT * REAL_WIDTH * sizeof(float));
        pshadowbuffer = shadowbuffer;
 
        float c_yaw = 0.0f;
        float c_pitch = 0.0f;
        float c_movementspeed = 2.0f;
        float c_mouse_sensitivity = 0.7f;
        float c_lastX = SCREEN_WIDTH >> 1, c_lastY = SCREEN_HEIGHT >> 1;
        bool firstMouse = true;
        
        memset(screen_keys, 0, sizeof(int) * 512);

        dirLight = (dirlight_t){{0.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f, 0.0f}, false};
        dirLight = (dirlight_t){{0.0f, -1.0f, 1.0f, 0.0f}, {0.3f, 0.3f, 0.3f, 1.0f}, {0.8f, 0.8f, 0.8f, 1.0f}, {0.3f, 0.3f, 0.3f, 1.0f}, true};
        if(dirLight.shadow == true)
        {
            // 影子摄像机
            camera_t *camera = &cameras[camera_count];
            camera->pos = (vector_t){0.0f, 3.0f, -3.0f, 1.0f};
            camera->front = dirLight.dir;
            camera->worldup = (vector_t){0.0f, 1.0f, 0.0f, 0.0f};
            camera->fovy = 3.1415926 * 0.5f;
            camera->zn = 0.1f;
            camera->zf = 15.0f;
            camera->width = REAL_WIDTH;
            camera->height = REAL_HEIGHT;
            camera->aspect = (float)REAL_WIDTH / (float)REAL_HEIGHT;
            camera->projection = orthographic;
            
            camera->left = -3.0f;
            camera->right = 3.0f;
            camera->bottom = -3.0f;
            camera->top = 3.0f;
            camera->dirty = true;
            camera_init_projection(camera);
            camera_count++;
        }

//        for(int i = 0; i < 4; i++)
//            pointLights[pointlight_cnt++] = (pointlight_t){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, false};
        
        pointLights[pointlight_cnt++] = (pointlight_t){{0.0f, 6.0f, -1.0f, 1.0f}, 1.0f, 0.09f, 0.032f, {0.6f, 0.6f, 0.6f, 1.0f}, {0.8f, 0.8f, 0.8f, 1.0f}, {0.7f, 0.7f, 0.7f, 1.0f}, false};
        pointLights[pointlight_cnt++] = (pointlight_t){{0.0f, 6.0f, 2.0f, 1.0f}, 1.0f, 0.09f, 0.032f, {0.6f, 0.6f, 0.6f, 1.0f}, {0.8f, 0.8f, 0.8f, 1.0f}, {0.6f, 0.6f, 0.6f, 1.0f}, false};
//        pointLights[pointlight_cnt++] = (pointlight_t){{0.0f, 6.0f, -1.0f, 1.0f}, 1.0f, 0.09f, 0.032f, {0.6f, 0.6f, 0.6f, 1.0f}, {0.8f, 0.8f, 0.8f, 1.0f}, {0.3f, 0.3f, 0.3f, 1.0f}, false};
//        pointLights[pointlight_cnt++] = (pointlight_t){{0.0f, 6.0f, -1.0f, 1.0f}, 1.0f, 0.09f, 0.032f, {0.6f, 0.6f, 0.6f, 1.0f}, {0.8f, 0.8f, 0.8f, 1.0f}, {0.3f, 0.3f, 0.3f, 1.0f}, false};
        
        // 主摄像机
        camera_t *camera = &cameras[camera_count];
        camera->main = true;
        camera->pos = (vector_t){0.0f, 2.0f, -2.5f, 1.0f};
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
        
        device_set_framebuffer(&device, framebuffer);
        device_set_zbuffer(&device, zbuffer);
        device_set_shadowbuffer(&device, shadowbuffer);
        
        device_set_background(&device, 0x55555555);
        
        device_set_camera(&device, camera);
        transform_update(&device.transform);

        // init object
        // ground
        object_t *ground = &objects[object_count++];
        ground->pos = (point_t){0, 0, 0, 1};
        ground->scale = (vector_t){20, 1, 20, 0};
        ground->axis = (vector_t){0, 0, 0, 1};
        ground->theta = 0.0f;
        ground->mesh = ground_mesh;
        ground->mesh_num = 6;
        ground->material_ids = NULL;
        ground->texture_id = 1;
        ground->shadow = false;
        ground->dirty = true;
        
        // box
        object_t *box = &objects[object_count++];
        box->pos = (point_t){-1, 0, 0, 1};
        box->scale = (vector_t){0.1, 0.1, 0.1, 0};
        box->axis = (vector_t){0, 1, 0, 1};
        box->theta = 0.0f;
        box->mesh = mesh_nan;
        box->mesh_num = mesh_num_nan;
        box->material_ids = material_ids_nan;
        box->texture_id = 1;
        box->shadow = true;
        box->dirty = true;
        
        // box
        object_t *box1 = &objects[object_count++];
        box1->pos = (point_t){0, 2, -1, 1};
        box1->scale = (vector_t){0.5, 0.5, 0.5, 0};
        box1->axis = (vector_t){1, 0, 1, 1};
        box1->theta = 0.0f;
        box1->mesh = box_mesh;
        box1->mesh_num = 36;
        box1->material_ids = NULL;
        box1->texture_id = 0;
        box1->shadow = false;
        box1->dirty = true;
        
        object_t *controlObj = box1;

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
                controlObj->theta -= 0.04f;
                controlObj->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_E]) {
                controlObj->theta += 0.04f;
                controlObj->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_UP]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = {0.0f, 1.0f, 0.0f, 0.0f};
                vector_scale(&temp, velocity);
                vector_add(&controlObj->pos, &controlObj->pos, &temp);
                controlObj->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_LEFT]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = {-1.0f, 0.0f, 0.0f, 0.0f};
                vector_scale(&temp, velocity);
                vector_add(&controlObj->pos, &controlObj->pos, &temp);
                controlObj->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_DOWN]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = {0.0f, -1.0f, 0.0f, 0.0f};
                vector_scale(&temp, velocity);
                vector_add(&controlObj->pos, &controlObj->pos, &temp);
                controlObj->dirty = true;
            }
            if (screen_keys[SDL_SCANCODE_RIGHT]) {
                float velocity = c_movementspeed * deltaTime;
                vector_t temp = {1.0f, 0.0f, 0.0f, 0.0f};
                vector_scale(&temp, velocity);
                vector_add(&controlObj->pos, &controlObj->pos, &temp);
                controlObj->dirty = true;
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
            device_set_shadowbuffer(&device, shadowbuffer);
            
            device_clear(&device);
            
            if(camera->dirty)
                camera_init_by_euler(camera, c_yaw, c_pitch);
            
            for(int i = 0; i < camera_count; i++)
            {
                camera_t *camera = &cameras[i];
                if(camera->main == true) {
                    device.cull = 1;
                    device_set_framebuffer(&device, framebuffer);
                    device_set_zbuffer(&device, zbuffer);
                    device_set_shadowbuffer(&device, NULL);
                } else {
                    device.cull = 2;
                    device_set_framebuffer(&device, NULL);
                    device_set_zbuffer(&device, NULL);
                    device_set_shadowbuffer(&device, shadowbuffer);
                }
                
                if(camera->dirty == true) {
                    camera_update(camera);
                    camera->dirty = false;
                }
                device_set_camera(&device, camera);
                transform_update(&device.transform);
                draw_object(&device, objects, object_count);
            }
 
            // 渲染阴影 采用创建mesh方式，只能投影到平面上
//             draw_shadow(&device, objects, object_count);

            
            for(int y = 0; y < SCREEN_HEIGHT; y++)
            {
                for(int x = 0; x < SCREEN_WIDTH; x++)
                {
                    IUINT32 color = framebuffer[y * REAL_WIDTH + x];
                    SDL_SetRenderDrawColor(gRenderer, (0xff<<16&color)>>16, (0xff<<8&color)>>8, 0xff&color, (0xff<<24&color)>>24);
                    SDL_RenderDrawPoint(gRenderer, x, y);
                    
//                    float depth = shadowbuffer[y * REAL_WIDTH + x];
//                    SDL_SetRenderDrawColor(gRenderer, depth*255, depth*255, depth*255, 255);
//                    SDL_RenderDrawPoint(gRenderer, x, y);
                }
            }

            //Update screen
            SDL_RenderPresent( gRenderer );
        }
        free(mesh_nan);
        free(material_ids_nan);
        free_materials();
        free_textures();
        free(framebuffer);
        free(zbuffer);
        free(shadowbuffer);
	}
	sdl_close();
	return 0;
}

