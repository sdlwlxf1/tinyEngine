/*This source code copyrighted by Lazy Foo' Productions (2004-2015)
and may not be redistributed without written permission.*/

//Using SDL, SDL_image, standard IO, math, and strings
#include "SDL2/SDL.h"
#include <stdio.h>
#include <string>
#include <cmath>
//#include "mini3d.c"
#include "tiny3D.h"

#define PI 3.141592653
#define angle_to_radian(X) ((X)/180*PI)
#define radian_to_angle(X) ((X)/PI*180)

//Screen dimension constants
const int SCREEN_WIDTH = 600;
const int SCREEN_HEIGHT = 600;

//Starts up SDL and creates window
bool init();

//Loads media
bool loadMedia();

//Frees media and shuts down SDL
void close();

//Loads individual image as texture
//SDL_Texture* loadTexture( std::string path );

//The window we'll be rendering to
SDL_Window* gWindow = NULL;

//The window renderer
SDL_Renderer* gRenderer = NULL;

bool init(int width, int height, const char *title)
{
	//Initialization flag
	bool success = true;

	//Initialize SDL
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

//				//Initialize PNG loading
//				int imgFlags = IMG_INIT_PNG;
//				if( !( IMG_Init( imgFlags ) & imgFlags ) )
//				{
//					printf( "SDL_image could not initialize! SDL_image Error: %s\n", IMG_GetError() );
//					success = false;
//				}
			}
		}
	}

	return success;
}

bool loadMedia()
{
	//Loading success flag
	bool success = true;

	//Nothing to load
	return success;
}

void close()
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

//SDL_Texture* loadTexture( std::string path )
//{
//	//The final texture
//	SDL_Texture* newTexture = NULL;
//
//	//Load image at specified path
//	SDL_Surface* loadedSurface = IMG_Load( path.c_str() );
//	if( loadedSurface == NULL )
//	{
//		printf( "Unable to load image %s! SDL_image Error: %s\n", path.c_str(), IMG_GetError() );
//	}
//	else
//	{
//		//Create texture from surface pixels
//        newTexture = SDL_CreateTextureFromSurface( gRenderer, loadedSurface );
//		if( newTexture == NULL )
//		{
//			printf( "Unable to create texture from %s! SDL Error: %s\n", path.c_str(), SDL_GetError() );
//		}
//
//		//Get rid of old loaded surface
//		SDL_FreeSurface( loadedSurface );
//	}
//
//	return newTexture;
//}

//=====================================================================
// 主程序
//=====================================================================
vertex_t mesh[36] = {
//    { {  -1, -1, -1, 1 }, { 0, 0 }, { 1.0f, 0.2f, 0.2f }, 1 , { -1, -1, -1, 0}},
//    { { 1, -1, -1, 1 }, { 0, 1 }, { 0.2f, 1.0f, 0.2f }, 1 ,{ 1, -1, -1, 0 }},
//    { { -1,  1, -1, 1 }, { 1, 1 }, { 0.2f, 0.2f, 1.0f }, 1 ,{ -1,  1, -1, 0 }},
//    { {  1,  1, -1, 1 }, { 1, 0 }, { 1.0f, 0.2f, 1.0f }, 1 ,{  1,  1, -1, 0 }},
//    { { -1, -1,  1, 1 }, { 0, 0 }, { 1.0f, 1.0f, 0.2f }, 1 ,{ -1, -1,  1, 0 }},
//    { { 1, -1,  1, 1 }, { 0, 1 }, { 0.2f, 1.0f, 1.0f }, 1 ,{ 1, -1,  1, 0 }},
//    { { -1,  1,  1, 1 }, { 1, 1 }, { 1.0f, 0.3f, 0.3f }, 1 ,{ -1,  1,  1, 0 }},
//    { { 1,  1,  1, 1 }, { 1, 0 }, { 0.2f, 1.0f, 0.3f }, 1 ,{ 1,  1,  1, 0 }},
    
    // Positions                  // Texture Coords  //color           //rhw // Normals
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  0.0f},{ 1.0f, 0.2f, 0.2f }, 1 , { 0.0f,  0.0f,-1.0f,0.0f}},
    {{-0.5f,  0.5f, -0.5f, 1.0f},{ 0.0f,  1.0f},{ 1.0f, 0.2f, 0.2f }, 1 , { 0.0f,  0.0f,-1.0f,0.0f}},
    {{0.5f,  0.5f, -0.5f, 1.0f}, {1.0f,  1.0f}, { 1.0f, 0.2f, 0.2f }, 1 , {0.0f,  0.0f,-1.0f ,0.0f}},
    {{0.5f,  0.5f, -0.5f, 1.0f}, { 1.0f,  1.0f}, { 1.0f, 0.2f, 0.2f }, 1 , {0.0f,  0.0f,-1.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f, 1.0f}, {1.0f,  0.0f}, { 1.0f, 0.2f, 0.2f }, 1 , {0.0f,  0.0f,-1.0f ,0.0f}},
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  0.0f},{ 1.0f, 0.2f, 0.2f }, 1 , { 0.0f,  0.0f,-1.0f,0.0f}},
    
    {{-0.5f, -0.5f,  0.5f, 1.0f},{ 0.0f,  0.0f},{ 0.2f, 1.0f, 0.2f }, 1 , { 0.0f,  0.0f, 1.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f, 1.0f},{ 1.0f,  0.0f}, { 0.2f, 1.0f, 0.2f }, 1 , {0.0f,  0.0f,  1.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},{ 1.0f,  1.0f}, { 0.2f, 1.0f, 0.2f }, 1 , {0.0f,  0.0f,  1.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},{ 1.0f,  1.0f}, { 0.2f, 1.0f, 0.2f }, 1 , {0.0f,  0.0f,  1.0f,0.0f}},
    {{-0.5f,  0.5f,  0.5f, 1.0f},{ 0.0f,  1.0f},{ 0.2f, 1.0f, 0.2f }, 1 , { 0.0f,  0.0f,  1.0f,0.0f}},
    {{-0.5f, -0.5f,  0.5f, 1.0f},{ 0.0f,  0.0f},{ 0.2f, 1.0f, 0.2f }, 1 , { 0.0f,  0.0f,  1.0f,0.0f}},
    
    {{-0.5f,  0.5f,  0.5f, 1.0f}, { 1.0f,  0.0f},{ 0.2f, 0.2f, 1.0f}, 1 , {-1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f, -0.5f, 1.0f},{ 1.0f,  1.0f},{ 0.2f, 0.2f, 1.0f }, 1 , { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  1.0f},{ 0.2f, 0.2f, 1.0f }, 1 , { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f, -0.5f, 1.0f},{ 0.0f,  1.0f},{ 0.2f, 0.2f, 1.0f }, 1 , { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f,  0.5f, 1.0f},{ 0.0f,  0.0f},{ 0.2f, 0.2f, 1.0f}, 1 , { -1.0f,  0.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f,  0.5f,1.0f},{ 1.0f,  0.0f},{ 0.2f, 0.2f, 1.0f }, 1 , { -1.0f,  0.0f,  0.0f,0.0f}},
    
    {{0.5f,  0.5f,  0.5f,1.0f}, { 1.0f,  0.0f}, { 1.0f, 0.2f, 1.0f }, 1 , {1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f,1.0f},{ 0.0f,  0.0f}, { 1.0f, 0.2f, 1.0f }, 1 , { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f,1.0f},{ 0.0f,  1.0f}, { 1.0f, 0.2f, 1.0f }, 1 , { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f,1.0f},{ 0.0f,  1.0f}, { 1.0f, 0.2f, 1.0f }, 1 , { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f, -0.5f,1.0f},{ 1.0f,  1.0f}, { 1.0f, 0.2f, 1.0f }, 1 , { 1.0f,  0.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f,1.0f},{ 1.0f,  0.0f}, { 1.0f, 0.2f, 1.0f }, 1 , { 1.0f,  0.0f,  0.0f,0.0f}},
    
    {{-0.5f, -0.5f, -0.5f,1.0f},{  0.0f,  1.0f},{ 1.0f, 1.0f, 0.2f }, 1 , {  0.0f, -1.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f, -0.5f,1.0f}, { 1.0f,  1.0f},{ 1.0f, 1.0f, 0.2f }, 1 ,  { 0.0f, -1.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f,1.0f}, { 1.0f,  0.0f}, { 1.0f, 1.0f, 0.2f }, 1 , { 0.0f, -1.0f,  0.0f,0.0f}},
    {{0.5f, -0.5f,  0.5f,1.0f}, { 1.0f,  0.0f},{ 1.0f, 1.0f, 0.2f }, 1 ,  { 0.0f, -1.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f,  0.5f,1.0f},{ 0.0f,  0.0f},{ 1.0f, 1.0f, 0.2f }, 1 , {  0.0f, -1.0f,  0.0f,0.0f}},
    {{-0.5f, -0.5f, -0.5f,1.0f},{ 0.0f,  1.0f},{ 1.0f, 1.0f, 0.2f }, 1 , {  0.0f, -1.0f,  0.0f,0.0f}},

    {{-0.5f,  0.5f, -0.5f, 1.0f}, {0.0f,  1.0f},{ 0.2f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f,  0.5f, 1.0f},  {0.0f,  0.0f},{ 0.2f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},  {1.0f,  0.0f},{ 0.2f, 1.0f, 1.0f }, 1 ,  { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f,  0.5f, 1.0f},  {1.0f,  0.0f}, { 0.2f, 1.0f, 1.0f }, 1 , { 0.0f,1.0f,  0.0f,0.0f}},
    {{0.5f,  0.5f, -0.5f, 1.0f},  {1.0f,  1.0f},{ 0.2f, 1.0f, 1.0f }, 1 ,  { 0.0f,1.0f,  0.0f,0.0f}},
    {{-0.5f,  0.5f, -0.5f, 1.0f},  {0.0f,  1.0f},{ 0.2f, 1.0f, 1.0f }, 1 , { 0.0f, 1.0f,  0.0f,0.0f}}

    
};

void draw_plane(device_t *device, int a, int b, int c, int d) {
    vertex_t p1 = mesh[a], p2 = mesh[b], p3 = mesh[c], p4 = mesh[d];
    p1.tc.u = 0, p1.tc.v = 0, p2.tc.u = 0, p2.tc.v = 1;
    p3.tc.u = 1, p3.tc.v = 1, p4.tc.u = 1, p4.tc.v = 0;
    device_draw_primitive(device, &p1, &p2, &p3);
    device_draw_primitive(device, &p3, &p4, &p1);
}

//void draw_box(device_t *device, vector_t *axis, float theta, float x, float y, float z) {
//    matrix_t m;
//    matrix_set_rotate(&m, axis, theta);
//    matrix_set_translate(&m, x, y, z);
//    device->transform.world = m;
//    transform_update(&device->transform);
////    draw_plane(device, 0, 2, 3, 1);
////    draw_plane(device, 0, 4, 6, 2);
////    draw_plane(device, 0, 1, 5, 4);
////    draw_plane(device, 4, 5, 7, 6);
////    draw_plane(device, 1, 3, 7, 5);
////    draw_plane(device, 2, 6, 7, 3);
//    for(int i = 0; i < 36; i+=3)
//        device_draw_primitive(device, &mesh[i], &mesh[i+1], &mesh[i+2]);
//}

void draw_box(device_t *device, const matrix_t *m) {
    device->transform.world = *m;
    transform_update(&device->transform);
    for(int i = 0; i < 36; i+=3)
        device_draw_primitive(device, &mesh[i], &mesh[i+1], &mesh[i+2]);
}


void camera_at_zero(device_t *device, const point_t *eye, const vector_t *at, const vector_t *up) {
    viewPos = *eye;
    matrix_set_lookat(&device->transform.view, eye, at, up);
    transform_update(&device->transform);
}

void init_texture(device_t *device) {
    static IUINT32 texture[256][256];
    int i, j;
    for (j = 0; j < 256; j++) {
        for (i = 0; i < 256; i++) {
            int x = i / 32, y = j / 32;
            texture[j][i] = ((x + y) & 1)? 0xffffff : 0x3fbcef;
        }
    }
    device_set_texture(device, texture, 256 * 4, 256, 256);
}

int screen_keys[512];	// 当前键盘按下状态
float deltaTime = 0.0f;
Uint32 lastFrame = 0;

int main( int argc, char* args[] )
{
    //Start up SDL and create window
    if( !init(SCREEN_WIDTH, SCREEN_HEIGHT, "lixuefeng") )
    {
        printf( "Failed to initialize!\n" );
    }
    else
    {
        //Load media
        if( !loadMedia() )
        {
            printf( "Failed to load media!\n" );
        }
        else
        {
            //Main loop flag
            bool quit = false;
            
            device_t device;
            int states[] = { RENDER_STATE_TEXTURE, RENDER_STATE_COLOR, RENDER_STATE_WIREFRAME };
            int indicator = 0;
            int kbhit = 0;
            float alpha = 0;
            bool box_dirty = true;
            
            matrix_t m;
            point_t pos = {0, 0, 0, 1};
            vector_t scale = {1, 1, 1, 1};
            vector_t axis = {0, 1, 0, 0};
            float theta = 0.0f;
            
            float c_yaw = 0.0f;
            float c_pitch = 0.0f;
            vector_t c_pos = {0.0f, 0.0f, -3.0f, 1.0f};
            vector_t c_front = {0.0f, 0.0f, 1.0f, 0.0f};
            vector_t c_up = {0.0f, 1.0f, 0.0f, 0.0f};
            vector_t c_right = {1.0f, 0.0f, 0.0f, 0.0f};
            vector_t c_worldup = {0.0f, 1.0f, 0.0f, 0.0f};
            float c_movementspeed = 2.0f;
            float c_mouse_sensitivity = 0.25f;
            float c_zoom = 45.0f;
            bool c_dirty = true;
            
            memset(screen_keys, 0, sizeof(int) * 512);
            device_init(&device, SCREEN_WIDTH, SCREEN_HEIGHT, NULL);
            
            init_texture(&device);
            device.render_state = RENDER_STATE_TEXTURE;
            
            material = {0.2f, 0.2f, 0.2f, 0.5f, 0.5f, 0.5f, 1.0f, 1.0f, 1.0f, 32.0f};
            
            //dirLight = {0.0f, 0.0f, 1.0f, 1.0f, 0.05f, 0.05f, 0.05f, 0.4f, 0.4f, 0.4f, 0.5f, 0.5f, 0.5f};
            
            dirLight = {0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
            
//            int i = 0;
//            for(i = 0; i < NR_POINT_LIGHTS; i++)
//            {
//                pointLights[i] = {0.0f, 2.0f, 3.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
//            }
            pointLights[0] = {{0.0f, 0.0f, -1.0f, 1.0f}, 1.0f, 0.09f, 0.032f, {0.05f, 0.05f, 0.05f}, {0.4f, 0.4f, 0.4f}, {0.5f, 0.5f, 0.5f}};
            pointLights[1] = {-1.0f, 0.0f, -1.0f, 1.0f, 1.0f, 0.09f, 0.032f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
            pointLights[2] = {7.0f, -1.0f, -6.0f, 1.0f, 1.0f, 0.09f, 0.032f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
            pointLights[3] = {0.0f, 0.0f, -1.0f, 1.0f, 1.0f, 0.09f, 0.032f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

            
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
                        
                    }
                }
                
                if (screen_keys[SDL_SCANCODE_UP]) {
                    theta += 0.04f;
                    box_dirty = true;
                }
                if (screen_keys[SDL_SCANCODE_DOWN]) {
                    theta -= 0.04f;
                    box_dirty = true;
                }
                if (screen_keys[SDL_SCANCODE_LEFT]) {
                    alpha -= 0.04f;
                    box_dirty = true;
                }
                if (screen_keys[SDL_SCANCODE_RIGHT]) {
                    alpha += 0.04f;
                    box_dirty = true;
                }
                if (screen_keys[SDL_SCANCODE_W]) {
                    float velocity = c_movementspeed * deltaTime;
                    vector_t temp = c_front;
                    vector_scale(&temp, velocity);
                    vector_add(&c_pos, &c_pos, &temp);
                    c_dirty = true;
                }
                if (screen_keys[SDL_SCANCODE_A]) {
                    float velocity = c_movementspeed * deltaTime;
                    vector_t temp = c_right;
                    vector_scale(&temp, velocity);
                    vector_sub(&c_pos, &c_pos, &temp);
                    c_dirty = true;
                }
                if (screen_keys[SDL_SCANCODE_S]) {
                    float velocity = c_movementspeed * deltaTime;
                    vector_t temp = c_front;
                    vector_scale(&temp, velocity);
                    vector_sub(&c_pos, &c_pos, &temp);
                    c_dirty = true;
                }
                if (screen_keys[SDL_SCANCODE_D]) {
                    float velocity = c_movementspeed * deltaTime;
                    vector_t temp = c_right;
                    vector_scale(&temp, velocity);
                    vector_add(&c_pos, &c_pos, &temp);
                    c_dirty = true;
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
                
                //Clear screen
                SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
                SDL_RenderClear( gRenderer );
                
                device_clear(&device, 1);
                
                if(c_dirty == true) {
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
                    c_dirty = false;
                }
                
                pos.x = theta;
                if(box_dirty == true) {
                    matrix_set_rotate_translate_scale(&m, &axis, alpha, &pos, &scale);
                    box_dirty = false;
                }
                
                draw_box(&device, &m);

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
		}
	}

	//Free resources and close SDL
	close();

	return 0;
}

