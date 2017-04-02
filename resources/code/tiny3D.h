// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#ifndef tiny3D_h
#define tiny3D_h

#include <stdbool.h>

typedef unsigned int IUINT32;
#define PI 3.141592653
#define angle_to_radian(X) ((X)/180*PI)
#define radian_to_angle(X) ((X)/PI*180)
#define CMID(x, min, max) (((x) < (min)) ? (min) : (((x) > (max)) ? (max) : (x)))
#define interp(a, b, t) ((a) + ((b) - (a)) * (t))
extern int logbase2ofx(int n);

// modules:
// 1. math library
//  1). vector
typedef struct {
    float x, y, z, w;
} vector_t;
typedef vector_t point_t;

//      1)). length
float vector_length(vector_t *v);
//      2)). add
void vector_add(vector_t *c, const vector_t *a, const vector_t *b);
//      3)). sub
void vector_sub(vector_t *c, const vector_t *a, const vector_t *b);
//      4)). scale
void vector_scale(vector_t *v, float k);
//      4)). inverse
void vector_inverse(vector_t *v);
//      4)). dotproduct
float vector_dotproduct(const vector_t *a, const vector_t *b);
//      5)). crossproduct
void vector_crossproduct(vector_t *c, const vector_t *a, const vector_t *b);
//      6)). interp t[0, 1]
void vector_interp(vector_t *c, const vector_t *a, const vector_t *b, float t);
void vector_clone(vector_t *dest, const vector_t *src);
//      6)). reflect
void vector_reflect(vector_t *r, const vector_t *v, const vector_t *n);
//      7)). normalize
void vector_normalize(vector_t *v);
//      8)). interpolating
void vector_interpolating(vector_t *dest, const vector_t *src1, const vector_t *src2, const vector_t *src3, float a, float b, float c);

typedef struct {
    float m[4][4];
} matrix_t;

//      1)). add
void matrix_add(matrix_t *c, const matrix_t *a, const matrix_t *b);
//      2)). sub
void matrix_sub(matrix_t *c, const matrix_t *a, const matrix_t *b);
//      3)). mul
void matrix_mul(matrix_t *c, const matrix_t *a, const matrix_t *b);
//      4)). scale
void matrix_scale(matrix_t *m, float k);
//      5)). inverse
void matrix_inverse(matrix_t *m);
//      6)). transpose
void matrix_transpose(matrix_t *m);
//      5)). apply v = v * m
void matrix_apply(vector_t *y, const vector_t *x, const matrix_t *m);
//      7)). clone 
void matrix_clone(matrix_t *dest, const matrix_t *src);
//      6)). set_identity
void matrix_set_identity(matrix_t *m);
//      7)). set_zero
void matrix_set_zero(matrix_t *m);
//      8)). set_translate
void matrix_set_translate(matrix_t *m, float dx, float dy, float dz);
//      9)). set_scale
void matrix_set_scale(matrix_t *m, float sx, float sy, float sz);
//      10)).set_rotate m, x, y, z, theta
void matrix_set_rotate(matrix_t *m, const vector_t *v, float theta);
void matrix_set_rotate_translate_scale(matrix_t *m, const vector_t *axis, float theta, const point_t *pos, const vector_t *scale);
void matrix_set_axis(matrix_t *m, const vector_t *xaxis, const vector_t *yaxis, const vector_t *zaxis, const point_t *pos);
//      11)).set_lookat m, eye, at, up
// zaxis = normal(At - Eye)
// xaxis = normal(cross(Up, zaxis))
// yaxis = cross(zaxis, xaxis)

// xaxis.x           yaxis.x           zaxis.x          0
// xaxis.y           yaxis.y           zaxis.y          0
// xaxis.z           yaxis.z           zaxis.z          0
//-dot(xaxis, eye)  -dot(yaxis, eye)  -dot(zaxis, eye)  1
void matrix_set_lookat(matrix_t *m, const vector_t *eye, const vector_t *at, const vector_t *up);
//      12)).set_perspective m, fovy, aspect, zn, zf
// zoom = 1 / tan(fov/2)
// zoomy = 1 / tan(fovy/2)
// zoomx = zoomy * aspect
// zoomx    0       0               0
// 0        zoomy   0               0
// 0        0       zf/(zf-zn)      1
// 0        0       zn*zf/(zn-zf)   0
void matrix_set_perspective(matrix_t *m, float fovy, float aspect, float zn, float zf);
//      13)).set_ortho m, left, right, bottom, top, near, far
// 2/(r-l)      0            0           0
// 0            2/(t-b)      0           0
// 0            0            1/(zf-zn)   0
// (l+r)/(l-r)  (t+b)/(b-t)  zn/(zn-zf)  1
void matrix_set_ortho(matrix_t *m, float l, float r, float b, float t, float zn, float zf);
// I have realized the math library. Then I should realize the core data structure of computer graphics.
// 2. transform_t
typedef struct {
    matrix_t model;
    matrix_t view;
    matrix_t view_r;       // view reverse
    matrix_t projection;
    matrix_t vp;           // view * projection
    matrix_t mv;           // model * view
    matrix_t mvp;          // model * view * projection
} transform_t;

//  1). transform_update (world * view * projection
void transform_update(transform_t *ts);
//  3). transform_apply
void transform_apply(const transform_t *ts, vector_t *y, const vector_t *x);
//  4). transform_check_cvv(v)
int transform_check_cvv(const vector_t *v);
//  5). transform_homogenize(ts, y, x)
void transform_homogenize(vector_t *y, const vector_t *x, float width, float height);
//  6). transform_homogenize(ts, y, x)
void transform_homogenize_reverse(vector_t *y, const vector_t *x, float w, float width, float height);

// 3. geometry calculation (vertex, scanline, border check, rect and so on)
typedef struct { float r, g, b, a; } color_t;
void color_init(color_t *c);
//      4)). product
void color_product(color_t *c, const color_t *a, const color_t *b);
void color_scale(color_t *c, float k);
void color_add(color_t *c, const color_t *a, const color_t *b);
void color_sub(color_t *c, const color_t *a, const color_t *b);
void color_interpolating(color_t *dest, const color_t *src1, const color_t *src2, const color_t *src3, float a, float b, float c);


typedef struct {
    char *name;
    
    color_t ambient;
    color_t diffuse;
    color_t specular;
    color_t transmittance;
    color_t emission;
    float shininess;
    float ior;      /* index of refraction */
    float dissolve;
    int illum;
    
    int pad0;
    
    char *ambient_texname;            /* map_Ka */
    int ambient_tex_id;
    char *diffuse_texname;            /* map_Kd */
    int diffuse_tex_id;
    char *specular_texname;           /* map_Ks */
    int specular_tex_id;
    char *specular_highlight_texname; /* map_Ns */
    int specular_highlight_tex_id;
    char *bump_texname;               /* map_bump, bump */
    int bump_tex_id;
    char *displacement_texname;       /* disp */
    int displacement_tex_id;
    char *alpha_texname;              /* map_d */
    int alpha_tex_id;
} material_t;
#define NUM_MATERIAL 100
extern material_t materials[NUM_MATERIAL];
extern int material_cnt;

void free_material(material_t *material);

typedef struct {
    point_t pos;
    float constant;
    float linear;
    float quadratic;
    color_t ambi;
    color_t diff;
    color_t spec;
    bool shadow;
} pointlight_t;
#define NR_POINT_LIGHTS 100
extern pointlight_t pointLights[NR_POINT_LIGHTS];
extern int pointlight_cnt;

typedef struct {
    vector_t dir;
    color_t ambi;
    color_t diff;
    color_t spec;
    bool shadow;
} dirlight_t;
extern dirlight_t dirLight;

typedef enum {
    perspective,
    orthographic
} PROJECTION;
typedef struct {
    // public
    vector_t pos;
    vector_t front;
    vector_t worldup;
    matrix_t view_matrix;
    matrix_t projection_matrix;
    matrix_t view_matrix_r;
    int width;
    int height;
    float fovy;
    float zn;
    float zf;
    
    float left;
    float right;
    float bottom;
    float top;
    
    bool dirty;
    PROJECTION projection;
    bool main;
    
    // private
    float aspect;
} camera_t;
#define MAX_NUM_CAMERA 10
extern camera_t cameras[MAX_NUM_CAMERA];
extern int camera_count;

// 利用欧拉角原理来实现摄像机旋转
void camera_init_by_euler(camera_t *camera, float yaw, float pitch);

void camera_init_projection(camera_t *camera);

void camera_update(camera_t *camera);

typedef struct { float u, v; } texcoord_t;
void texcoord_add(texcoord_t *c, texcoord_t *a, const texcoord_t *b);
void texcoord_scale(texcoord_t *t, float k);
void texcoord_interpolating(texcoord_t *dest, const texcoord_t *src1, const texcoord_t *src2, const texcoord_t *src3, float a, float b, float c);

typedef struct { float a, b, c, d; } storage_t; // 插值寄存器
void storage_add(storage_t *c, storage_t *a, const storage_t *b);
void storage_scale(storage_t *t, float k);
void storage_interpolating(storage_t *dest, const storage_t *src1, const storage_t *src2, const storage_t *src3, float a, float b, float c);

typedef struct { point_t pos; texcoord_t tc; color_t color; vector_t normal;} vertex_t; // 提供4个额外的插值寄存器

typedef struct { vertex_t v, v1, v2; } edge_t;
typedef struct { float top, bottom; edge_t left, right; } trapezoid_t;
typedef struct { vertex_t v, step; int x, y, w; } scanline_t;

// 注意是除坐标意外颜色和纹理索引除以w
void vertex_rhw_init(vertex_t *v);

void vertex_interp(vertex_t *y, const vertex_t *x1, const vertex_t *x2, float k);

void vertex_division(vertex_t *y, const vertex_t *x1, const vertex_t *x2, float w);

void vertex_add(vertex_t *y, const vertex_t *x);

extern float *pshadowbuffer;

typedef struct {
    int render_state;           // 渲染状态
	transform_t transform;      // 坐标变换器
    material_t material;        // 当前材质
	IUINT32 *framebuffer;       // 像素缓存
	float *zbuffer;             // 深度缓存
    float *shadowbuffer;        // 阴影缓存
	IUINT32 background;         // 背景颜色
	IUINT32 foreground;         // 线框颜色
    camera_t *camera;           // 摄像机
    bool blend;                 // 是否开启混合
    float blend_sfactor;        // 混合源因子
    float blend_dfactor;        // 混合目标因子
    int cull;                   // 0:不裁剪;1:裁剪反面;2:裁剪正面
}	device_t;

#define RENDER_STATE_WIREFRAME      1		// 渲染线框
#define RENDER_STATE_TEXTURE        2		// 渲染纹理
#define RENDER_STATE_COLOR          4		// 渲染颜色

void device_init(device_t *device);

void device_set_background(device_t *device, IUINT32 color);

void device_set_framebuffer(device_t *device, IUINT32 *framebuffer);

void device_set_zbuffer(device_t *device, float *zbuffer);

void device_set_shadowbuffer(device_t *device, float *shadowbuffer);

void device_set_camera(device_t *device, camera_t *camera);

void device_pixel(device_t *device, int x, int y, IUINT32 color);

void device_clear(device_t *device);

// entity object
typedef struct {
    vertex_t *mesh;
    unsigned long mesh_num;
    int* material_ids;
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
extern object_t objects[MAX_NUM_OBJECT];
extern int object_count;

// texture
typedef struct {
    IUINT32 **datas;            // data
    IUINT32 datas_len;
    bool use_mipmap;            // able mipmap
    IUINT32 width;
    IUINT32 height;
} texture_t;
#define MAX_NUM_TEXTURE 100
extern texture_t textures[MAX_NUM_TEXTURE];
extern int texture_count;

void clip_polys(device_t *device, vertex_t *v1, vertex_t *v2, vertex_t *v3, bool world);

typedef struct {
    vector_t pos;
    color_t color;
    vector_t normal;
    vector_t tangent;
    vector_t binormal;
    texcoord_t texcoord;
} a2v;

typedef struct {
    vector_t pos;
    texcoord_t texcoord;
    color_t color;
    vector_t normal;
    vector_t storage0;
    vector_t storage1;
    vector_t storage2;
} v2f;

void vert_shader(device_t *device, a2v *av, v2f *vf);

void frag_shader(device_t *device, v2f *vf, color_t *color);


#endif /* tiny3D_h */
