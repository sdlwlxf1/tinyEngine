// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#include "tiny3D.h"

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void vector_interpolating(vector_t *dest, const vector_t *src1, const vector_t *src2, const vector_t *src3, float a, float b, float c) {
    dest->x = dest->y = dest->z = dest->w = 0.0f;
    vector_t each = *src1;
    vector_scale(&each, a);
    vector_add(dest, dest, &each);
    each = *src2;
    vector_scale(&each, b);
    vector_add(dest, dest, &each);
    each = *src3;
    vector_scale(&each, c);
    vector_add(dest, dest, &each);
}

void texcoord_add(texcoord_t *c, texcoord_t *a, const texcoord_t *b) {
    c->u = a->u + b->u;
    c->v = a->v + b->v;
}
void texcoord_scale(texcoord_t *t, float k) {
    t->u *= k;
    t->v *= k;
}

void texcoord_interpolating(texcoord_t *dest, const texcoord_t *src1, const texcoord_t *src2, const texcoord_t *src3, float a, float b, float c) {
    dest->u = dest->v = 0.0f;
    texcoord_t each = *src1;
    texcoord_scale(&each, a);
    texcoord_add(dest, dest, &each);
    each = *src2;
    texcoord_scale(&each, b);
    texcoord_add(dest, dest, &each);
    each = *src3;
    texcoord_scale(&each, c);
    texcoord_add(dest, dest, &each);
}

void storage_add(storage_t *c, storage_t *a, const storage_t *b) {
    c->a = a->a + b->a;
    c->b = a->b + b->b;
    c->c = a->c + b->c;
    c->d = a->d + b->d;
}

void storage_scale(storage_t *t, float k) {
    t->a *= k;
    t->b *= k;
    t->c *= k;
    t->d *= k;
}

void storage_interpolating(storage_t *dest, const storage_t *src1, const storage_t *src2, const storage_t *src3, float a, float b, float c) {
    dest->a = dest->b = dest->c = dest->d = 0.0f;
    storage_t each = *src1;
    storage_scale(&each, a);
    storage_add(dest, dest, &each);
    each = *src2;
    storage_scale(&each, b);
    storage_add(dest, dest, &each);
    each = *src3;
    storage_scale(&each, c);
    storage_add(dest, dest, &each);
}

void v2f_interpolating(v2f *dest, const v2f *src1, const v2f *src2, const v2f *src3, float a, float b, float c) {
    vector_interpolating(&dest->pos, &src1->pos, &src2->pos, &src3->pos, a, b, c);
    color_interpolating(&dest->color, &src1->color, &src2->color, &src3->color, a, b, c);
    texcoord_interpolating(&dest->texcoord, &src1->texcoord, &src2->texcoord, &src3->texcoord, a, b, c);
    vector_interpolating(&dest->normal, &src1->normal, &src2->normal, &src3->normal, a, b, c);
    vector_interpolating(&dest->storage0, &src1->storage0, &src2->storage0, &src3->storage0, a, b, c);
    vector_interpolating(&dest->storage1, &src1->storage1, &src2->storage1, &src3->storage1, a, b, c);
    vector_interpolating(&dest->storage2, &src1->storage2, &src2->storage2, &src3->storage2, a, b, c);
}

int logbase2ofx(int n) {
    if(n <= 0) return 0;
    int r = 0;
    while(n != 0)
    {
        n = n >> 1;
        r++;
    }
    return r-1;
}

//      1)). length
float vector_length(vector_t *v) {
    return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}
//      2)). add
void vector_add(vector_t *c, const vector_t *a, const vector_t *b) {
    c->x = a->x + b->x; 
    c->y = a->y + b->y; 
    c->z = a->z + b->z;
}
//      3)). sub
void vector_sub(vector_t *c, const vector_t *a, const vector_t *b) {
    c->x = a->x - b->x; 
    c->y = a->y - b->y; 
    c->z = a->z - b->z;
}
//      4)). scale
void vector_scale(vector_t *v, float k) {
    v->x *= k;
    v->y *= k;
    v->z *= k;
    v->w *= k;
}
//      4)). inverse
void vector_inverse(vector_t *v) {
    v->x = -v->x;
    v->y = -v->y;
    v->z = -v->z;
}
//      4)). dotproduct
float vector_dotproduct(const vector_t *a, const vector_t *b) {
    return a->x * b->x + a->y * b->y + a->z * b->z;
}
//      5)). crossproduct
void vector_crossproduct(vector_t *c, const vector_t *a, const vector_t *b) {
    c->x = a->y*b->z - a->z*b->y;
    c->y = a->z*b->x - a->x*b->z;
    c->z = a->x*b->y - a->y*b->x;
}
//      6)). interp t[0, 1]
void vector_interp(vector_t *c, const vector_t *a, const vector_t *b, float t) {
    c->x = interp(a->x, b->x, t);
    c->y = interp(a->y, b->y, t);
    c->z = interp(a->z, b->z, t);
    c->w = interp(a->w, b->w, t);
}
void vector_clone(vector_t *dest, const vector_t *src) {
    dest->x = src->x;
    dest->y = src->y;
    dest->z = src->z;
    dest->w = src->w;
}
//      6)). reflect
void vector_reflect(vector_t *r, const vector_t *v, const vector_t *n) {
    vector_clone(r, n);
    vector_scale(r,  -2 * vector_dotproduct(v, n));
    vector_add(r, r, v);
}
//      7)). normalize
void vector_normalize(vector_t *v) {
    float len = vector_length(v);
    if(fabsf(len) < 1e-6) return;
    float k = 1 / len;
    v->x *= k;
    v->y *= k;
    v->z *= k;
}

//      1)). add
void matrix_add(matrix_t *c, const matrix_t *a, const matrix_t *b) {
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            c->m[i][j] = a->m[i][j] + b->m[i][j];
}
//      2)). sub
void matrix_sub(matrix_t *c, const matrix_t *a, const matrix_t *b) {
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            c->m[i][j] = a->m[i][j] - b->m[i][j];
}
//      3)). mul
void matrix_mul(matrix_t *c, const matrix_t *a, const matrix_t *b) {
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            c->m[i][j] = a->m[i][0]*b->m[0][j] + a->m[i][1]*b->m[1][j] + a->m[i][2]*b->m[2][j] + a->m[i][3]*b->m[3][j];
}
//      4)). scale
void matrix_scale(matrix_t *m, float k) {
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            m->m[i][j] *= k;
}
//      5)). inverse
void matrix_inverse(matrix_t *m) {
    float t[3][6];
    int i, j, k;
    float f;

    for(i = 0; i < 3; i++)
        for(j = 0; j < 6; j++) {
            if(j < 3)
                t[i][j] = m->m[i][j];
            else if(j == i + 3)
                t[i][j] = 1;
            else
                t[i][j] = 0;
        }

    for(i = 0; i < 3; i++) {
        f = t[i][i];
        for(j = 0; j < 6; j++)
            t[i][j] /= f;
        for(j = 0; j < 3; j++) {
            if(j != i) {
                f = t[j][i];
                for(k = 0; k < 6; k++)
                    t[j][k] = t[j][k] - t[i][k] * f;
            }
        }
    }

    for(i = 0; i < 3; i++)
        for(j = 3; j < 6; j++)
            m->m[i][j-3] = t[i][j];
    
    m->m[3][0] = -m->m[3][0];
    m->m[3][1] = -m->m[3][1];
    m->m[3][2] = -m->m[3][2];
}

//      6)). transpose
void matrix_transpose(matrix_t *m) {
    for(int i = 0; i < 3; i++)
        for(int j = i+1; j < 3; j++)
        {
            float temp = m->m[i][j];
            m->m[i][j] = m->m[j][i];
            m->m[j][i] = temp;
        }
}
//      5)). apply v = v * m
void matrix_apply(vector_t *y, const vector_t *x, const matrix_t *m) {
    float X = x->x, Y = x->y, Z = x->z, W = x->w;
    y->x = X * m->m[0][0] + Y * m->m[1][0] + Z * m->m[2][0] + W * m->m[3][0];
    y->y = X * m->m[0][1] + Y * m->m[1][1] + Z * m->m[2][1] + W * m->m[3][1];
    y->z = X * m->m[0][2] + Y * m->m[1][2] + Z * m->m[2][2] + W * m->m[3][2];
    y->w = X * m->m[0][3] + Y * m->m[1][3] + Z * m->m[2][3] + W * m->m[3][3];
}
//      7)). clone 
void matrix_clone(matrix_t *dest, const matrix_t *src) {
    int i, j;
    for(i = 0; i < 4; i++)
        for(j = 0; j < 4; j++)
            dest->m[i][j] = src->m[i][j];
}
//      6)). set_identity
void matrix_set_identity(matrix_t *m) {
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            if(i == j)
                m->m[i][j] = 1.0f;
            else
                m->m[i][j] = 0.0f;
}
//      7)). set_zero
void matrix_set_zero(matrix_t *m) {
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            m->m[i][j] = 0.0f;
}
//      8)). set_translate
void matrix_set_translate(matrix_t *m, float dx, float dy, float dz) {
    matrix_set_identity(m);
    m->m[3][0] = dx;
    m->m[3][1] = dy;
    m->m[3][2] = dz;
}
//      9)). set_scale
void matrix_set_scale(matrix_t *m, float sx, float sy, float sz) {
    matrix_set_identity(m);
    m->m[0][0] = sx;
    m->m[1][1] = sy;
    m->m[2][2] = sz;
}
//      10)).set_rotate m, x, y, z, theta
// x*x*(1-cos0)+cos0    x*y*(1-cos0)+z*sin0   x*z*(1-cos0)-y*sin0
// x*y*(1-cos0)-z*sin0  y*y*(1-cos0)+cos0     y*z*(1-cos0)+x*sin0
// x*z*(1-cos0)+y*sin0  y*z*(1-cos0)-x*sin0   z*z*(1-cos0)+cos0
void matrix_set_rotate(matrix_t *m, const vector_t *v, float theta) {
//    float s = sin(theta);
//    float c = cos(theta);
//
//    float a = 1.0f - c;
//    float ax = a * v->x;
//    float ay = a * v->y;
//    float az = a * v->z;
//
//    m->m[0][0] = ax * v->x + c;
//    m->m[0][1] = ax * v->y + v->z * s;
//    m->m[0][2] = ax * v->z - v->y * s;
//
//    m->m[1][0] = ay * v->x - v->z * s;
//    m->m[1][1] = ay * v->y + c;
//    m->m[1][2] = ay * v->z + v->x * s;
//
//    m->m[2][0] = az * v->x + v->y * s;
//    m->m[2][1] = az * v->y - v->x * s;
//    m->m[2][2] = az * v->z + c;
//
//    m->m[3][0] = m->m[3][1] = m->m[3][2] = 0.0f;
//    m->m[0][3] = m->m[1][3] = m->m[2][3] = 0.0f;
//    m->m[3][3] = 1.0f;
    float x = v->x, y = v->y, z = v->z;
    float qsin = (float)sin(theta * 0.5f);
    float qcos = (float)cos(theta * 0.5f);
    vector_t vec = { x, y, z, 1.0f };
    float w = qcos;
    vector_normalize(&vec);
    x = vec.x * qsin;
    y = vec.y * qsin;
    z = vec.z * qsin;
    m->m[0][0] = 1 - 2 * y * y - 2 * z * z;
    m->m[1][0] = 2 * x * y - 2 * w * z;
    m->m[2][0] = 2 * x * z + 2 * w * y;
    m->m[0][1] = 2 * x * y + 2 * w * z;
    m->m[1][1] = 1 - 2 * x * x - 2 * z * z;
    m->m[2][1] = 2 * y * z - 2 * w * x;
    m->m[0][2] = 2 * x * z - 2 * w * y;
    m->m[1][2] = 2 * y * z + 2 * w * x;
    m->m[2][2] = 1 - 2 * x * x - 2 * y * y;
    m->m[0][3] = m->m[1][3] = m->m[2][3] = 0.0f;
    m->m[3][0] = m->m[3][1] = m->m[3][2] = 0.0f;
    m->m[3][3] = 1.0f;
}
void matrix_set_rotate_translate_scale(matrix_t *m, const vector_t *axis, float theta, const point_t *pos, const vector_t *scale) {
    matrix_set_scale(m, scale->x, scale->y, scale->z);
    matrix_t r, t = *m;
    matrix_set_rotate(&r, axis, theta);
    matrix_mul(m, &t, &r);
    m->m[3][0] = pos->x;
    m->m[3][1] = pos->y;
    m->m[3][2] = pos->z;
}
void matrix_set_axis(matrix_t *m, const vector_t *xaxis, const vector_t *yaxis, const vector_t *zaxis, const point_t *pos) {
    m->m[0][0] = xaxis->x;
    m->m[0][1] = xaxis->y;
    m->m[0][2] = xaxis->z;
    m->m[1][0] = yaxis->x;
    m->m[1][1] = yaxis->y;
    m->m[1][2] = yaxis->z;
    m->m[2][0] = zaxis->x;
    m->m[2][1] = zaxis->y;
    m->m[2][2] = zaxis->z;
    m->m[3][0] = pos->x;
    m->m[3][1] = pos->y;
    m->m[3][2] = pos->z;
    m->m[0][3] = m->m[1][3] = m->m[2][3] = 0.0f;
    m->m[3][3] = 1.0f;
}
//      11)).set_lookat m, eye, at, up
// zaxis = normal(At - Eye)
// xaxis = normal(cross(Up, zaxis))
// yaxis = cross(zaxis, xaxis)

// xaxis.x           yaxis.x           zaxis.x          0
// xaxis.y           yaxis.y           zaxis.y          0
// xaxis.z           yaxis.z           zaxis.z          0
//-dot(xaxis, eye)  -dot(yaxis, eye)  -dot(zaxis, eye)  1
void matrix_set_lookat(matrix_t *m, const vector_t *eye, const vector_t *at, const vector_t *up) {
    vector_t zaxis, xaxis, yaxis;
    vector_sub(&zaxis, at, eye);
    vector_normalize(&zaxis);
    vector_crossproduct(&xaxis, up, &zaxis);
    vector_normalize(&xaxis);
    vector_crossproduct(&yaxis, &zaxis, &xaxis);

    m->m[0][0] = xaxis.x;
    m->m[0][1] = yaxis.x;
    m->m[0][2] = zaxis.x;
    m->m[1][0] = xaxis.y;
    m->m[1][1] = yaxis.y;
    m->m[1][2] = zaxis.y;
    m->m[2][0] = xaxis.z;
    m->m[2][1] = yaxis.z;
    m->m[2][2] = zaxis.z;
    m->m[3][0] = -vector_dotproduct(&xaxis, eye);
    m->m[3][1] = -vector_dotproduct(&yaxis, eye);
    m->m[3][2] = -vector_dotproduct(&zaxis, eye);
    m->m[0][3] = m->m[1][3] = m->m[2][3] = 0.0f;
    m->m[3][3] = 1.0f;
}
//      12)).set_perspective m, fovy, aspect, zn, zf
// zoom = 1 / tan(fov/2)
// zoomy = 1 / tan(fovy/2)
// zoomx = zoomy * aspect
// zoomx    0       0               0
// 0        zoomy   0               0
// 0        0       zf/(zf-zn)      1
// 0        0       zn*zf/(zn-zf)   0
void matrix_set_perspective(matrix_t *m, float fovy, float aspect, float zn, float zf) {
    float zoomy = 1.0f / (float)tan(fovy * 0.5f);
    float zoomx = zoomy * aspect;
    m->m[0][0] = zoomx;
    m->m[1][1] = zoomy;
    m->m[2][2] = zf / (zf - zn);
    m->m[3][2] = -m->m[2][2] * zn;
    m->m[2][3] = 1;
    m->m[0][1] = m->m[0][2] = m->m[0][3] = m->m[1][0] = m->m[1][2] = m->m[1][3] = 0.0f;
    m->m[2][0] = m->m[2][1] = m->m[3][0] = m->m[3][1] = m->m[3][3] = 0.0f;
}
//      13)).set_ortho m, left, right, bottom, top, near, far
// 2/(r-l)      0            0           0
// 0            2/(t-b)      0           0
// 0            0            1/(zf-zn)   0
// (l+r)/(l-r)  (t+b)/(b-t)  zn/(zn-zf)  1
void matrix_set_ortho(matrix_t *m, float l, float r, float b, float t, float zn, float zf) {
    m->m[0][0] = 2.0f / (r - l);
    m->m[1][1] = 2.0f / (t - b);
    m->m[2][2] = 1.0f / (zf - zn);
    m->m[3][0] = (l+r)/(l-r);
    m->m[3][1] = (t+b)/(b-t);
    m->m[3][2] = zn / (zn-zf);
    m->m[3][3] = 1.0f;
    m->m[0][1] = m->m[0][2] = m->m[0][3] = m->m[1][0] = m->m[1][2] = m->m[1][3] = 0.0f;
    m->m[2][0] = m->m[2][1] = m->m[2][3] = 0.0f;
}

//  1). transform_update (world * view * projection
void transform_update(transform_t *ts) {
    matrix_mul(&ts->mv, &ts->model, &ts->view);
    matrix_mul(&ts->vp, &ts->view, &ts->projection);
    matrix_mul(&ts->mvp, &ts->mv, &ts->projection);
}
//  3). transform_apply
void transform_apply(const transform_t *ts, vector_t *y, const vector_t *x) {
    matrix_apply(y, x, &ts->mvp);
}
//  4). transform_check_cvv(v)
int transform_check_cvv(const vector_t *v) {
    float w = v->w;
    int check = 0;
    if(v->z < 0.0f) check |= 1;
    if(v->z > w) check |= 2;
    if(v->x < -w) check |= 4;
    if(v->x > w) check |= 8;
    if(v->y < -w) check |= 16;
    if(v->y > w) check |= 32;
    return check;
}
//  5). transform_homogenize(ts, y, x)
void transform_homogenize(vector_t *y, const vector_t *x, float width, float height) {
    float rhw = 1.0f / x->w;
    y->x = (x->x * rhw + 1.0f) * width * 0.5f;
    y->y = (1.0f - x->y * rhw) * height * 0.5f;
    y->z = x->z * rhw;
    y->w = rhw;
}
//  6). transform_homogenize(ts, y, x)
void transform_homogenize_reverse(vector_t *y, const vector_t *x, float w, float width, float height) {
    y->x = (x->x * 2 / width - 1.0f) * w;
    y->y = (1.0f - x->y * 2 / height) * w;
    y->z = x->z * w;
    y->w = w;
}

void color_init(color_t *c) {
    c->r = 0.0f;
    c->g = 0.0f;
    c->b = 0.0f;
    c->a = 0.0f;
}
//      4)). product
void color_product(color_t *c, const color_t *a, const color_t *b) {
    c->r = a->r * b->r;
    c->g = a->g * b->g;
    c->b = a->b * b->b;
    c->a = a->a * b->a;
}
//      5)). product
void color_product_array(color_t *c, const color_t *a, const float *b) {
    c->r = a->r * b[0];
    c->g = a->g * b[1];
    c->b = a->b * b[2];
}
void color_scale(color_t *c, float k) {
    c->r *= k;
    c->g *= k;
    c->b *= k;
    c->a *= k;
}
void color_add(color_t *c, const color_t *a, const color_t *b) {
    c->r = a->r + b->r;
    c->g = a->g + b->g;
    c->b = a->b + b->b;
    c->a = a->a + b->a;
}
void color_sub(color_t *c, const color_t *a, const color_t *b) {
    c->r = a->r - b->r;
    c->g = a->g - b->g;
    c->b = a->b - b->b;
    c->a = a->a - b->a;
}

void color_interpolating(color_t *dest, const color_t *src1, const color_t *src2, const color_t *src3, float a, float b, float c) {
    dest->r = dest->g = dest->b = dest->a = 0.0f;
    color_t each = *src1;
    color_scale(&each, a);
    color_add(dest, dest, &each);
    each = *src2;
    color_scale(&each, b);
    color_add(dest, dest, &each);
    each = *src3;
    color_scale(&each, c);
    color_add(dest, dest, &each);
}

object_t objects[MAX_NUM_OBJECT];
int object_count = 0;

texture_t textures[MAX_NUM_TEXTURE];
int texture_count = 0;

material_t materials[NUM_MATERIAL];
int material_cnt;

void free_material(material_t *material) {
    free(material->name);
    free(material->bump_texname);
    free(material->alpha_texname);
    free(material->ambient_texname);
    free(material->diffuse_texname);
    free(material->specular_texname);
    free(material->displacement_texname);
    free(material->specular_highlight_texname);
}

pointlight_t pointLights[NR_POINT_LIGHTS];
int pointlight_cnt;

dirlight_t dirLight;

camera_t cameras[MAX_NUM_CAMERA];
int camera_count = 0;

// 利用欧拉角原理来实现摄像机旋转
void camera_init_by_euler(camera_t *camera, float yaw, float pitch) {
    camera->front.x = sin(angle_to_radian(yaw)) * cos(angle_to_radian(pitch));
    camera->front.y = -sin(angle_to_radian(pitch));
    camera->front.z = cos(angle_to_radian(yaw)) * cos(angle_to_radian(pitch));
    vector_normalize(&camera->front);
}

void camera_init_projection(camera_t *camera)
{
    if(camera->projection == perspective)
        matrix_set_perspective(&camera->projection_matrix, camera->fovy, camera->aspect, camera->zn , camera->zf);
    else if(camera->projection == orthographic)
        matrix_set_ortho(&camera->projection_matrix, camera->left, camera->right, camera->bottom, camera->top, camera->zn, camera->zf);
}

void camera_update(camera_t *camera) {
    vector_t right, at, up, front = camera->front;
    vector_crossproduct(&right, &camera->worldup, &camera->front);
    vector_normalize(&right);
    vector_crossproduct(&up, &camera->front, &right);
    vector_normalize(&up);
    vector_add(&at, &camera->pos, &camera->front);
    
    matrix_set_lookat(&camera->view_matrix, &camera->pos,  &at,  &up);
    vector_normalize(&front);
    matrix_set_axis(&camera->view_matrix_r, &right, &up, &front, &camera->pos);
}


// 注意是除坐标意外颜色和纹理索引除以w
void vertex_rhw_init(vertex_t *v) {
    float rhw = 1.0f / v->pos.w;
    v->tc.u *= rhw;
    v->tc.v *= rhw;
    v->color.r *= rhw;
    v->color.g *= rhw;
    v->color.b *= rhw;
}

void vertex_interp(vertex_t *y, const vertex_t *x1, const vertex_t *x2, float k) {
    vector_interp(&y->pos, &x1->pos, &x2->pos, k);
    y->tc.u = interp(x1->tc.u, x2->tc.u, k);
    y->tc.v = interp(x1->tc.v, x2->tc.v, k);
    y->color.r = interp(x1->color.r, x2->color.r, k);
    y->color.g = interp(x1->color.g, x2->color.g, k);
    y->color.b = interp(x1->color.b, x2->color.b, k);
}

void vertex_division(vertex_t *y, const vertex_t *x1, const vertex_t *x2, float w) {
    float inv = 1.0f / w;
    y->pos.x = (x2->pos.x - x1->pos.x) * inv;
    y->pos.y = (x2->pos.y - x1->pos.y) * inv;
    y->pos.z = (x2->pos.z - x1->pos.z) * inv;
    y->pos.w = (x2->pos.w - x1->pos.w) * inv;
    y->tc.u = (x2->tc.u - x1->tc.u) * inv;
    y->tc.v = (x2->tc.v - x1->tc.v) * inv;
    y->color.r = (x2->color.r - x1->color.r) * inv;
    y->color.g = (x2->color.g - x1->color.g) * inv;
    y->color.b = (x2->color.b - x1->color.b) * inv;
}

void vertex_add(vertex_t *y, const vertex_t *x) {
    y->pos.x += x->pos.x;
    y->pos.y += x->pos.y;
    y->pos.z += x->pos.z;
    y->pos.w += x->pos.w;
    y->tc.u += x->tc.u;
    y->tc.v += x->tc.v;
    y->color.r += x->color.r;
    y->color.g += x->color.g;
    y->color.b += x->color.b;
}

// 根据三角形生成 0-2 个梯形，并且返回合法梯形的数量
int trapezoid_init_triangle(trapezoid_t *trap, const vertex_t *p1, const vertex_t *p2, const vertex_t *p3) {
	const vertex_t *p;
	float k, x;

	if (p1->pos.y > p2->pos.y) p = p1, p1 = p2, p2 = p;
	if (p1->pos.y > p3->pos.y) p = p1, p1 = p3, p3 = p;
	if (p2->pos.y > p3->pos.y) p = p2, p2 = p3, p3 = p;
	if (p1->pos.y == p2->pos.y && p1->pos.y == p3->pos.y) return 0;
	if (p1->pos.x == p2->pos.x && p1->pos.x == p3->pos.x) return 0;

	if (p1->pos.y == p2->pos.y) {	// triangle down
		if (p1->pos.x > p2->pos.x) p = p1, p1 = p2, p2 = p;
		trap[0].top = p1->pos.y;
		trap[0].bottom = p3->pos.y;
		trap[0].left.v1 = *p1;
		trap[0].left.v2 = *p3;
		trap[0].right.v1 = *p2;
		trap[0].right.v2 = *p3;
		return (trap[0].top < trap[0].bottom)? 1 : 0;
	}

	if (p2->pos.y == p3->pos.y) {	// triangle up
		if (p2->pos.x > p3->pos.x) p = p2, p2 = p3, p3 = p;
		trap[0].top = p1->pos.y;
		trap[0].bottom = p3->pos.y;
		trap[0].left.v1 = *p1;
		trap[0].left.v2 = *p2;
		trap[0].right.v1 = *p1;
		trap[0].right.v2 = *p3;
		return (trap[0].top < trap[0].bottom)? 1 : 0;
	}

	trap[0].top = p1->pos.y;
	trap[0].bottom = p2->pos.y;
	trap[1].top = p2->pos.y;
	trap[1].bottom = p3->pos.y;

	k = (p3->pos.y - p1->pos.y) / (p2->pos.y - p1->pos.y);
	x = p1->pos.x + (p2->pos.x - p1->pos.x) * k;

	if (x <= p3->pos.x) {		// triangle left
		trap[0].left.v1 = *p1;
		trap[0].left.v2 = *p2;
		trap[0].right.v1 = *p1;
		trap[0].right.v2 = *p3;
		trap[1].left.v1 = *p2;
		trap[1].left.v2 = *p3;
		trap[1].right.v1 = *p1;
		trap[1].right.v2 = *p3;
	}	else {					// triangle right
		trap[0].left.v1 = *p1;
		trap[0].left.v2 = *p3;
		trap[0].right.v1 = *p1;
		trap[0].right.v2 = *p2;
		trap[1].left.v1 = *p1;
		trap[1].left.v2 = *p3;
		trap[1].right.v1 = *p2;
		trap[1].right.v2 = *p3;
	}

	return 2;
}

// 按照 Y 坐标计算出左右两条边纵坐标等于 Y 的顶点
void trapezoid_edge_interp(trapezoid_t *trap, float y) {
	float s1 = trap->left.v2.pos.y - trap->left.v1.pos.y;
	float s2 = trap->right.v2.pos.y - trap->right.v1.pos.y;
	float t1 = (y - trap->left.v1.pos.y) / s1;
	float t2 = (y - trap->right.v1.pos.y) / s2;
	vertex_interp(&trap->left.v, &trap->left.v1, &trap->left.v2, t1);
	vertex_interp(&trap->right.v, &trap->right.v1, &trap->right.v2, t2);
}

// 根据左右两边的端点，初始化计算出扫描线的起点和步长
void trapezoid_init_scan_line(const trapezoid_t *trap, scanline_t *scanline, int y) {
	float width = trap->right.v.pos.x - trap->left.v.pos.x;
	scanline->x = (int)(trap->left.v.pos.x + 0.5f);
	scanline->w = (int)(trap->right.v.pos.x + 0.5f) - scanline->x;
    if (trap->left.v.pos.x >= trap->right.v.pos.x)
        scanline->w = 0;
    
    scanline->y = y;
    scanline->v = trap->left.v;

    vertex_division(&scanline->step, &trap->left.v, &trap->right.v, width);
}

float *pshadowbuffer = NULL;

void device_init(device_t *device)
{
    device->background = 0xFFFFFF;
    device->foreground = 0;
    device->render_state = RENDER_STATE_TEXTURE;
    matrix_set_identity(&device->transform.model);
    matrix_set_identity(&device->transform.view);
}

void device_set_background(device_t *device, IUINT32 color)
{
    device->background = color;
}

void device_set_framebuffer(device_t *device, IUINT32 *framebuffer)
{
    device->framebuffer = framebuffer;
}

void device_set_zbuffer(device_t *device, float *zbuffer)
{
    device->zbuffer = zbuffer;
}

void device_set_shadowbuffer(device_t *device, float *shadowbuffer)
{
    device->shadowbuffer = shadowbuffer;
}

void device_set_camera(device_t *device, camera_t *camera)
{
    device->camera = camera;
    device->transform.view = camera->view_matrix;
    device->transform.view_r = camera->view_matrix_r;
    device->transform.projection = camera->projection_matrix;
}

void device_pixel(device_t *device, int x, int y, IUINT32 color) {
    if (x >= 0 && x < device->camera->width && y >= 0 && y < device->camera->height) {
        device->framebuffer[y * device->camera->width + x] = color;
    }
}

void device_clear(device_t *device) {
    if(device->framebuffer != NULL) {
        for(int y = 0; y < device->camera->height; y++)
            for(int x = 0; x < device->camera->width; x++)
                device->framebuffer[y * device->camera->width + x] = device->background;
    }
    // memset(device->framebuffer, 0xff, device->camera->width * device->camera->height * sizeof(IUINT32));
    if(device->zbuffer != NULL)
        memset(device->zbuffer, 0, device->camera->width * device->camera->height * sizeof(float));
    if(device->shadowbuffer != NULL) {
        for(int y = 0; y < device->camera->height; y++)
            for(int x = 0; x < device->camera->width; x++)
                device->shadowbuffer[y * device->camera->width + x] = 1.0f;
    }
}


///*
//    m = dy / dx
//    如果(ε + m) < 0.5 (或表示为2*(ε + m) < 1)
//        ε = ε + m
//    其他情况
//        ε = ε + m – 1
//将上述公式都乘以dx, 并将ε*dx用新符号ξ表示，可得
//    如果2*(ξ + dy) < dx
//        ξ = ξ + dy
//    其他情况
//        ξ = ξ + dy – dx
//*/
void device_draw_line(device_t *device, int x1, int y1, int x2, int y2, IUINT32 c) {
    int dx = x2 - x1;
    int dy = y2 - y1;
    int ux = ((dx > 0) << 1) - 1;
    int uy = ((dy > 0) << 1) - 1;
    int x = x1, y = y1, eps;

    eps = 0; dx = abs(dx); dy = abs(dy);
    if(dx > dy)
    {
        for(x = x1; x != x2 + ux; x += ux)
        {
            device_pixel(device, x, y, c);
            eps += dy;
            if((eps << 1) >= dx) {
                y += uy;
                eps -= dx;
            }
        }
    } else {
        for(y = y1; y != y2 + uy; y += uy)
        {
            device_pixel(device, x, y, c);
            eps += dx;
            if((eps << 1) >= dy) {
                x += ux;
                eps -= dy;
            }
        }
    }
}

IUINT32 texture_value_read(const texture_t *texture, float u, float v) {
    IUINT32* data = texture->datas[0];
    int width, height;
    width = texture->width;
    height = texture->height;
    u = (u-(int)u) * (width-1);
    v = (v-(int)v) * (height-1);
    int uint = (int)u;
    int vint = (int)v;
    IUINT32 res = data[vint*width+uint];
    return res;
}

// 双线性插值和mipmap
color_t texture_read(const texture_t *texture, float u, float v, float z, float maxz) {
    color_t color;
    IUINT32* data = texture->datas[0];
    int width, height;
    width = texture->width;
    height = texture->height;
    if(texture->use_mipmap) {
        int tmiplevels = logbase2ofx(width);
        int miplevel = tmiplevels * (z / maxz);
        if(miplevel > tmiplevels) miplevel = tmiplevels;
        data = (IUINT32*)texture->datas[miplevel];
        for(int ts = 0; ts < miplevel; ts++) {
            width = width >> 1;
            height = height >> 1;
        }
    }
    // wrap 方式
    u = (u-(int)u) * (width-1);
    v = (v-(int)v) * (height-1);
    int uint = (int)u;
    int vint = (int)v;
    int uint_pls_1 = uint+1;
    int vint_pls_1 = vint+1;
    uint_pls_1 = CMID(uint_pls_1, 0, width - 1);
    vint_pls_1 = CMID(vint_pls_1, 0, height - 1);
    
    int textel00, textel10, textel01, textel11;

    textel00 = data[(vint+0)*width + (uint+0)];
    textel10 = data[(vint_pls_1)*width + (uint+0)];
    textel01 = data[(vint+0)*width + (uint_pls_1)];
    textel11 = data[(vint_pls_1)*width + (uint_pls_1)];

    int textel00_a = (textel00 >> 24) & 0xff;
    int textel00_r = (textel00 >> 16) & 0xff;
    int textel00_g = (textel00 >> 8) & 0xff;
    int textel00_b = textel00 & 0xff;
    
    int textel10_a = (textel10 >> 24) & 0xff;
    int textel10_r = (textel10 >> 16) & 0xff;
    int textel10_g = (textel10 >> 8) & 0xff;
    int textel10_b = textel10 & 0xff;
    
    int textel01_a = (textel01 >> 24) & 0xff;
    int textel01_r = (textel01 >> 16) & 0xff;
    int textel01_g = (textel01 >> 8) & 0xff;
    int textel01_b = textel01 & 0xff;
    
    int textel11_a = (textel11 >> 24) & 0xff;
    int textel11_r = (textel11 >> 16) & 0xff;
    int textel11_g = (textel11 >> 8) & 0xff;
    int textel11_b = textel11 & 0xff;
    
    float dtu = u - (float)uint;
    float dtv = v - (float)vint;
    
    float one_minus_dtu = 1.0 - dtu;
    float one_minus_dtv = 1.0 - dtv;

    float one_minus_dtu_x_one_minus_dtv = (one_minus_dtu) * (one_minus_dtv);
    float dtu_x_one_minus_dtv = (dtu) * (one_minus_dtv);
    float dtu_x_dtv = (dtu) * (dtv);
    float one_minus_dtu_x_dtv = (one_minus_dtu) * (dtv);

    color.a = one_minus_dtu_x_one_minus_dtv * textel00_a +
        dtu_x_one_minus_dtv * textel01_a +
        dtu_x_dtv * textel11_a +
        one_minus_dtu_x_dtv * textel10_a;
    
    color.r = one_minus_dtu_x_one_minus_dtv * textel00_r +
        dtu_x_one_minus_dtv * textel01_r +
        dtu_x_dtv * textel11_r + 
        one_minus_dtu_x_dtv * textel10_r;

    color.g = one_minus_dtu_x_one_minus_dtv * textel00_g + 
        dtu_x_one_minus_dtv * textel01_g +
        dtu_x_dtv * textel11_g + 
        one_minus_dtu_x_dtv * textel10_g;

    color.b = one_minus_dtu_x_one_minus_dtv * textel00_b + 
        dtu_x_one_minus_dtv * textel01_b +
        dtu_x_dtv * textel11_b + 
        one_minus_dtu_x_dtv * textel10_b;
    color_scale(&color, 1.0f/255);
    return color;
}

bool computeBarycentricCoords3d(point_t *res, const point_t *p0, const point_t *p1,const point_t *p2,const point_t *p) {
    vector_t d1, d2, n;
    vector_sub(&d1, p1, p0);
    vector_sub(&d2, p2, p1);
    vector_crossproduct(&n, &d1, &d2);
    float u1, u2, u3, u4;
    float v1, v2, v3, v4;
    if((fabs(n.x) >= fabs(n.y)) && (fabs(n.x) >= fabs(n.z))) {
        u1 = p0->y - p2->y;
        u2 = p1->y - p2->y;
        u3 = p->y - p0->y;
        u4 = p->y - p2->y;
        v1 = p0->z - p2->z;
        v2 = p1->z - p2->z;
        v3 = p->z - p0->z;
        v4 = p->z - p2->z;
    } else if(fabs(n.y) >= fabs(n.z)) {
        u1 = p0->z - p2->z;
        u2 = p1->z - p2->z;
        u3 = p->z - p0->z;
        u4 = p->z - p2->z;
        v1 = p0->x - p2->x;
        v2 = p1->x - p2->x;
        v3 = p->x - p0->x;
        v4 = p->x - p2->x;
    } else {
        u1 = p0->x - p2->x;
        u2 = p1->x - p2->x;
        u3 = p->x - p0->x;
        u4 = p->x - p2->x;
        v1 = p0->y - p2->y;
        v2 = p1->y - p2->y;
        v3 = p->y - p0->y;
        v4 = p->y - p2->y;
    }
    
    float denom = v1 * u2 - v2 * u1;
    if (fabsf(denom) < 1e-6) {
        return false;
    }
    float oneOverDenom = 1.0f / denom;
    res->x = (v4 * u2 - v2 * u4) * oneOverDenom;
    res->y = (v1 * u3 - v3 * u1) * oneOverDenom;
    res->z = 1.0f - res->x - res->y;
    return true;
}

// now is the core realize for rendering, so just do it.
void device_draw_scanline(device_t *device, scanline_t *scanline, point_t *points, v2f *vfs) {
    int y = scanline->y;
    int x = scanline->x;
    int width = device->camera->width;
    int render_state = device->render_state;
    int count = scanline->w;
    for(; count > 0 && x < width; x++, count--) {
        if(x >= 0) {
            if(device->shadowbuffer != NULL) {
                float z = scanline->v.pos.z;
                if(z <= device->shadowbuffer[y*width+x]) {
                    device->shadowbuffer[y*width+x] = z;
                }
            }
            
            float rhw = scanline->v.pos.w;
            if(device->zbuffer == NULL || rhw >= device->zbuffer[y*width+x]) {
                if(device->zbuffer != NULL)
                    device->zbuffer[y*width+x] = rhw;
                
                if(device->framebuffer != NULL) {
                    color_t color = {0.0f, 0.0f, 0.0f, 1.0f};
                    v2f vf;
                    float w = 1.0f / scanline->v.pos.w;
                    point_t barycenter = {0.0f, 0.0f, 0.0f, 1.0f};
                    point_t interpos = scanline->v.pos;
                    transform_homogenize_reverse(&interpos, &interpos, w, device->camera->width, device->camera->height);
                    computeBarycentricCoords3d(&barycenter, &points[0], &points[1], &points[2], &interpos);
                    
                    
                    v2f_interpolating(&vf, &vfs[0], &vfs[1], &vfs[2], barycenter.x, barycenter.y, barycenter.z);
                    vf.pos.w = w;
                    vector_normalize(&vf.normal);
                    
                    frag_shader(device, &vf, &color);

                    float a = 1.0f;
                    float r = 0.0f;
                    float g = 0.0f;
                    float b = 0.0f;
                    
                    if(render_state & RENDER_STATE_COLOR) {
                        a = vf.color.a;
                        r = vf.color.r;
                        g = vf.color.g;
                        b = vf.color.b;
                    }
                    if(render_state & RENDER_STATE_TEXTURE) {
                        a = color.a;
                        r = color.r;
                        g = color.g;
                        b = color.b;
                    }
                    
                    int A = CMID((int)(a * 255.0f), 0, 255);
                    int R = CMID((int)(r * 255.0f), 0, 255);
                    int G = CMID((int)(g * 255.0f), 0, 255);
                    int B = CMID((int)(b * 255.0f), 0, 255);
                    
                    device->framebuffer[y*width+x] = (R << 16) | (G << 8) | B;
                }
            }
        }
        vertex_add(&scanline->v, &scanline->step);
    }
}

// core render function
void device_render_trap(device_t *device, trapezoid_t *trap, point_t *points, v2f *v2fs) {
    scanline_t scanline;
    int j, top, bottom;
    top = (int)(trap->top + 0.5f);
    bottom = (int)(trap->bottom + 0.5f);
    for(j = top; j < bottom; j++) {
        if(j >= 0 && j < device->camera->height) {
            trapezoid_edge_interp(trap, (float)j + 0.5f);
            trapezoid_init_scan_line(trap, &scanline, j);
            device_draw_scanline(device, &scanline, points, v2fs);
        }
        if(j >= device->camera->height)
            break;
    }
}
/*
void CalculateTangentArray(vector_t *tangent, vector_t *binormal, const vector_t *vec1, const vector_t *vec2, const vector_t *vec3, float u1, float v1, float u2, float v2, float u3, float v3)
{
    float x1 = vec2->x - vec1->x;
    float x2 = vec3->x - vec1->x;
    float y1 = vec2->y - vec1->y;
    float y2 = vec3->y - vec1->y;
    float z1 = vec2->z - vec1->z;
    float z2 = vec3->z - vec1->z;
    float s1 = u2 - u1;
    float s2 = u3 - u1;
    float t1 = v2 - v1;
    float t2 = v3 - v1;
    float r = 1.0f / (s1 * t2 - s2 * t1);
    vector_t sdir = {(t2 * x1 - t1 * x2) * r, (t2 * y1 - t1 * y2) * r, (t2 * z1 - t1 * z2) * r};
    vector_t tdir = {(s1 * x2 - s2 * x1) * r, (s1 * y2 - s2 * y1) * r, (s1 * z2 - s2 * z1) * r};

    const Vector3D& n = normal[a]; const Vector3D& t = tan1[a];
    // Gram-Schmidt orthogonalize.
    tangent[a] = (t - n * Dot(n, t)).Normalize();
    // Calculate handedness.
    tangent[a].w = (Dot(Cross(n, t), tan2[a]) < 0.0F) ? -1.0F : 1.0F;
}
*/
// http://www.cnblogs.com/ThreeThousandBigWorld/archive/2012/07/16/2593892.html
// http://blog.chinaunix.net/uid-26651460-id-3083223.html
// http://stackoverflow.com/questions/5255806/how-to-calculate-tangent-and-binormal
void calculate_tangent_and_binormal(vector_t *tangent, vector_t *binormal, const vector_t *position1, const vector_t *position2, const vector_t *position3,
                                          float u1, float v1, float u2, float v2, float u3, float v3)
{
    //side0 is the vector along one side of the triangle of vertices passed in,
    //and side1 is the vector along another side. Taking the cross product of these returns the normal.
    vector_t side0 = {0.0f, 0.0f, 0.0f, 1.0f};
    vector_sub(&side0, position1, position2);
    vector_t side1 = {0.0f, 0.0f, 0.0f, 1.0f};
    vector_sub(&side1, position3, position1);
    //Calculate face normal
    vector_t normal = {0.0f, 0.0f, 0.0f, 0.0f};
    vector_crossproduct(&normal, &side1, &side0);
    vector_normalize(&normal);
    //Now we use a formula to calculate the tangent.
    float deltaV0 = v1 - v2;
    float deltaV1 = v3 - v1;
    *tangent = side0;
    vector_scale(tangent, deltaV1);
    vector_t temp = side1;
    vector_scale(&temp, deltaV0);
    vector_sub(tangent, tangent, &temp);
    vector_normalize(tangent);
    //Calculate binormal
    float deltaU0 = u1 - u2;
    float deltaU1 = u3 - u1;
    *binormal = side0;
    vector_scale(binormal, deltaU1);
    temp = side1;
    vector_scale(&temp, deltaU0);
    vector_sub(binormal, binormal, &temp);
    vector_normalize(binormal);
    //Now, we take the cross product of the tangents to get a vector which
    //should point in the same direction as our normal calculated above.
    //If it points in the opposite direction (the dot product between the normals is less than zero),
    //then we need to reverse the s and t tangents.
    //This is because the triangle has been mirrored when going from tangent space to object space.
    //reverse tangents if necessary
    vector_t tangentCross;
    vector_crossproduct(&tangentCross, tangent, binormal);
    if (vector_dotproduct(&tangentCross, &normal) < 0.0f)
    {
        tangent->w = -1;
        //vector_inverse(tangent);
        //vector_inverse(binormal);
    }
    
    //binormal->w = 0.0f;
}

void device_draw_primitive(device_t *device, vertex_t *t1, vertex_t *t2, vertex_t *t3) {
    vertex_t *vertice[3] = {t1, t2, t3};
    point_t points[3];
    
    // 使用法向量背面剔除
    //    matrix_apply(&c1, &t1->pos, &device->transform.world);
    //    matrix_apply(&c2, &t2->pos, &device->transform.world);
    //    matrix_apply(&c3, &t3->pos, &device->transform.world);
    
    //    matrix_apply(&p1, &c1, &device->transform.view);
    
    //    matrix_apply(&p2, &c2, &device->transform.view);
    //    matrix_apply(&p3, &c3, &device->transform.view);
    
    //    vector_t t21, t32;
    //    vector_sub(&t21, &p2, &p1);
    //    vector_sub(&t32, &p3, &p2);
    //    vector_t normal;
    //    vector_crossproduct(&normal, &t21, &t32);
    //    if(vector_dotproduct(&normal, &p2) >= 0)
    //        return;
    
    // 正规矩阵
    matrix_t nm;
    matrix_clone(&nm, &device->transform.model);
    matrix_inverse(&nm);
    matrix_transpose(&nm);
    
//    简单的cvv裁剪，整个三角形裁减掉
//    if (transform_check_cvv(&c1) != 0) return;
//    if (transform_check_cvv(&c2) != 0) return;
//    if (transform_check_cvv(&c3) != 0) return;
    
    a2v a2vs[3];
    v2f v2fs[3];
    for(int i = 0; i < 3; i++) {
        vertex_t *vertex = vertice[i];
        a2v *av = &a2vs[i];
        
        av->pos = vertex->pos; // 世界空间的pos
        int a = 0, b = 0;
        if(i == 0) a = 1, b = 2;
        if(i == 1) a = 0, b = 2;
        if(i == 2) a = 0, b = 1;
        calculate_tangent_and_binormal(&av->tangent, &av->binormal, &vertex->pos, &vertice[a]->pos, &vertice[b]->pos, vertex->tc.u, vertex->tc.v, vertice[a]->tc.u, vertice[a]->tc.v, vertice[b]->tc.u, vertice[b]->tc.v);
  
        matrix_apply(&av->tangent, &av->tangent, &device->transform.model);
        vector_crossproduct(&av->binormal, &av->normal, &av->tangent);
        vector_scale(&av->binormal, av->tangent.w);
        //matrix_apply(&av->binormal, &av->binormal, &device->transform.model);
        matrix_apply(&vertex->pos, &vertex->pos, &device->transform.vp);
        points[i] = vertex->pos; // 透视空间的pos
        
        matrix_apply(&vertex->normal, &vertex->normal, &nm); // 法向量乘正规矩阵
        vector_normalize(&vertex->normal);
        av->normal = vertex->normal; // 世界空间的normal
        av->color = vertex->color;
        av->texcoord = vertex->tc;
    
        vert_shader(device, av, &v2fs[i]); // 顶点着色器
//        vertex_rhw_init(vertex);
        
        transform_homogenize(&vertex->pos, &vertex->pos, device->camera->width, device->camera->height);
    }
    
    // 背面剔除
    if(device->cull > 0) {
        vector_t t21, t32;
        vector_sub(&t21, &t2->pos, &t1->pos);
        vector_sub(&t32, &t3->pos, &t2->pos);
        if(device->cull == 1) {
            if(t21.x * t32.y - t32.x * t21.y <= 0)    // 计算叉积
                return;
        }
        else if(device->cull == 2) {
            if(t21.x * t32.y - t32.x * t21.y > 0)     // 计算叉积
                return;
        }
    }
    
    if (device->render_state & (RENDER_STATE_TEXTURE | RENDER_STATE_COLOR)) {
        trapezoid_t traps[2];
        int n = trapezoid_init_triangle(traps, t1, t2, t3);
        if(n >= 1) device_render_trap(device, &traps[0], points, v2fs);
        if(n >= 2) device_render_trap(device, &traps[1], points, v2fs);
    }
    
    if((device->render_state & RENDER_STATE_WIREFRAME) && device->framebuffer != NULL) {
        device_draw_line(device, (int)t1->pos.x, (int)t1->pos.y, (int)t2->pos.x, (int)t2->pos.y, device->foreground);
        device_draw_line(device, (int)t1->pos.x, (int)t1->pos.y, (int)t3->pos.x, (int)t3->pos.y, device->foreground);
        device_draw_line(device, (int)t3->pos.x, (int)t3->pos.y, (int)t2->pos.x, (int)t2->pos.y, device->foreground);
    }
}

void clip_polys(device_t *device, vertex_t *v1, vertex_t *v2, vertex_t *v3, bool world) {
    #define CLIP_CODE_GZ    0x0001
    #define CLIP_CODE_LZ    0x0002
    #define CLIP_CODE_IZ    0x0004

    #define CLIP_CODE_GX    0x0001
    #define CLIP_CODE_LX    0x0002
    #define CLIP_CODE_IX    0x0004

    #define CLIP_CODE_GY    0x0001
    #define CLIP_CODE_LY    0x0002
    #define CLIP_CODE_IY    0x0004

    #define CLIP_CODE_NULL  0x0000
    
    int vertex_ccodes[3];
    int num_verts_in = 0;

    float z_factor_x, z_factor_y, z_factor, z_test;
    
    float xi, yi, x01i, y01i, x02i, y02i, t1, t2, ui, vi, u01i, v01i, u02i, v02i;
    
    bool cliped = false;

    vector_t v;
    
    vertex_t p1 = *v1, p2 = *v2, p3 = *v3;
    
    if(world == false) {
        matrix_apply(&p1.pos, &p1.pos, &device->transform.mv);
        matrix_apply(&p2.pos, &p2.pos, &device->transform.mv);
        matrix_apply(&p3.pos, &p3.pos, &device->transform.mv);
    } else {
        matrix_apply(&p1.pos, &p1.pos, &device->transform.view);
        matrix_apply(&p2.pos, &p2.pos, &device->transform.view);
        matrix_apply(&p3.pos, &p3.pos, &device->transform.view);
    }
    
    
    z_factor_y = tan(device->camera->fovy*0.5);
    z_factor_x = z_factor_y / device->camera->aspect;
    z_factor = z_factor_x;
    z_test = z_factor * p1.pos.z;

    if (p1.pos.x > z_test)
        vertex_ccodes[0] = CLIP_CODE_GX;
    else if (p1.pos.x < -z_test)
        vertex_ccodes[0] = CLIP_CODE_LX;
    else
        vertex_ccodes[0] = CLIP_CODE_IX;

    z_test = z_factor * p2.pos.z;

    if (p2.pos.x > z_test)
        vertex_ccodes[1] = CLIP_CODE_GX;
    else if (p2.pos.x < -z_test)
        vertex_ccodes[1] = CLIP_CODE_LX;
    else
        vertex_ccodes[1] = CLIP_CODE_IX;

    z_test = z_factor * p3.pos.z;

    if (p3.pos.x > z_test)
        vertex_ccodes[2] = CLIP_CODE_GX;
    else if (p3.pos.x < -z_test)
        vertex_ccodes[2] = CLIP_CODE_LX;
    else
        vertex_ccodes[2] = CLIP_CODE_IX;

    if (((vertex_ccodes[0] == CLIP_CODE_GX) && (vertex_ccodes[1] == CLIP_CODE_GX) && (vertex_ccodes[2] == CLIP_CODE_GX))
     || ((vertex_ccodes[0] == CLIP_CODE_LX) && (vertex_ccodes[1] == CLIP_CODE_LX) && (vertex_ccodes[2] == CLIP_CODE_LX))) {
        return;
    }

    z_factor = z_factor_y;
    z_test = z_factor * p1.pos.z;

    if (p1.pos.y > z_test)
        vertex_ccodes[0] = CLIP_CODE_GY;
    else if (p1.pos.y < -z_test)
        vertex_ccodes[0] = CLIP_CODE_LY;
    else
        vertex_ccodes[0] = CLIP_CODE_IY;

    z_test = z_factor * p2.pos.z;

    if (p2.pos.y > z_test)
        vertex_ccodes[1] = CLIP_CODE_GY;
    else if (p2.pos.y < -z_test)
        vertex_ccodes[1] = CLIP_CODE_LY;
    else
        vertex_ccodes[1] = CLIP_CODE_IY;

    z_test = z_factor * p3.pos.z;

    if (p3.pos.y > z_test)
        vertex_ccodes[2] = CLIP_CODE_GY;
    else if (p3.pos.y < -z_test)
        vertex_ccodes[2] = CLIP_CODE_LY;
    else
        vertex_ccodes[2] = CLIP_CODE_IY;

    if (((vertex_ccodes[0] == CLIP_CODE_GY) && (vertex_ccodes[1] == CLIP_CODE_GY) && (vertex_ccodes[2] == CLIP_CODE_GY))
     || ((vertex_ccodes[0] == CLIP_CODE_LY) && (vertex_ccodes[1] == CLIP_CODE_LY) && (vertex_ccodes[2] == CLIP_CODE_LY))) {
        return;
    }

    // 是否完全在远裁剪面或近裁剪面外侧
    if (p1.pos.z > device->camera->zf)
        vertex_ccodes[0] = CLIP_CODE_GZ;
    else if (p1.pos.z < device->camera->zn)
        vertex_ccodes[0] = CLIP_CODE_LZ;
    else {
        vertex_ccodes[0] = CLIP_CODE_IZ;
        num_verts_in++;
    }

    if (p2.pos.z > device->camera->zf)
        vertex_ccodes[1] = CLIP_CODE_GZ;
    else if (p2.pos.z < device->camera->zn)
        vertex_ccodes[1] = CLIP_CODE_LZ;
    else {
        vertex_ccodes[1] = CLIP_CODE_IZ;
        num_verts_in++;
    }
    
        
    if (p3.pos.z > device->camera->zf)
        vertex_ccodes[2] = CLIP_CODE_GZ;
    else if (p3.pos.z < device->camera->zn)
        vertex_ccodes[2] = CLIP_CODE_LZ;
    else {
        vertex_ccodes[2] = CLIP_CODE_IZ;
        num_verts_in++;
    }
    

    if (((vertex_ccodes[0] == CLIP_CODE_GZ) && (vertex_ccodes[1] == CLIP_CODE_GZ) && (vertex_ccodes[2] == CLIP_CODE_GZ))
     || ((vertex_ccodes[0] == CLIP_CODE_LZ) && (vertex_ccodes[1] == CLIP_CODE_LZ) && (vertex_ccodes[2] == CLIP_CODE_LZ))) {
        return;
    }

    // 判断是否有顶点在近裁剪面外侧
    if (((vertex_ccodes[0] | vertex_ccodes[1] | vertex_ccodes[2]) & CLIP_CODE_LZ))
    {
        vertex_t temp;
        //num_verts_in = 0;
        // 三角形有1个顶点在近裁剪面内侧，2个顶点在外侧
        // 三角形有2个顶点在近裁剪面内侧，1个顶点在外侧
        if (num_verts_in == 1)
        {
            if (vertex_ccodes[0] == CLIP_CODE_IZ) {
            } else if (vertex_ccodes[1] == CLIP_CODE_IZ) {
                temp = p1;
                p1 = p2;
                p2 = p3;
                p3 = temp;
            } else {
                temp = p1;
                p1 = p3;
                p3 = p2;
                p2 = temp;
            }
            // 对每条边进行裁剪
            // 创建参数化方程p = v0 + v01 * t
            vector_sub(&v, &p2.pos, &p1.pos);
            t1 = (device->camera->zn - p1.pos.z) / v.z;

            xi = p1.pos.x + v.x * t1;
            yi = p1.pos.y + v.y * t1;

            // 用交点覆盖原来的顶点
            p2.pos.x = xi;
            p2.pos.y = yi;
            p2.pos.z = device->camera->zn;

            // 对三角形边v0->v2进行裁剪
            vector_sub(&v, &p3.pos, &p1.pos);
            t2 = (device->camera->zn - p1.pos.z) / v.z;

            xi = p1.pos.x + v.x * t2;
            yi = p1.pos.y + v.y * t2;

            // 用交点覆盖原来的顶点
            p3.pos.x = xi;
            p3.pos.y = yi;
            p3.pos.z = device->camera->zn;

            // 对纹理进行裁剪
            ui = p1.tc.u + (p2.tc.u - p1.tc.u) * t1;
            vi = p1.tc.v + (p2.tc.v - p1.tc.v) * t1;

            p2.tc.u = ui;
            p2.tc.v = vi;

            ui = p1.tc.u + (p3.tc.u - p1.tc.u) * t2;
            vi = p1.tc.v + (p3.tc.v - p1.tc.v) * t2;

            p3.tc.u = ui;
            p3.tc.v = vi;
            
            cliped = true;
        } else if (num_verts_in == 2) {
            // 外侧的点
            if (vertex_ccodes[0] == CLIP_CODE_LZ) {
            } else if (vertex_ccodes[1] == CLIP_CODE_LZ) {
                temp = p1;
                p1 = p2;
                p2 = p3;
                p3 = temp;
            } else {
                temp = p1;
                p1 = p3;
                p3 = p2;
                p2 = temp;
            }
            
            vertex_t np1 = p1, np2 = p2, np3 = p3;
            // 对每条边进行裁剪
            // 创建参数化方程p = v0 + v01 * t
            vector_sub(&v, &p2.pos, &p1.pos);
            t1 = (device->camera->zn - p1.pos.z) / v.z;

            x01i = p1.pos.x + v.x * t1;
            y01i = p1.pos.y + v.y * t1;

            // 对三角形边v0->v2进行裁剪
            vector_sub(&v, &p3.pos, &p1.pos);
            t2 = (device->camera->zn - p1.pos.z) / v.z;

            x02i = p1.pos.x + v.x * t2;
            y02i = p1.pos.y + v.y * t2;

            // 用交点覆盖原来的顶点
            p1.pos.x = x01i;
            p1.pos.y = y01i;
            p1.pos.z = device->camera->zn;
            
            np2.pos.x = x01i;
            np2.pos.y = y01i;
            np2.pos.z = device->camera->zn;

            np1.pos.x = x02i;
            np1.pos.y = y02i;
            np1.pos.z = device->camera->zn;

            // 对纹理进行裁剪
            u01i = p1.tc.u + (p2.tc.u - p1.tc.u) * t1;
            v01i = p1.tc.v + (p2.tc.v - p1.tc.v) * t1;

            u02i = p1.tc.u + (p3.tc.u - p1.tc.u) * t2;
            v02i = p1.tc.v + (p3.tc.v - p1.tc.v) * t2;

            p1.tc.u = u01i;
            p1.tc.v = v01i;

            np1.tc.u = u02i;
            np1.tc.v = v02i;
            np2.tc.u = u01i;
            np2.tc.v = v01i;
            
            matrix_apply(&np1.pos, &np1.pos, &device->transform.view_r);
            matrix_apply(&np2.pos, &np2.pos, &device->transform.view_r);
            matrix_apply(&np3.pos, &np3.pos, &device->transform.view_r);
            device_draw_primitive(device, &np1, &np2, &np3);
            
            cliped = true;
        }
    }
    matrix_apply(&p1.pos, &p1.pos, &device->transform.view_r);
    matrix_apply(&p2.pos, &p2.pos, &device->transform.view_r);
    matrix_apply(&p3.pos, &p3.pos, &device->transform.view_r);
    device_draw_primitive(device, &p1, &p2, &p3);
}

void vert_shader(device_t *device, a2v *av, v2f *vf) {
    vf->pos = av->pos;
    vf->normal = av->normal;
    vf->color = av->color;
    vf->texcoord = av->texcoord;
    
    vf->storage0 = (vector_t){av->tangent.x, av->binormal.x, av->normal.x};
    vf->storage1 = (vector_t){av->tangent.y, av->binormal.y, av->normal.y};
    vf->storage2 = (vector_t){av->tangent.z, av->binormal.z, av->normal.z};
    
    //vf->storage0 = av->tangent;
    //vf->storage1 = av->binormal;
    //vf->storage2 = av->normal;
}

void frag_shader(device_t *device, v2f *vf, color_t *color) {
    material_t *material = &device->material;
    vector_t viewdir, viewPos = device->camera->pos;
    vector_sub(&viewdir, &viewPos, &vf->pos);
    vector_normalize(&viewdir);
    texcoord_t *tex = &vf->texcoord;
    vector_t normal = vf->normal;
    
    vector_t lightDir = dirLight.dir;
    vector_inverse(&lightDir);
    vector_normalize(&lightDir);
    
//    if(material->bump_tex_id != -1) {
//        color_t color = texture_read(&textures[material->bump_tex_id], tex->u, tex->v, vf->pos.w, 15);
//        vector_t bump = {color.r, color.g, color.b, color.a};
//        bump.x = bump.x * 2 - 1.0f;
//        bump.y = bump.y * 2 - 1.0f;
//        vector_scale(&bump, 1.0f);
//        float n = bump.x * bump.x + bump.y * bump.y;
//        if( n < 0.0f ) n = 0.0f; if( n > 1.0f ) n = 1.0f;
//        bump.z = sqrtf(1.0f - n);
//        normal = (vector_t){vector_dotproduct(&vf->storage0, &bump), vector_dotproduct(&vf->storage1, &bump), vector_dotproduct(&vf->storage2, &bump)};
//        vector_normalize(&normal);
//    }
    
    float diff = fmaxf(vector_dotproduct(&normal, &lightDir), 0.0f);
    lightDir = dirLight.dir;
    vector_normalize(&lightDir);
    vector_t vec;
    vector_reflect(&vec, &lightDir, &normal);
    float shininess = material->shininess * (material->specular_highlight_tex_id == -1 ? 1 : texture_value_read(&textures[material->specular_highlight_tex_id], tex->u, tex->v));
    float spec = powf(fmaxf(vector_dotproduct(&viewdir, &vec), 0.0f), shininess);
    
    color_t temp = {0.0f, 0.0f ,0.0f, 1.0f};
    color_t temp2 = material->ambient_tex_id == -1 ? material->ambient : texture_read(&textures[material->ambient_tex_id], tex->u, tex->v, vf->pos.w, 15);
    color_product(&temp, &dirLight.ambi, &temp2);
    color_add(color, color, &temp);
    temp2 = material->diffuse_tex_id == -1 ? material->diffuse : texture_read(&textures[material->diffuse_tex_id], tex->u, tex->v, vf->pos.w, 15);
    color_product(&temp, &dirLight.diff, &temp2);
    color_scale(&temp, diff);
    color_add(color, color, &temp);
    temp2 = material->specular_tex_id == -1 ? material->specular : texture_read(&textures[material->specular_tex_id], tex->u, tex->v, vf->pos.w, 15);
    color_product(&temp, &dirLight.spec, &temp2);
    color_scale(&temp, spec);
    color_add(color, color, &temp);
    
    // 计算阴影
    if(dirLight.shadow)
    {
        // fragPos -> 灯的坐标系 -> 灯的透视矩阵 -> 求得z坐标比较
        point_t tempPos = vf->pos;
        tempPos.w = 1.0f;
        camera_t *camera = &cameras[0];
        matrix_apply(&tempPos, &tempPos, &camera->view_matrix);
        matrix_apply(&tempPos, &tempPos, &camera->projection_matrix);
        transform_homogenize(&tempPos, &tempPos, camera->width, camera->height);
        int y = (int)(tempPos.y+0.5);
        int x = (int)(tempPos.x+0.5);
        
        vector_t tempNormal = vf->normal;
        matrix_apply(&tempNormal, &tempNormal, &camera->view_matrix);
        vector_inverse(&tempNormal);
        float dot = vector_dotproduct(&tempNormal, &camera->front);
        if(dot > 0) {
            float bias = 0.015 * (1.0 - dot);
            if(bias < 0.002f) bias = 0.001;
            if(y >= 0 && x >= 0 && y < camera->height && x < camera->width) {
                float shadow = 0.0;
                for(int i = -1; i <= 1; ++i)
                {
                    for(int j = -1; j <= 1; ++j)
                    {
                        if(y+j < 0 || y+j >= camera->height || x+i < 0 || x+i >= camera->width)
                            continue;
                        float pcfDepth = pshadowbuffer[(y+j)*camera->width+(x+i)];
                        shadow +=  tempPos.z - bias > pcfDepth ? 1.0 : 0.0;
                    }
                }
                shadow /= 9.0;
                
                color_t temp = {0.3f,0.3f,0.3f,0.3f};
                color_scale(&temp, shadow);
                color_sub(color, color, &temp);
            }
        }
    }
    
    
    int i = 0;
    for(i = 0; i < pointlight_cnt; i++)
    {
        temp = (color_t){0.0f, 0.0f ,0.0f, 1.0f};
        pointlight_t *pointlight = &pointLights[i];
        
        vector_t lightDir = {0.0f, 0.0f, 0.0f, 0.0f};
        vector_sub(&lightDir, &pointlight->pos, &vf->pos);
        float distance = vector_length(&lightDir);
        vector_normalize(&lightDir);
        float diff = fmaxf(vector_dotproduct(&normal, &lightDir), 0.0f);
        vector_t vec;
        vector_inverse(&lightDir);
        vector_reflect(&vec, &lightDir, &normal);
        float shininess = material->shininess * (material->specular_highlight_tex_id == -1 ? 1 : texture_value_read(&textures[material->specular_highlight_tex_id], tex->u, tex->v));
        float spec = powf(fmaxf(vector_dotproduct(&viewdir, &vec), 0.0f), shininess);
        float num = pointlight->constant + pointlight->linear * distance + pointlight->quadratic * (distance * distance);
        float attenuation = 0;
        if(num != 0)
            attenuation = 1.0f / num;
        
        color_t c = (color_t){0.0f, 0.0f ,0.0f, 1.0f};
        color_t c2 = material->ambient_tex_id == -1 ? material->ambient : texture_read(&textures[material->ambient_tex_id], tex->u, tex->v, vf->pos.w, 15);
        color_product(&c, &pointlight->ambi, &c2);
        color_scale(&c, attenuation);
        color_add(&temp, &temp, &c);
        c2 = material->diffuse_tex_id == -1 ? material->diffuse : texture_read(&textures[material->diffuse_tex_id], tex->u, tex->v, vf->pos.w, 15);
        color_product(&c, &pointlight->diff, &c2);
        color_scale(&c, diff * attenuation);
        color_add(&temp, &temp, &c);
        c2 = material->specular_tex_id == -1 ? material->specular : texture_read(&textures[material->specular_tex_id], tex->u, tex->v, vf->pos.w, 15);
        color_product(&c, &pointlight->spec, &c2);
        color_scale(&c, spec * attenuation);
        color_add(&temp, &temp, &c);
        
        color_add(color, color, &temp);
    }
}


