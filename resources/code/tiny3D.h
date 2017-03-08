// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!
#include "math.h"

typedef unsigned int IUINT32;

// modules:
// 1. math library
//  1). vector
typedef struct {
    float x, y, z, w;
} vector_t;
typedef vector_t point_t;

int CMID(int x, int min, int max) { return (x < min) ? min : ((x > max) ? max : x); }

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

float interp(float a, float b, float t) {
    return a + (b - a) * t;
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
    c->w = 1.0f;
}
//      3)). sub
void vector_sub(vector_t *c, const vector_t *a, const vector_t *b) {
    c->x = a->x - b->x; 
    c->y = a->y - b->y; 
    c->z = a->z - b->z; 
    c->w = 1.0f;
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
    c->w = 1.0f;
}
//      6)). interp t[0, 1]
void vector_interp(vector_t *c, const vector_t *a, const vector_t *b, float t) {
    c->x = interp(a->x, b->x, t);
    c->y = interp(a->y, b->y, t);
    c->z = interp(a->z, b->z, t);
    c->w = 1.0f;
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
    if(abs(len) < 1e-6) return;
    float k = 1 / len;
    v->x *= k;
    v->y *= k;
    v->z *= k;
}

typedef struct {
    float m[4][4];
} matrix_t;

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
#define N 3
float getA(float arcs[N][N], int n)//按第一行展开计算|A|
{
    if(n==1)
    {
        return arcs[0][0];
    }
    float ans = 0;
    float temp[N][N];
    int i,j,k;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n-1;j++)
        {
            for(k=0;k<n-1;k++)
            {
                temp[j][k] = arcs[j+1][(k>=i)?k+1:k];
                
            }
        }
        float t = getA(temp,n-1);
        if(i%2==0)
        {
            ans += arcs[0][i]*t;
        }
        else
        {
            ans -= arcs[0][i]*t;
        }
    }
    return ans;
}
void getAStart(float arcs[N][N], int n, float ans[N][N])//计算每一行每一列的每个元素所对应的余子式，组成A*
{
    if(n==1)
    {
        ans[0][0] = 1;
        return;
    }
    int i,j,k,t;
    float temp[N][N];
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<n-1;k++)
            {
                for(t=0;t<n-1;t++)
                {
                    temp[k][t] = arcs[k>=i?k+1:k][t>=j?t+1:t];
                }
            }
            
            
            ans[j][i] = getA(temp,n-1);
            if((i+j)%2 == 1)
            {
                ans[j][i] = -ans[j][i];
            }
        }
    }  
}


//      5)). inverse
void matrix_inverse(matrix_t *m) {
//    float arcs[3][3];
//    float astar[3][3];
//    int i, j;
//
//    for(i = 0; i < 3; i++)
//        for(j = 0; j < 3; j++)
//            arcs[i][j] = m->m[i][j];
//    
//    float a = getA(arcs, 3);
//    if(a==0)
//    {
//        printf("can not transform!\n");
//    }
//    else
//    {
//        getAStart(arcs, 3, astar);
//        for(i = 0; i < 3; i++)
//        {
//            for(j = 0; j <3 ; j++)
//            {
//                m->m[i][j] = astar[i][j]/a;
//            }
//        }
//    }
    
    
//    matrix_t res;
//    double determinant =    +m->m[0][0] * (m->m[1][1]*m->m[2][2] - m->m[2][1]*m->m[1][2])
//    - m->m[0][1] * (m->m[1][0]*m->m[2][2] - m->m[1][2]*m->m[2][0])
//    + m->m[0][2] * (m->m[1][0]*m->m[2][1] - m->m[1][1]*m->m[2][0]);
//    double invdet = 1/determinant;
//    res.m[0][0] =  (m->m[1][1]*m->m[2][2]-m->m[2][1]*m->m[1][2])*invdet;
//    res.m[1][0] = -(m->m[0][1]*m->m[2][2]-m->m[0][2]*m->m[2][1])*invdet;
//    res.m[2][0] =  (m->m[0][1]*m->m[1][2]-m->m[0][2]*m->m[1][1])*invdet;
//    res.m[0][1] = -(m->m[1][0]*m->m[2][2]-m->m[1][2]*m->m[2][0])*invdet;
//    res.m[1][1] =  (m->m[0][0]*m->m[2][2]-m->m[0][2]*m->m[2][0])*invdet;
//    res.m[2][1] = -(m->m[0][0]*m->m[1][2]-m->m[1][0]*m->m[0][2])*invdet;
//    res.m[0][2] =  (m->m[1][0]*m->m[2][1]-m->m[2][0]*m->m[1][1])*invdet;
//    res.m[1][2] = -(m->m[0][0]*m->m[2][1]-m->m[2][0]*m->m[0][1])*invdet;
//    res.m[2][2] =  (m->m[0][0]*m->m[1][1]-m->m[1][0]*m->m[0][1])*invdet;
//    
//    int i, j;
//    for(i = 0; i < 3; i++)
//        for(j = 0; j < 3; j++)
//            m->m[i][j-3] = res.m[i][j];
    
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

// I have realized the math library. Then I should realize the core data structure of computer graphics.
// 2. transform_t
typedef struct {
    matrix_t world;
    matrix_t view;
    matrix_t projection;
    matrix_t transform_wv;
    matrix_t transform_wv_r;
    matrix_t transform;
    float w, h;
} transform_t;

//  1). transform_update (world * view * projection
void transform_update(transform_t *ts) {
    matrix_mul(&ts->transform_wv, &ts->world, &ts->view);
    matrix_mul(&ts->transform, &ts->transform_wv, &ts->projection);
    matrix_clone(&ts->transform_wv_r, &ts->transform_wv);
    matrix_inverse(&ts->transform_wv_r);
}
//  2). transform_init (ts, width, height)
void transform_init(transform_t *ts, int width, int height, float fovy, float zn, float zf) {
    float aspect = (float)width / (float)height;
    matrix_set_identity(&ts->world);
    matrix_set_identity(&ts->view);
    matrix_set_perspective(&ts->projection, fovy, aspect, zn , zf);
    ts->w = (float)width;
    ts->h = (float)height;
    transform_update(ts);
}
//  3). transform_apply
void transform_apply(const transform_t *ts, vector_t *y, const vector_t *x) {
    matrix_apply(y, x, &ts->transform);
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
void transform_homogenize(const transform_t *ts, vector_t *y, const vector_t *x) {
    float rhw = 1.0f / x->w;
    y->x = (x->x * rhw + 1.0f) * ts->w * 0.5f;
    y->y = (1.0f - x->y * rhw) * ts->h * 0.5f;
    y->z = x->z * rhw;
    y->w = 1.0f;
}

// 3. geometry calculation (vertex, scanline, border check, rect and so on)
typedef struct { float r, g, b, a; } color_t;
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

//point_t viewPos;

typedef struct { float u, v; } texcoord_t;
typedef struct { point_t pos; texcoord_t tc; color_t color; float rhw; vector_t normal; } vertex_t;

typedef struct { vertex_t v, v1, v2; } edge_t;
typedef struct { float top, bottom; edge_t left, right; } trapezoid_t;
typedef struct { vertex_t v, step; int x, y, w; } scanline_t;

typedef struct {
    color_t ambi;
    color_t diff;
    color_t spec;
    float shininess;
} material_t;
#define NUM_MATERIAL 4
material_t materials[NUM_MATERIAL];
int material_cnt;

typedef struct {
    point_t pos;
    point_t vpos;
    float constant;
    float linear;
    float quadratic;
    color_t ambi;
    color_t diff;
    color_t spec;
} pointlight_t;
#define NR_POINT_LIGHTS 100
pointlight_t pointLights[NR_POINT_LIGHTS];
int pointlight_cnt;

void calc_pointlight(color_t *color, const material_t *material, const pointlight_t *light, const vector_t *normal, const vector_t *fpos, const vector_t *viewdir) {
    vector_t lightDir;
    vector_sub(&lightDir, &light->vpos, fpos);
    float distance = vector_length(&lightDir);
    vector_normalize(&lightDir);
    float diff = fmaxf(vector_dotproduct(normal, &lightDir), 0.0f);
    vector_t vec;
    vector_inverse(&lightDir);
    vector_reflect(&vec, &lightDir, normal);
    float spec = powf(fmaxf(vector_dotproduct(viewdir, &vec), 0.0f), material->shininess);
    float temp = light->constant + light->linear * distance + light->quadratic * (distance * distance);
    float attenuation = 0;
    if(temp != 0)
     attenuation = 1.0f / temp;
    
    color_t c;
    color_product(&c, &light->ambi, &material->ambi);
    color_scale(&c, attenuation);
    color_add(color, color, &c);
    color_product(&c, &light->diff, &material->diff);
    color_scale(&c, diff * attenuation);
    color_add(color, color, &c);
    color_product(&c, &light->spec, &material->spec);
    color_scale(&c, spec * attenuation);
    color_add(color, color, &c);
}

typedef struct {
    vector_t dir;
    vector_t vdir;
    color_t ambi;
    color_t diff;
    color_t spec;
} dirlight_t;
dirlight_t dirLight;

void calc_dirlight(color_t *color, const material_t *material, const dirlight_t *light, const vector_t *normal, const vector_t *viewdir) {
    vector_t lightDir = light->vdir;
    vector_inverse(&lightDir);
    vector_normalize(&lightDir);
    float diff = fmaxf(vector_dotproduct(normal, &lightDir), 0.0f);
    lightDir = light->vdir;
    vector_normalize(&lightDir);
    vector_t vec;
    vector_reflect(&vec, &lightDir, normal);
    //vector_inverse(&vec);
    float spec = powf(fmaxf(vector_dotproduct(viewdir, &vec), 0.0f), material->shininess);
    
    color_t temp;
    color_product(&temp, &light->ambi, &material->ambi);
    color_add(color, color, &temp);
    color_product(&temp, &light->diff, &material->diff);
    color_scale(&temp, diff);
    color_add(color, color, &temp);
    color_product(&temp, &light->spec, &material->spec);
    color_scale(&temp, spec);
    color_add(color, color, &temp);
}

// 注意是除坐标意外颜色和纹理索引除以w
void vertex_rhw_init(vertex_t *v) {
    float rhw = 1.0f / v->pos.w;
    v->rhw = rhw;
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
    y->rhw = interp(x1->rhw, x2->rhw, k);
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
    y->rhw = (x2->rhw - x1->rhw) * inv;
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
    y->rhw += x->rhw;
}

// 根据三角形生成 0-2 个梯形，并且返回合法梯形的数量
int trapezoid_init_triangle(trapezoid_t *trap, const vertex_t *p1, 
	const vertex_t *p2, const vertex_t *p3) {
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
	scanline->y = y;
	scanline->v = trap->left.v;
	if (trap->left.v.pos.x >= trap->right.v.pos.x) scanline->w = 0;
	vertex_division(&scanline->step, &trap->left.v, &trap->right.v, width);
}

typedef struct {
	transform_t transform;      // 坐标变换器
	int width;                  // 窗口宽度
	int height;                 // 窗口高度
    float aspect;               // 宽高比
    material_t material;        // 当前材质
	IUINT32 **framebuffer;      // 像素缓存：framebuffer[y] 代表第 y行
	float **zbuffer;            // 深度缓存：zbuffer[y] 为第 y行指针
	IUINT32 **texture;          // 纹理
    bool use_mipmap;            // 是否开启mipmap
	int tex_width;              // 纹理宽度
	int tex_height;             // 纹理高度
	float max_u;                // 纹理最大宽度：tex_width - 1
	float max_v;                // 纹理最大高度：tex_height - 1
	int render_state;           // 渲染状态
	IUINT32 background;         // 背景颜色
	IUINT32 foreground;         // 线框颜色
    float fovy;                 // fov
    float zn;                   // 近裁剪面
    float zf;                   // 远裁剪面
    bool blend;                 // 是否开启混合
    float blend_sfactor;        // 混合源因子
    float blend_dfactor;        // 混合目标因子
}	device_t;

#define RENDER_STATE_WIREFRAME      1		// 渲染线框
#define RENDER_STATE_TEXTURE        2		// 渲染纹理
#define RENDER_STATE_COLOR          4		// 渲染颜色

// 设备初始化，fb为外部帧缓存，非 NULL 将引用外部帧缓存（每行 4字节对齐）
void device_init(device_t *device, int width, int height, float fovy, float zn, float zf, void *fb) {
	int need = sizeof(void*) * (height * 2 + 1024 + 1) + width * height * 8;
	char *ptr = (char*)malloc(need + 64);
	char *framebuf, *zbuf;
	int j;
	device->framebuffer = (IUINT32**)ptr;
    ptr += sizeof(void*) * height;
 
	device->zbuffer = (float**)ptr;
	ptr += sizeof(void*) * height;

    device->use_mipmap = false;
    device->texture = (IUINT32**)ptr;
    ptr += sizeof(IUINT32*) * 1;

    // framebuf和zbuf的内存数据直接初始化
	framebuf = (char*)ptr;
	zbuf = (char*)ptr + width * height * 4;
	ptr += width * height * 8;
	if (fb != NULL) framebuf = (char*)fb;
	for (j = 0; j < height; j++) {
		device->framebuffer[j] = (IUINT32*)(framebuf + width * 4 * j);
		device->zbuffer[j] = (float*)(zbuf + width * 4 * j);
	}
    // texture暂时初始化一下，以后还会重新赋值
	device->texture[0] = (IUINT32*)ptr;
	device->texture[1] = (IUINT32*)(ptr + 16);
	memset(device->texture[0], 0, 64);
	device->tex_width = 2;
	device->tex_height = 2;
	device->max_u = 1.0f;
	device->max_v = 1.0f;
    device->use_mipmap = false;
	device->width = width;
	device->height = height;
    device->aspect = (float)width / (float)height;
	device->background = 0xc0c0c0;
	device->foreground = 0;
    device->fovy = fovy;
    device->zn = zn;
    device->zf = zf;
	transform_init(&device->transform, width, height, fovy, zn, zf);
	device->render_state = RENDER_STATE_WIREFRAME;
}

void device_destroy(device_t *device) {
	if (device->framebuffer) 
		free(device->framebuffer);
	device->framebuffer = NULL;
	device->zbuffer = NULL;
	device->texture = NULL;
}

void device_set_texture(device_t *device, IUINT32 **texture, int w, int h, bool use_mipmap) {
    device->texture = texture;
	device->tex_width = w;
	device->tex_height = h;
	device->max_u = (float)(w - 1);
	device->max_v = (float)(h - 1);
    device->use_mipmap = use_mipmap;
}

void device_clear(device_t *device, int mode) {
	int y, x, height = device->height;
	for (y = 0; y < device->height; y++) {
		IUINT32 *dst = device->framebuffer[y];
		IUINT32 cc = (height - 1 - y) * 230 / (height - 1);
		cc = (cc << 16) | (cc << 8) | cc;
		if (mode == 0) cc = device->background;
		for (x = device->width; x > 0; dst++, x--) dst[0] = 0xffffffff;
	}
	for (y = 0; y < device->height; y++) {
		float *dst = device->zbuffer[y];
		for (x = device->width; x > 0; dst++, x--) dst[0] = 0.0f;
	}
}

void device_pixel(device_t *device, int x, int y, IUINT32 color) {
	if (((IUINT32)x) < (IUINT32)device->width && ((IUINT32)y) < (IUINT32)device->height) {
		device->framebuffer[y][x] = color;
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

//IUINT32 device_texture_read(const device_t *device, float u, float v) {
//    int x, y;
//    u = u * device->max_u;
//    v = v * device->max_v;
//    x = (int)(u + 0.5f);
//    y = (int)(v + 0.5f);
//    x = CMID(x, 0, device->tex_width - 1);
//    y = CMID(y, 0, device->tex_height - 1);
//    return device->texture[y][x];
//}

// 双线性插值和mipmap
color_t device_texture_read(const device_t *device, float u, float v, float z, float maxz) {
    color_t color;
    IUINT32* texture = device->texture[0];
    int width, height;
    width = device->tex_width;
    height = device->tex_height;
    if(device->use_mipmap) {
        int tmiplevels = logbase2ofx(device->tex_width);
        int miplevel = tmiplevels * (z / maxz);
        if(miplevel > tmiplevels) miplevel = tmiplevels;
        texture = (IUINT32*)device->texture[miplevel];
        for(int ts = 0; ts < miplevel; ts++) {
            width = width >> 1;
            height = height >> 1;
        }
    }
    u = u * (width-1);
    v = v * (height-1);
    int uint = (int)u;
    int vint = (int)v;
    int uint_pls_1 = uint+1;
    int vint_pls_1 = vint+1;
    uint_pls_1 = CMID(uint_pls_1, 0, width - 1);
    vint_pls_1 = CMID(vint_pls_1, 0, height - 1);
    
    int textel00, textel10, textel01, textel11;

    textel00 = texture[(vint+0)*width + (uint+0)];
    textel10 = texture[(vint_pls_1)*width + (uint+0)];
    textel01 = texture[(vint+0)*width + (uint_pls_1)];
    textel11 = texture[(vint_pls_1)*width + (uint_pls_1)];

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
    
//    color.r = textel00_r;
//    color.g = textel00_g;
//    color.b = textel00_b;

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
    if (abs(denom) < 1e-6) {
        return false;
    }
    float oneOverDenom = 1.0f / denom;
    res->x = (v4 * u2 - v2 * u4) * oneOverDenom;
    res->y = (v1 * u3 - v3 * u1) * oneOverDenom;
    res->z = 1.0f - res->x - res->y;
    return true;
}

// now is the core realize for rendering, so just do it.
void device_draw_scanline(device_t *device, scanline_t *scanline, const vertex_t *t1, const vertex_t *t2, const vertex_t *t3, const point_t *p1, const point_t *p2, const point_t *p3) {
    IUINT32 *framebuffer = device->framebuffer[scanline->y];
    float *zbuffer = device->zbuffer[scanline->y];
    int x = scanline->x;
    int width = device->width;
    int render_state = device->render_state;
    int count = scanline->w;
    for(; count > 0 && x < width; x++, count--) {
        if(x >= 0) {
            float rhw = scanline->v.rhw;
            if(rhw >= zbuffer[x]) {
                float w = 1.0f / rhw;
                zbuffer[x] = rhw;
                color_t color = {0.0f, 0.0f, 0.0f, 1.0f};
                
                point_t barycenter = {0.0f, 0.0f, 0.0f, 1.0f};
                computeBarycentricCoords3d(&barycenter, &t1->pos, &t2->pos, &t3->pos, &scanline->v.pos);
                
                point_t fragPos = {0.0f, 0.0f, 0.0f, 1.0f};
                point_t ptemp = *p1;
                vector_scale(&ptemp, barycenter.x);
                vector_add(&fragPos, &fragPos, &ptemp);
                ptemp = *p2;
                vector_scale(&ptemp, barycenter.y);
                vector_add(&fragPos, &fragPos, &ptemp);
                ptemp = *p3;
                vector_scale(&ptemp, barycenter.z);
                vector_add(&fragPos, &fragPos, &ptemp);
                
                vector_t normal = {0.0f, 0.0f, 0.0f, 0.0f};
                ptemp = t1->normal;
                vector_scale(&ptemp, barycenter.x);
                vector_add(&normal, &normal, &ptemp);
                ptemp = t2->normal;
                vector_scale(&ptemp, barycenter.y);
                vector_add(&normal, &normal, &ptemp);
                ptemp = t3->normal;
                vector_scale(&ptemp, barycenter.z);
                vector_add(&normal, &normal, &ptemp);
                vector_normalize(&normal);
                //printf("%f, %f, %f\n", normal.x, normal.y, normal.z);
                vector_t viewdir, viewPos = {0.0f, 0.0f ,0.0f};
                vector_sub(&viewdir, &viewPos, &fragPos);
                vector_normalize(&viewdir);
                calc_dirlight(&color, &device->material, &dirLight, &normal, &viewdir);
                color_t temp = {0.0f, 0.0f ,0.0f, 1.0f};
                int i = 0;
                for(i = 0; i < pointlight_cnt; i++)
                {
                    calc_pointlight(&temp, &device->material, &pointLights[i], &normal, &fragPos, &viewdir);
                    color_add(&color, &color, &temp);
                }
                
                float a = 1.0f;
                float r = 0.0f;
                float g = 0.0f;
                float b = 0.0f;
                
                if(render_state & RENDER_STATE_COLOR) {
                    a = scanline->v.color.a * w * 255.0f;
                    r = scanline->v.color.r * w * 255.0f;
                    g = scanline->v.color.g * w * 255.0f;
                    b = scanline->v.color.b * w * 255.0f;
                }
                if(render_state & RENDER_STATE_TEXTURE) {
                    float u = scanline->v.tc.u * w;
                    float v = scanline->v.tc.v * w;
                    color_t cc = device_texture_read(device, u, v, w, 15);
                    a = cc.a;
                    r = cc.r;
                    g = cc.g;
                    b = cc.b;
                }
                a *= color.a;
                r *= color.r;
                g *= color.g;
                b *= color.b;
                int A = CMID((int)a, 0, 255);
                int R = CMID((int)r, 0, 255);
                int G = CMID((int)g, 0, 255);
                int B = CMID((int)b, 0, 255);
                
                framebuffer[x] = (R << 16) | (G << 8) | B;
            }
        }
        vertex_add(&scanline->v, &scanline->step);
    }
}

// core render function
void device_render_trap(device_t *device, trapezoid_t *trap, const vertex_t *t1, const vertex_t *t2, const vertex_t *t3, const point_t *p1, const point_t *p2, const point_t *p3) {
    scanline_t scanline;
    int j, top, bottom;
    top = (int)(trap->top + 0.5f);
    bottom = (int)(trap->bottom + 0.5f);
    for(j = top; j < bottom; j++) {
        if(j >= 0 && j < device->height) {
            trapezoid_edge_interp(trap, (float)j + 0.5f);
            trapezoid_init_scan_line(trap, &scanline, j);
            device_draw_scanline(device, &scanline, t1, t2, t3, p1, p2, p3);
        }
        if(j >= device->height)
            break;
    }
}

void device_draw_primitive(device_t *device, vertex_t *t1, vertex_t *t2, vertex_t *t3) {
    point_t p1, p2, p3;
    int render_state = device->render_state;
    
    p1 = t1->pos;
    p2 = t2->pos;
    p3 = t3->pos;
    // 使用法向量
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
    
    matrix_apply(&t1->pos, &t1->pos, &device->transform.projection);
    matrix_apply(&t2->pos, &t2->pos, &device->transform.projection);
    matrix_apply(&t3->pos, &t3->pos, &device->transform.projection);
    
    // 法向量乘正规矩阵
//    vector_t normals[3];
    matrix_t nm;
    matrix_clone(&nm, &device->transform.transform_wv);
    matrix_inverse(&nm);
    matrix_transpose(&nm);
    matrix_apply(&t1->normal, &t1->normal, &nm);
    matrix_apply(&t2->normal, &t2->normal, &nm);
    matrix_apply(&t3->normal, &t3->normal, &nm);
    vector_normalize(&t1->normal);
    vector_normalize(&t2->normal);
    vector_normalize(&t3->normal);
    
//    if (transform_check_cvv(&c1) != 0) return;
//    if (transform_check_cvv(&c2) != 0) return;
//    if (transform_check_cvv(&c3) != 0) return;
    
    vertex_rhw_init(t1);
    vertex_rhw_init(t2);
    vertex_rhw_init(t3);
    
    transform_homogenize(&device->transform, &t1->pos, &t1->pos);
    transform_homogenize(&device->transform, &t2->pos, &t2->pos);
    transform_homogenize(&device->transform, &t3->pos, &t3->pos);
    
    // 背面剔除
    vector_t t21, t32;
    vector_sub(&t21, &t2->pos, &t1->pos);
    vector_sub(&t32, &t3->pos, &t2->pos);
    // 计算叉积
    if(t21.x * t32.y - t32.x * t21.y <= 0)
        return;
    
    if (render_state & (RENDER_STATE_TEXTURE | RENDER_STATE_COLOR)) {
//        vertex_t t[3];
//        point_t p[3];
//        t[0] = *t1, t[1] = *t2, t[2] = *t3;
        trapezoid_t traps[2];
        int n;
        
//        p[0] = p1;
//        p[1] = p2;
//        p[2] = p3;
//        t[0].pos = p1;
//        t[1].pos = p2;
//        t[2].pos = p3;
//        t[0].pos.w = c1.w;
//        t[1].pos.w = c2.w;
//        t[2].pos.w = c3.w;
        
        n = trapezoid_init_triangle(traps, t1, t2, t3);
        
//        matrix_apply(&t[0].pos, &t1->pos, &device->transform.world);
//        matrix_apply(&t[1].pos, &t2->pos, &device->transform.world);
//        matrix_apply(&t[2].pos, &t3->pos, &device->transform.world);
        
        if(n >= 1) device_render_trap(device, &traps[0], t1, t2, t3, &p1, &p2, &p3);
        if(n >= 2) device_render_trap(device, &traps[1], t1, t2, t3, &p1, &p2, &p3);
    }
    
    if(render_state & RENDER_STATE_WIREFRAME) {
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
        matrix_apply(&p1.pos, &p1.pos, &device->transform.transform_wv);
        matrix_apply(&p2.pos, &p2.pos, &device->transform.transform_wv);
        matrix_apply(&p3.pos, &p3.pos, &device->transform.transform_wv);
    } else {
        matrix_apply(&p1.pos, &p1.pos, &device->transform.view);
        matrix_apply(&p2.pos, &p2.pos, &device->transform.view);
        matrix_apply(&p3.pos, &p3.pos, &device->transform.view);
    }
    
    
    z_factor_y = tan(device->fovy*0.5);
    z_factor_x = z_factor_y / device->aspect;
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
    if (p1.pos.z > device->zf)
        vertex_ccodes[0] = CLIP_CODE_GZ;
    else if (p1.pos.z < device->zn)
        vertex_ccodes[0] = CLIP_CODE_LZ;
    else {
        vertex_ccodes[0] = CLIP_CODE_IZ;
        num_verts_in++;
    }

    if (p2.pos.z > device->zf)
        vertex_ccodes[1] = CLIP_CODE_GZ;
    else if (p2.pos.z < device->zn)
        vertex_ccodes[1] = CLIP_CODE_LZ;
    else {
        vertex_ccodes[1] = CLIP_CODE_IZ;
        num_verts_in++;
    }
    
        
    if (p3.pos.z > device->zf)
        vertex_ccodes[2] = CLIP_CODE_GZ;
    else if (p3.pos.z < device->zn)
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
            t1 = (device->zn - p1.pos.z) / v.z;

            xi = p1.pos.x + v.x * t1;
            yi = p1.pos.y + v.y * t1;

            // 用交点覆盖原来的顶点
            p2.pos.x = xi;
            p2.pos.y = yi;
            p2.pos.z = device->zn;

            // 对三角形边v0->v2进行裁剪
            vector_sub(&v, &p3.pos, &p1.pos);
            t2 = (device->zn - p1.pos.z) / v.z;

            xi = p1.pos.x + v.x * t2;
            yi = p1.pos.y + v.y * t2;

            // 用交点覆盖原来的顶点
            p3.pos.x = xi;
            p3.pos.y = yi;
            p3.pos.z = device->zn;

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
            t1 = (device->zn - p1.pos.z) / v.z;

            x01i = p1.pos.x + v.x * t1;
            y01i = p1.pos.y + v.y * t1;

            // 对三角形边v0->v2进行裁剪
            vector_sub(&v, &p3.pos, &p1.pos);
            t2 = (device->zn - p1.pos.z) / v.z;

            x02i = p1.pos.x + v.x * t2;
            y02i = p1.pos.y + v.y * t2;

            // 用交点覆盖原来的顶点
            p1.pos.x = x01i;
            p1.pos.y = y01i;
            p1.pos.z = device->zn;
            
            np2.pos.x = x01i;
            np2.pos.y = y01i;
            np2.pos.z = device->zn;

            np1.pos.x = x02i;
            np1.pos.y = y02i;
            np1.pos.z = device->zn;

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

            device_draw_primitive(device, &np1, &np2, &np3);
            
            cliped = true;
        }
    }

    device_draw_primitive(device, &p1, &p2, &p3);
}


