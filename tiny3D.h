// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!
// modules:
// 1. math library
//  1). vector
typedef struct {
    float x, y, z, w;
} vector_t;
typedef vector_t point_t;

float interp(float a, float b, float t) {
    return b + (a - b) * t;
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
//      7)). normalize
void vector_normalize(vector_t *v) {
    float k = 1 / vector_length(v);
    v->x *= k;
    v->y *= k;
    v->z *= k;
}
//  2). matrix
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
//      5)). apply v = v * m
void matrix_apply(vector_t *v1, const vector_t *v2, const matrix_t *m) {
    for(int j = 0; j < 4; j++)
        v1[j] = v2[0] * m->m[0][j] + v2[1] * m->m[1][j] + v2[2] * m->[2][j] + v2[3] * m->m[3][j];
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
    m->m[0][0] = dx;
    m->m[1][1] = dy;
    m->m[2][2] = dz;
}
//      10)).set_rotate m, x, y, z, theta
// x*x*(1-cos0)+cos0    x*y*(1-cos0)+z*sin0   x*z*(1-cos0)-y*sin0
// x*y*(1-cos0)-z*sin0  y*y*(1-cos0)+cos0     y*z*(1-cos0)+x*sin0
// x*z*(1-cos0)+y*sin0  y*z*(1-cos0)-x*sin0   z*z*(1-cos0)+cos0
void matrix_set_rotate(matrix_t *m, const vector_t *v, float theta) {
    assert(fabs(v*v - 1.0f) < .01f);

    float s = sin(theta);
    float c = cos(theta);

    float a = 1.0f - c;
    float ax = a * v->x;
    float ay = a * v->y;
    float az = a * v->z;

    m->m[0][0] = ax * v->x + c;
    m->m[0][1] = ax * v->y + v->z * s;
    m->m[0][2] = ax * v->z - v->y * s;

    m->m[1][0] = ay * v->x - v->z * s;
    m->m[1][1] = ay * v->y + c;
    m->m[1][2] = ay * v->z + v->x * s;

    m->m[2][0] = az * v->x + v->y * s;
    m->m[2][1] = az * v->y - v->x * s;
    m->m[2][2] = az * v->z + c;

    m->m[3][0] = m->m[3][1] = m->m[3][2] = 0.0f;
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
    float zoomx = zoomy / aspect;
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
    matrix_t transform;
    float w, h;
} transform_t;

//  1). transform_update (world * view * projection
void transform_update(transform_t *ts) {
    matrix_t m;
    matrix_mul(&m, &ts->world, &ts->view);
    matrix_mul(&ts->transform, &m, &ts->projection);
}
//  2). transform_init (ts, width, height)
void transform_init(transform_t *ts, int width, int height) {
    float aspect = (float)width / (float)height;
    matrix_identity(&ts->world);
    matrix_identity(&ts->view);
    matrix_set_perspective(&ts->projection, 3.1415926 * 0.5f, aspect, 1.0f, 500.0f);
    ts->w = (float)width;
    ts->h = (float)height;
    transform_update(ts);
}
//  3). transform_apply
void transform_apply(transform_t *ts, vector_t *y, vector_t *x) {
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
void transform_homogenize(const transform_t *ts, vector *y, const vector *x) {
    float rhw = 1.0f / x->w;
    y->x = (x->x * rhw + 1.0f) * 0.5f * ts->w;
    y->y = (1.0f - x->y * rhw) * 0.5f * ts->h;
    y->z = x->z * rhw;
    y->w = 1.0f;
}

// 3. geometry calculation (vertex, scanline, border check, rect and so on)
typedef struct { float r, g, b; } color_t;
typedef struct { float u, v; } texcoord_t;
typedef struct { point_t pos; texcoord_t tc; color_t color; float rhw; } vertex_t;

typedef struct { vertex_t v, v1, v2; } edge_t;
typedef struct { float top, bottom; edge_t left, right; } trapezoid_t;
typedef struct { vertex_t v, step; int x, y, w; } scanline_t;

void vertex_rhw_init(vertex_t *v) {
    float rhw = 1.0f / v->pos->w;
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


