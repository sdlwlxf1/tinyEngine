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
    float k = 1 / vector_length(v);
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
}

//      6)). transpose
void matrix_transpose(matrix_t *m) {
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            m->m[i][j] = m->m[j][i];
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
    matrix_set_identity(&ts->world);
    matrix_set_identity(&ts->view);
    matrix_set_perspective(&ts->projection, 3.1415926 * 0.5f, aspect, 1.0f, 500.0f);
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
    y->x = (x->x * rhw + 1.0f) * 0.5f * ts->w;
    y->y = (1.0f - x->y * rhw) * 0.5f * ts->h;
    y->z = x->z * rhw;
    y->w = 1.0f;
}

// 3. geometry calculation (vertex, scanline, border check, rect and so on)
typedef struct { float r, g, b; } color_t;
void color_init(color_t *c) {
    c->r = 0.0f;
    c->g = 0.0f;
    c->b = 0.0f;
}
//      4)). product
void color_product(color_t *c, const color_t *a, const color_t *b) {
    c->r = a->r * b->r;
    c->g = a->g * b->g;
    c->b = a->b * b->b;
}
void color_scale(color_t *c, float k) {
    c->r *= k;
    c->g *= k;
    c->b *= k;
}
void color_add(color_t *c, const color_t *a, const color_t *b) {
    c->r = a->r + b->r;
    c->g = a->g + b->g;
    c->b = a->b + b->b;
}

point_t viewPos;

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
material_t material;

typedef struct {
    point_t pos;
    float constant;
    float linear;
    float quadratic;
    color_t ambi;
    color_t diff;
    color_t spec;
} pointlight_t;
#define NR_POINT_LIGHTS 4
pointlight_t pointLights[NR_POINT_LIGHTS];

void calc_pointlight(color_t *color, const pointlight_t *light, const vector_t *normal, const vector_t *fpos, const vector_t *viewdir) {
    vector_t lightDir;
    vector_sub(&lightDir, &light->pos, fpos);
    float distance = vector_length(&lightDir);
    vector_normalize(&lightDir);
    float diff = fmaxf(vector_dotproduct(normal, &lightDir), 0.0f);
    vector_t vec;
    vector_inverse(&lightDir);
    vector_reflect(&vec, &lightDir, normal);
    float spec = powf(fmaxf(vector_dotproduct(viewdir, &vec), 0.0f), material.shininess);
    float attenuation = 1.0f / (light->constant + light->linear * distance + light->quadratic * (distance * distance));
    
    color_t c;
    color_product(&c, &light->ambi, &material.ambi);
    color_scale(&c, attenuation);
    color_add(color, color, &c);
    color_product(&c, &light->diff, &material.diff);
    color_scale(&c, diff * attenuation);
    color_add(color, color, &c);
    color_product(&c, &light->spec, &material.spec);
    color_scale(&c, spec * attenuation);
    color_add(color, color, &c);
}

typedef struct {
    vector_t dir;
    color_t ambi;
    color_t diff;
    color_t spec;
} dirlight_t;
dirlight_t dirLight;

void calc_dirlight(color_t *color, const dirlight_t *light, const vector_t *normal, const vector_t *viewdir) {
    vector_t lightDir = light->dir;
    vector_inverse(&lightDir);
    vector_normalize(&lightDir);
    float diff = fmaxf(vector_dotproduct(normal, &lightDir), 0.0f);
    vector_t vec;
    vector_inverse(&lightDir);
    vector_reflect(&vec, &lightDir, normal);
    float spec = powf(fmaxf(vector_dotproduct(viewdir, &vec), 0.0f), material.shininess);
    
    color_t temp;
    color_product(&temp, &light->ambi, &material.ambi);
    color_add(color, color, &temp);
    color_product(&temp, &light->diff, &material.diff);
    color_scale(&temp, diff);
    color_add(color, color, &temp);
    color_product(&temp, &light->spec, &material.spec);
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
	IUINT32 **framebuffer;      // 像素缓存：framebuffer[y] 代表第 y行
	float **zbuffer;            // 深度缓存：zbuffer[y] 为第 y行指针
	IUINT32 **texture;          // 纹理：同样是每行索引
	int tex_width;              // 纹理宽度
	int tex_height;             // 纹理高度
	float max_u;                // 纹理最大宽度：tex_width - 1
	float max_v;                // 纹理最大高度：tex_height - 1
	int render_state;           // 渲染状态
	IUINT32 background;         // 背景颜色
	IUINT32 foreground;         // 线框颜色
}	device_t;

#define RENDER_STATE_WIREFRAME      1		// 渲染线框
#define RENDER_STATE_TEXTURE        2		// 渲染纹理
#define RENDER_STATE_COLOR          4		// 渲染颜色

// 设备初始化，fb为外部帧缓存，非 NULL 将引用外部帧缓存（每行 4字节对齐）
void device_init(device_t *device, int width, int height, void *fb) {
	int need = sizeof(void*) * (height * 2 + 1024) + width * height * 8;
	char *ptr = (char*)malloc(need + 64);
	char *framebuf, *zbuf;
	int j;
	device->framebuffer = (IUINT32**)ptr;
    ptr += sizeof(void*) * height;
    
	device->zbuffer = (float**)ptr;
	ptr += sizeof(void*) * height;

	device->texture = (IUINT32**)ptr;
	ptr += sizeof(void*) * 1024;

	framebuf = (char*)ptr;
	zbuf = (char*)ptr + width * height * 4;
	ptr += width * height * 8;
	if (fb != NULL) framebuf = (char*)fb;
	for (j = 0; j < height; j++) {
		device->framebuffer[j] = (IUINT32*)(framebuf + width * 4 * j);
		device->zbuffer[j] = (float*)(zbuf + width * 4 * j);
	}
	device->texture[0] = (IUINT32*)ptr;
	device->texture[1] = (IUINT32*)(ptr + 16);
	memset(device->texture[0], 0, 64);
	device->tex_width = 2;
	device->tex_height = 2;
	device->max_u = 1.0f;
	device->max_v = 1.0f;
	device->width = width;
	device->height = height;
	device->background = 0xc0c0c0;
	device->foreground = 0;
	transform_init(&device->transform, width, height);
	device->render_state = RENDER_STATE_WIREFRAME;
}

void device_destroy(device_t *device) {
	if (device->framebuffer) 
		free(device->framebuffer);
	device->framebuffer = NULL;
	device->zbuffer = NULL;
	device->texture = NULL;
}

void device_set_texture(device_t *device, void *bits, long pitch, int w, int h) {
	char *ptr = (char*)bits;
	int j;
	for (j = 0; j < h; ptr += pitch, j++) 	// 重新计算每行纹理的指针
		device->texture[j] = (IUINT32*)ptr;
	device->tex_width = w;
	device->tex_height = h;
	device->max_u = (float)(w - 1);
	device->max_v = (float)(h - 1);
}

void device_clear(device_t *device, int mode) {
	int y, x, height = device->height;
	for (y = 0; y < device->height; y++) {
		IUINT32 *dst = device->framebuffer[y];
		IUINT32 cc = (height - 1 - y) * 230 / (height - 1);
		cc = (cc << 16) | (cc << 8) | cc;
		if (mode == 0) cc = device->background;
		for (x = device->width; x > 0; dst++, x--) dst[0] = cc;
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

color_t device_texture_read(const device_t *device, float u, float v) {
    color_t color;
    int x, y;
    u = u * device->max_u;
    v = v * device->max_v;
    x = (int)(u + 0.5f);
    y = (int)(v + 0.5f);
    x = CMID(x, 0, device->tex_width - 1);
    y = CMID(y, 0, device->tex_height - 1);
    IUINT32 c = device->texture[y][x];
    color.r = c & 0xff;
    color.g = (c >> 8) & 0xff;
    color.b = (c >> 16) & 0xff;
    return color;
}

bool computeBarycentricCoords3d(point_t *res, const point_t *v, const point_t *p) {
    vector_t d1, d2, n;
    vector_sub(&d1, &v[1], &v[0]);
    vector_sub(&d2, &v[2], &v[1]);
    vector_crossproduct(&n, &d1, &d2);
    float u1, u2, u3, u4;
    float v1, v2, v3, v4;
    if((fabs(n.x) >= fabs(n.y)) && (fabs(n.x) >= fabs(n.z))) {
        u1 = v[0].y - v[2].y;
        u2 = v[1].y - v[2].y;
        u3 = p->y - v[0].y;
        u4 = p->y - v[2].y;
        v1 = v[0].z - v[2].z;
        v2 = v[1].z - v[2].z;
        v3 = p->z - v[0].z;
        v4 = p->z - v[2].z;
    } else if(fabs(n.y) >= fabs(n.z)) {
        u1 = v[0].z - v[2].z;
        u2 = v[1].z - v[2].z;
        u3 = p->z - v[0].z;
        u4 = p->z - v[2].z;
        v1 = v[0].x - v[2].x;
        v2 = v[1].x - v[2].x;
        v3 = p->x - v[0].x;
        v4 = p->x - v[2].x;
    } else {
        u1 = v[0].x - v[2].x;
        u2 = v[1].x - v[2].x;
        u3 = p->x - v[0].x;
        u4 = p->x - v[2].x;
        v1 = v[0].y - v[2].y;
        v2 = v[1].y - v[2].y;
        v3 = p->y - v[0].y;
        v4 = p->y - v[2].y;
    }
    
    float denom = v1 * u2 - v2 * u1;
    if (denom == 0.0f) {
        return false;
    }
    float oneOverDenom = 1.0f / denom;
    res->x = (v4 * u2 - v2 * u4) * oneOverDenom;
    res->y = (v1 * u3 - v3 * u1) * oneOverDenom;
    res->z = 1.0f - res->x - res->y;
    return true;
}

// now is the core realize for rendering, so just do it.
void device_draw_scanline(device_t *device, scanline_t *scanline, const point_t *t, const vertex_t *vt, const vector_t *n) {
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
                color_t color = {0.0f, 0.0f, 0.0f};
                
                point_t barycenter = {0.0f, 0.0f, 0.0f, 1.0f};
                computeBarycentricCoords3d(&barycenter, t, &scanline->v.pos);
                
                point_t fragPos = {0.0f, 0.0f, 0.0f, 1.0f};
                point_t ptemp = vt[0].pos;
                vector_scale(&ptemp, barycenter.x);
                vector_add(&fragPos, &fragPos, &ptemp);
                ptemp = vt[1].pos;
                vector_scale(&ptemp, barycenter.y);
                vector_add(&fragPos, &fragPos, &ptemp);
                ptemp = vt[2].pos;
                vector_scale(&ptemp, barycenter.z);
                vector_add(&fragPos, &fragPos, &ptemp);
                
                vector_t normal = {0.0f, 0.0f, 0.0f, 0.0f};
                ptemp = n[0];
                vector_scale(&ptemp, barycenter.x);
                vector_add(&normal, &normal, &ptemp);
                ptemp = n[1];
                vector_scale(&ptemp, barycenter.y);
                vector_add(&normal, &normal, &ptemp);
                ptemp = n[2];
                vector_scale(&ptemp, barycenter.z);
                vector_add(&normal, &normal, &ptemp);
                vector_normalize(&normal);
                //printf("%f, %f, %f\n", normal.x, normal.y, normal.z);
                vector_t viewdir;
                vector_sub(&viewdir, &viewPos, &fragPos);
                vector_normalize(&viewdir);
                calc_dirlight(&color, &dirLight, &normal, &viewdir);
                color_t temp = {0.0f, 0.0f ,0.0f};
                int i = 0;
                for(i = 0; i < NR_POINT_LIGHTS; i++)
                {
                    calc_pointlight(&temp, &pointLights[i], &normal, &fragPos, &viewdir);
                    color_add(&color, &color, &temp);
                }
                
                float r = 0.0f;
                float g = 0.0f;
                float b = 0.0f;
                
                if(render_state & RENDER_STATE_COLOR) {
                    r = scanline->v.color.r * w * 255.0f;
                    g = scanline->v.color.g * w * 255.0f;
                    b = scanline->v.color.b * w * 255.0f;
                }
                if(render_state & RENDER_STATE_TEXTURE) {
                    float u = scanline->v.tc.u * w;
                    float v = scanline->v.tc.v * w;
                    color_t cc = device_texture_read(device, u, v);
                    r = cc.r;
                    g = cc.g;
                    b = cc.b;
                }
                r = color.r * 255.0f;
                g = color.g * 255.0f;
                b = color.b * 255.0f;
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
void device_render_trap(device_t *device, trapezoid_t *trap, const point_t *p, const vertex_t *t, const vector_t *n) {
    scanline_t scanline;
    int j, top, bottom;
    top = (int)(trap->top + 0.5f);
    bottom = (int)(trap->bottom + 0.5f);
    for(j = top; j < bottom; j++) {
        if(j >= 0 && j < device->height) {
            trapezoid_edge_interp(trap, (float)j + 0.5f);
            trapezoid_init_scan_line(trap, &scanline, j);
            device_draw_scanline(device, &scanline, p, t, n);
        }
        if(j >= device->height)
            break;
    }
}

void device_draw_primitive(device_t *device, vertex_t *v1,
    vertex_t *v2, vertex_t *v3) {
    point_t p1, p2, p3, c1, c2, c3;
    int render_state = device->render_state;
    
    // 使用法向量
//    matrix_apply(&c1, &v1->pos, &device->transform.world);
//    matrix_apply(&c2, &v2->pos, &device->transform.world);
//    matrix_apply(&c3, &v3->pos, &device->transform.world);
    
//    matrix_apply(&p1, &c1, &device->transform.view);
//    matrix_apply(&p2, &c2, &device->transform.view);
//    matrix_apply(&p3, &c3, &device->transform.view);
    
//    vector_t v21, v32;
//    vector_sub(&v21, &p2, &p1);
//    vector_sub(&v32, &p3, &p2);
//    vector_t normal;
//    vector_crossproduct(&normal, &v21, &v32);
//    if(vector_dotproduct(&normal, &p2) >= 0)
//        return;

    transform_apply(&device->transform, &c1, &v1->pos);
    transform_apply(&device->transform, &c2, &v2->pos);
    transform_apply(&device->transform, &c3, &v3->pos);
    
    // 法向量乘正规矩阵
    vector_t normals[3];
    matrix_t nm;
    matrix_clone(&nm, &device->transform.world);
    matrix_inverse(&nm);
    matrix_transpose(&nm);
    matrix_apply(&normals[0], &v1->normal, &nm);
    matrix_apply(&normals[1], &v2->normal, &nm);
    matrix_apply(&normals[2], &v3->normal, &nm);
    vector_normalize(&normals[0]);
    vector_normalize(&normals[1]);
    vector_normalize(&normals[2]);

    if (transform_check_cvv(&c1) != 0) return;
    if (transform_check_cvv(&c2) != 0) return;
    if (transform_check_cvv(&c3) != 0) return;
    
    transform_homogenize(&device->transform, &p1, &c1);
    transform_homogenize(&device->transform, &p2, &c2);
    transform_homogenize(&device->transform, &p3, &c3);
    
    // 背面剔除
    vector_t v21, v32;
    vector_sub(&v21, &p2, &p1);
    vector_sub(&v32, &p3, &p2);
    // 计算叉积
    if(v21.x * v32.y - v32.x * v21.y <= 0)
        return;
    
    if (render_state & (RENDER_STATE_TEXTURE | RENDER_STATE_COLOR)) {
        vertex_t t[3];
        point_t p[3];
        t[0] = *v1, t[1] = *v2, t[2] = *v3;
        trapezoid_t traps[2];
        int n;

        p[0] = p1;
        p[1] = p2;
        p[2] = p3;
        t[0].pos = p1;
        t[1].pos = p2;
        t[2].pos = p3;
        t[0].pos.w = c1.w;
        t[1].pos.w = c2.w;
        t[2].pos.w = c3.w;

        vertex_rhw_init(&t[0]);
        vertex_rhw_init(&t[1]);
        vertex_rhw_init(&t[2]);

        n = trapezoid_init_triangle(traps, &t[0], &t[1], &t[2]);
        
        matrix_apply(&t[0].pos, &v1->pos, &device->transform.world);
        matrix_apply(&t[1].pos, &v2->pos, &device->transform.world);
        matrix_apply(&t[2].pos, &v3->pos, &device->transform.world);

        if(n >= 1) device_render_trap(device, &traps[0], p, t, normals);
        if(n >= 2) device_render_trap(device, &traps[1], p, t, normals);
    }

    if(render_state & RENDER_STATE_WIREFRAME) {
        device_draw_line(device, (int)p1.x, (int)p1.y, (int)p2.x, (int)p2.y, device->foreground);
        device_draw_line(device, (int)p1.x, (int)p1.y, (int)p3.x, (int)p3.y, device->foreground);
        device_draw_line(device, (int)p3.x, (int)p3.y, (int)p2.x, (int)p2.y, device->foreground);
    }
}
