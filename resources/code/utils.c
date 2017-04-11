//
//  utils.h
//  tinyEngine
//
//  Created by lixuefeng on 2017/3/9.
//  Copyright © 2017年 sdlwlxf. All rights reserved.
//
#include "utils.h"
#include "png.h"
#include "tinyobj_loader_c.h"

#include <math.h>

#define FLT_MIN 1.175494351e-38F
#define FLT_MAX 3.402823466e+38F
#define O_RDONLY         00
#define O_WRONLY         01
#define O_RDWR           02

void CalcNormal(float N[3], float v0[3], float v1[3], float v2[3]) {
    float v10[3];
    float v20[3];
    float len2;
    
    v10[0] = v1[0] - v0[0];
    v10[1] = v1[1] - v0[1];
    v10[2] = v1[2] - v0[2];
    
    v20[0] = v2[0] - v0[0];
    v20[1] = v2[1] - v0[1];
    v20[2] = v2[2] - v0[2];
    
    N[0] = v20[1] * v10[2] - v20[2] * v10[1];
    N[1] = v20[2] * v10[0] - v20[0] * v10[2];
    N[2] = v20[0] * v10[1] - v20[1] * v10[0];
    
    len2 = N[0] * N[0] + N[1] * N[1] + N[2] * N[2];
    if (len2 > 0.0f) {
        float len = (float)sqrt((double)len2);
        
        N[0] /= len;
        N[1] /= len;
    }
}

int make_mesh_and_material_by_obj(vertex_t **mesh, unsigned long *mesh_num, int **material_ids, unsigned long *material_ids_num, const char *name) {
    tinyobj_attrib_t attrib;
    tinyobj_shape_t* shapes = NULL;
    size_t num_shapes;
    tinyobj_material_t* tmaterials = NULL;
    size_t num_materials;
    int start_material_cnt = material_cnt;
    
    FILE * pFile;
    long lSize;
    size_t data_len;
    char * buffer;
    size_t result;
    const char *path = getFilePath(name, "obj");
    pFile = fopen(path, "r");
    fseek(pFile, 0, SEEK_END);
    lSize = ftell (pFile);
    rewind (pFile);
    
    // allocate memory to contain the whole file:
    buffer = (char*) malloc (sizeof(char)*lSize);
    if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}
    
    // copy the file into the buffer:
    result = fread (buffer,1,lSize,pFile);
    if (result > lSize) {fputs ("Reading error",stderr); exit (3);}
    
    
    data_len = (size_t)lSize;
    if (buffer == NULL) {
        exit(-1);
        /* return 0; */
    }
    //printf("filesize: %d\n", (int)data_len);
    
    {
//        printf("path %s\n", path);
        char prefix_path[1000];
        long len = strrchr(path, '/')-path+1;
        strncpy(prefix_path, path, len);
        prefix_path[len] = '\0';
//        printf("prefix_path : %s\n", prefix_path);
        unsigned int flags = TINYOBJ_FLAG_TRIANGULATE;
        int ret = tinyobj_parse_obj(&attrib, &shapes, &num_shapes, &tmaterials,
                                    &num_materials, buffer, data_len, flags, prefix_path);
        if (ret != TINYOBJ_SUCCESS) {
            return 0;
        }
        
        //printf("# of shapes    = %d\n", (int)num_shapes);
        //printf("# of materiasl = %d\n", (int)num_materials);
    }
    
    {
        size_t i;
        for (i = 0; i < num_materials; i++) {
            material_t m;
            tinyobj_material_t tm = tmaterials[i];
            if(tm.name != NULL) {
                m.name = (char*)malloc(sizeof(char) * (strlen(tm.name)+1));
                strcpy(m.name, tm.name);
            }
            m.ambient.r = tm.ambient[0];
            m.ambient.g = tm.ambient[1];
            m.ambient.b = tm.ambient[2];
            
            m.diffuse.r = tm.diffuse[0];
            m.diffuse.g = tm.diffuse[1];
            m.diffuse.b = tm.diffuse[2];
            
            m.specular.r = tm.specular[0];
            m.specular.g = tm.specular[1];
            m.specular.b = tm.specular[2];
            
            m.transmittance.r = tm.transmittance[0];
            m.transmittance.g = tm.transmittance[1];
            m.transmittance.b = tm.transmittance[2];
            
            m.emission.r = tm.emission[0];
            m.emission.g = tm.emission[1];
            m.emission.b = tm.emission[2];
            
            m.shininess = tm.shininess;
            m.ior = tm.ior;
            m.dissolve = tm.dissolve;
            m.illum = tm.illum;
            m.pad0 = tm.pad0;
            
            m.ambient_tex_id = -1;
            m.diffuse_tex_id = -1;
            m.specular_tex_id = -1;
            m.specular_highlight_tex_id = -1;
            m.bump_tex_id = -1;
            m.displacement_tex_id = -1;
            m.alpha_tex_id = -1;
            
            m.ambient_texname = NULL;
            m.diffuse_texname = NULL;
            m.specular_texname = NULL;
            m.specular_highlight_texname = NULL;
            m.bump_texname = NULL;
            m.displacement_texname = NULL;
            m.alpha_texname = NULL;
            
            if(tm.ambient_texname != NULL) {
                m.ambient_texname = (char*)malloc(sizeof(char) * (strlen(tm.ambient_texname)+1));
                strcpy(m.ambient_texname, tm.ambient_texname);
                m.ambient_tex_id = make_texture_by_png(m.ambient_texname, true);
            }
            if(tm.diffuse_texname != NULL) {
                m.diffuse_texname = (char*)malloc(sizeof(char) * (strlen(tm.diffuse_texname)+1));
                strcpy(m.diffuse_texname, tm.diffuse_texname);
                m.diffuse_tex_id = make_texture_by_png(m.diffuse_texname, true);
            }
            if(tm.specular_texname != NULL) {
                m.specular_texname = (char*)malloc(sizeof(char) * (strlen(tm.specular_texname)+1));
                strcpy(m.specular_texname, tm.specular_texname);
                m.specular_tex_id = make_texture_by_png(m.specular_texname, true);
            }
            if(tm.specular_highlight_texname != NULL) {
                m.specular_highlight_texname = (char*)malloc(sizeof(char) * (strlen(tm.specular_highlight_texname)+1));
                strcpy(m.specular_highlight_texname, tm.specular_highlight_texname);
                m.specular_highlight_tex_id = make_texture_by_png(m.specular_highlight_texname, true);
            }
            if(tm.bump_texname != NULL) {
                m.bump_texname = (char*)malloc(sizeof(char) * (strlen(tm.bump_texname)+1));
                strcpy(m.bump_texname, tm.bump_texname);
                m.bump_tex_id = make_texture_by_png(m.bump_texname, true);
            }
            if(tm.displacement_texname != NULL) {
                m.displacement_texname = (char*)malloc(sizeof(char) * (strlen(tm.displacement_texname)+1));
                strcpy(m.displacement_texname, tm.displacement_texname);
                m.displacement_tex_id = make_texture_by_png(m.displacement_texname, true);
            }
            if(tm.alpha_texname != NULL) {
                m.alpha_texname = (char*)malloc(sizeof(char) * (strlen(tm.alpha_texname)+1));
                strcpy(m.alpha_texname, tm.alpha_texname);
                m.alpha_tex_id = make_texture_by_png(m.alpha_texname, true);
            }
            materials[material_cnt++] = m;
        }
    }
    
    {
//        float* vb;
        size_t face_offset = 0;
        size_t i;
        
        /* Assume triangulated face. */
        size_t num_triangles = attrib.num_face_num_verts;
//        size_t stride = 9; /* 9 = pos(3float), normal(3float), color(3float) */
        *mesh_num = num_triangles * 3;
        *mesh = (vertex_t*)malloc(sizeof(vertex_t) * num_triangles * 3);
        
        *material_ids_num = num_triangles;
        *material_ids = (int*)malloc(sizeof(int) * (*material_ids_num));
        for (i = 0; i < attrib.num_face_num_verts; i++) {
            size_t f;
            assert(attrib.face_num_verts[i] % 3 == 0); /* assume all triangle faces. */
            (*material_ids)[i] = attrib.material_ids[i];
            if((*material_ids)[i] >= 0)
                (*material_ids)[i] += start_material_cnt;
            
            for (f = 0; f < (size_t)attrib.face_num_verts[i] / 3; f++) {
                size_t k;
                float v[3][3];
                float n[3][3];
                float c[3];
                float t[3][2];
                float len2;
                
                tinyobj_vertex_index_t idx0 = attrib.faces[face_offset + 3 * f + 0];
                tinyobj_vertex_index_t idx1 = attrib.faces[face_offset + 3 * f + 1];
                tinyobj_vertex_index_t idx2 = attrib.faces[face_offset + 3 * f + 2];
                
                for (k = 0; k < 3; k++) {
                    int f0 = idx0.v_idx;
                    int f1 = idx1.v_idx;
                    int f2 = idx2.v_idx;
                    assert(f0 >= 0);
                    assert(f1 >= 0);
                    assert(f2 >= 0);
                    
                    v[0][k] = attrib.vertices[3 * (size_t)f0 + k];
                    v[1][k] = attrib.vertices[3 * (size_t)f1 + k];
                    v[2][k] = attrib.vertices[3 * (size_t)f2 + k];
                }
                
                if (attrib.num_normals > 0) {
                    int f0 = idx0.vn_idx;
                    int f1 = idx1.vn_idx;
                    int f2 = idx2.vn_idx;
                    if (f0 >= 0 && f1 >= 0 && f2 >= 0) {
                        assert(f0 < (int)attrib.num_normals);
                        assert(f1 < (int)attrib.num_normals);
                        assert(f2 < (int)attrib.num_normals);
                        for (k = 0; k < 3; k++) {
                            n[0][k] = attrib.normals[3 * (size_t)f0 + k];
                            n[1][k] = attrib.normals[3 * (size_t)f1 + k];
                            n[2][k] = attrib.normals[3 * (size_t)f2 + k];
                        }
                    } else { /* normal index is not defined for this face */
                        /* compute geometric normal */
                        CalcNormal(n[0], v[0], v[1], v[2]);
                        n[1][0] = n[0][0];
                        n[1][1] = n[0][1];
                        n[1][2] = n[0][2];
                        n[2][0] = n[0][0];
                        n[2][1] = n[0][1];
                        n[2][2] = n[0][2];
                    }
                } else {
                    /* compute geometric normal */
                    CalcNormal(n[0], v[0], v[1], v[2]);
                    n[1][0] = n[0][0];
                    n[1][1] = n[0][1];
                    n[1][2] = n[0][2];
                    n[2][0] = n[0][0];
                    n[2][1] = n[0][1];
                    n[2][2] = n[0][2];
                }
                
                if (attrib.num_texcoords > 0) {
                    int t0 = idx0.vt_idx;
                    int t1 = idx1.vt_idx;
                    int t2 = idx2.vt_idx;
                    if (t0 >= 0 && t1 >= 0 && t2 >= 0) {
                        assert(t0 < (int)attrib.num_texcoords);
                        assert(t1 < (int)attrib.num_texcoords);
                        assert(t2 < (int)attrib.num_texcoords);
                        for (k = 0; k < 2; k++) {
                            t[0][k] = attrib.texcoords[2 * (size_t)t0 + k];
                            t[1][k] = attrib.texcoords[2 * (size_t)t1 + k];
                            t[2][k] = attrib.texcoords[2 * (size_t)t2 + k];
                        }
                    } else {
                    }
                }
                
                for (k = 0; k < 3; k++) {
                    (*mesh)[3 * i + k].pos = (vector_t){v[k][0], v[k][1], v[k][2], 1};
                    (*mesh)[3 * i + k].normal = (vector_t){n[k][0], n[k][1], n[k][2], 0};
                    
                    (*mesh)[3 * i + k].tc = (texcoord_t){t[k][0], t[k][1]};
                    /* Use normal as color. */
                    c[0] = n[k][0];
                    c[1] = n[k][1];
                    c[2] = n[k][2];
                    len2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
                    if (len2 > 0.0f) {
                        float len = (float)sqrt((double)len2);
                        
                        c[0] /= len;
                        c[1] /= len;
                        c[2] /= len;
                    }
                    (*mesh)[3 * i + k].color = (color_t){(c[0] * 0.5f + 0.5f), (c[1] * 0.5f + 0.5f), (c[2] * 0.5f + 0.5f), 1.0f};
                }
            }
            face_offset += (size_t)attrib.face_num_verts[i];
        }
    }
    
    tinyobj_attrib_free(&attrib);
    tinyobj_shapes_free(shapes, num_shapes);
    tinyobj_materials_free(tmaterials, num_materials);
    
    // terminate
    fclose (pFile);
    free (buffer);
    return 0;
}

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
int load_png_image( const char *name, unsigned int **bits, unsigned int *width, unsigned int *height)
{
    FILE *fp;
    png_structp png_ptr;
    png_infop info_ptr;
    png_bytep* row_pointers;
    char buf[PNG_BYTES_TO_CHECK];
    int w, h, x, y, temp, color_type;
    
    fp = fopen(getFilePath(name, "png"), "rb");
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
    
    *bits = (unsigned int*)malloc(sizeof(unsigned int) * w * h);
    
    /* 获取图像的所有行像素数据，row_pointers里边就是rgba数据 */
    row_pointers = png_get_rows( png_ptr, info_ptr );
    /* 根据不同的色彩类型进行相应处理 */
    switch( color_type ) {
        case PNG_COLOR_TYPE_RGB_ALPHA:
            for( y=0; y<h; ++y ) {
                for( x=0; x < w; ++x ) {
                    (*bits)[y*w+x] = 0;
                    /* 以下是RGBA数据，需要自己补充代码，保存RGBA数据 */
                    (*bits)[y*w+x] |= row_pointers[y][4*x+0] << 16; // red
                    (*bits)[y*w+x] |= row_pointers[y][4*x+1] << 8; // green
                    (*bits)[y*w+x] |= row_pointers[y][4*x+2]; // blue
                    (*bits)[y*w+x] |= row_pointers[y][4*x+3] << 24; // alpha
                }
            }
            break;
            
        case PNG_COLOR_TYPE_RGB:
            for( y=0; y<h; ++y ) {
                for( x=0; x<w; ++x ) {
                    (*bits)[y*w+x] = 0xff000000;
                    (*bits)[y*w+x] |= row_pointers[y][3*x+0] << 16; // red
                    (*bits)[y*w+x] |= row_pointers[y][3*x+1] << 8; // green
                    (*bits)[y*w+x] |= row_pointers[y][3*x+2]; // blue
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
    
    *width = w;
    *height = h;
    
    return 0;
}

// 返回-1读取失败，返回id号对应texture编号
int make_texture_by_png(const char *name, bool mipmap) {
    IUINT32 *data = NULL;
    texture_t *texture = &textures[texture_count];
    char trueName[100];
    char *findPos = NULL;
    if((findPos = strrchr(name, '.')) != NULL) {
        long len = findPos-name;
        strncpy(trueName, name, len);
        trueName[len] = '\0';
    }
    int res = load_png_image(findPos == NULL ? name : trueName, &data, &texture->width, &texture->height);
    if(res == 0) {
        texture->datas = (IUINT32**)malloc(1 * sizeof(IUINT32*));
        texture->datas[0] = data;
        if(mipmap) {
            texture->use_mipmap = true;
            generate_mipmaps(texture, 1.01);
        }
        texture_count++;
        return texture_count-1;
    }
    return -1;
}

